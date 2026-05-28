function results = slowWaveAnalysis_new(signal, lowPassOn, lowPassCutoff, lowPassOrder, fs, ...
    window, figtrue, folderpath, condition, blankIdx, edgeBufferSec)
% slowWaveAnalysis
%
% Extract slow waves from stomach EMG data. Operates on a single signal
% matrix; callers are expected to slice the signal themselves before
% calling (the previous 'stim_rec' struct mode was removed).
%
% INPUTS
%   signal        : [N x 3] stomach EMG matrix. Pass only the stomach
%                   channels (the function uses all columns as candidate
%                   slow-wave channels).
%   lowPassOn     : logical, apply low-pass filter
%   lowPassCutoff : scalar, low-pass cutoff frequency (Hz). Typical: 2 Hz
%   lowPassOrder  : scalar, low-pass filter order. Typical: 4
%                   Filter settling time = lowPassOrder/lowPassCutoff seconds.
%                   edgeBufferSec must exceed this value.
%   fs            : scalar, sampling frequency (Hz). Expected: 25000 Hz.
%   window        : scalar, gaussian smoothing window in seconds (default: 5s).
%                   5s preserves slow waves up to ~8 cpm; increase to attenuate more.
%   figtrue       : logical, generate and save diagnostic figures (default: false)
%   folderpath    : folder to save outputs and figures
%   condition     : base condition name for saving
%   mode          : 'single' or 'stim_rec'
%   blankIdx      : - single mode: [nSeg x 2] blanked sample index ranges
%                   - stim_recovery mode: struct with fields .stim and .recovery
%   edgeBufferSec : seconds to mask around blank boundaries and signal edges.
%                   Default = 3s (covers order=4, cutoff=2Hz settling of 2s + margin).
%
% OUTPUT
%   results : struct with slow wave time series, peak locations, rate series, QC fields
%
% SLOW WAVE RATE COMPUTATION (per channel, channel 1 plotted in figures)
%   Centered 60s window (rateWinSec, hardcoded) stepped every 1s.
%   For each timepoint t = slowWaveRateTime(ti):
%     - Skip (NaN) when the centered window cannot fit fully, i.e.
%       t < 30s or t > T - 30s. (Symmetrical to HR/BR script.)
%     - Skip (NaN) when the sample at t is inside invalidMask
%       (don't report a rate at a blanked timepoint).
%     - Otherwise, find longest contiguous valid+non-edge stretch in the
%       window. If stretch >= 30s AND >= 3 detected peaks in the stretch:
%       rate = nPeaks / (stretchDuration / 60) [cpm]. Else NaN.
%   Per-window inter-peak plausibility gate intentionally removed (matches
%   HR/BR pattern); rates > 8 cpm are flagged globally via
%   sw_implausibleFraction but not discarded.
%
% CHANNELS
%   Only channels 1-3 (stomach) are analyzed. Channels 4-5 (nerve) are ignored.
%
% DIAGNOSTIC FIGURES (channel 1 only, saved as .fig and .png):
%   Figure: Slow wave overview (3 panels)
%     1. Raw signal (channel 1) with invalid (gray) and edge-buffer (yellow) regions shaded
%     2. Lowpass + Gaussian-smoothed signal (channel 1) with detected
%        slow wave peaks marked.
%     3. Slow wave rate time series (channel 1)

    if nargin < 10 || isempty(blankIdx),      blankIdx      = []; end
    if nargin < 11 || isempty(edgeBufferSec), edgeBufferSec = 3.0; end

    results = analyzeOneSignal(signal, lowPassOn, lowPassCutoff, lowPassOrder, ...
        fs, window, figtrue, folderpath, condition, blankIdx, edgeBufferSec);
end


% =========================================================================
function out = analyzeOneSignal(signal, lowPassOn, lowPassCutoff, lowPassOrder, ...
    fs, window, figtrue, folderpath, condition, blankIdx, edgeBufferSec)

    % --- Input validation ---
    assert(fs >= 100 && fs <= 30000, ...
        'fs = %.1f Hz seems implausible. Check units.', fs);

    settlingTimeSec = lowPassOrder / lowPassCutoff;
    if edgeBufferSec < settlingTimeSec
        warning('slowWaveAnalysis: edgeBufferSec (%.2fs) < filter settling time (%.2fs). Consider increasing it.', ...
            edgeBufferSec, settlingTimeSec);
    end

    condition      = string(condition);
    condition      = condition(1);
    conditionLabel = char(condition);
    folderpath     = char(string(folderpath));

    validateattributes(signal, {'numeric'}, {'2d','nonempty'});
    N          = size(signal, 1);
    swChans    = 1:3;       % Stomach channels only
    plotChan   = 1;         % Channel shown in diagnostic figures
    t          = (0:N-1)' / fs;

    if isempty(blankIdx)
        blankIdx = zeros(0,2);
    end

    % ---------------- Invalid / edge masks ----------------
    % invalidMask is joint across all channels (artifacts are correlated)
    invalidMask = any(isnan(signal), 2);

    if ~isempty(blankIdx)
        blankIdx      = round(blankIdx);
        blankIdx(:,1) = max(blankIdx(:,1), 1);
        blankIdx(:,2) = min(blankIdx(:,2), N);
        for k = 1:size(blankIdx,1)
            invalidMask(blankIdx(k,1):blankIdx(k,2)) = true;
        end
    end

    edgeBufferSamp = round(edgeBufferSec * fs);

    % Edge mask: blank boundaries + signal edges
    edgeMask = false(N,1);
    for k = 1:size(blankIdx,1)
        s1 = max(1, blankIdx(k,1) - edgeBufferSamp);
        s2 = min(N, blankIdx(k,2) + edgeBufferSamp);
        edgeMask(s1:s2) = true;
    end
    edgeMask(1 : min(N, edgeBufferSamp))         = true;
    edgeMask(max(1, N - edgeBufferSamp + 1) : N) = true;

    % Pure edge-buffer region (not invalid) — for differential shading
    edgeOnlyMask = edgeMask & ~invalidMask;

    % ---------------- Filtering ----------------
    xFill = signal;
    for ch = 1:size(signal,2)
        xFill(:,ch) = fillmissing(xFill(:,ch), 'linear', 'EndValues', 'nearest');
    end

    if lowPassOn
        Wn       = lowPassCutoff / (fs/2);
        [z,p,k_] = butter(lowPassOrder, Wn, 'low');
        [sos,g]  = zp2sos(z, p, k_);
        filtSignal = filtfilt(sos, g, xFill);
    else
        filtSignal = xFill;
    end

    % Detrend always — slow wave analysis on drifting baseline produces spurious peaks
    filteredSignal = detrend(filtSignal);

    % ---------------- Smoothing ----------------
    % Gaussian smoothing of width `window` (seconds) is applied after the
    % lowpass+detrend stage. The smoothing suppresses residual transition-band
    % ripple that the order-2 lowpass leaves between true slow-wave peaks
    % (and that otherwise causes findpeaks to register doublets).
    windowlen = max(1, floor(window * fs));
    slowWaveTimeSeries = smoothdata(filteredSignal, 1, 'gaussian', windowlen);

    % Restore NaNs after smoothing for display/bookkeeping
    slowWaveTimeSeries(invalidMask, :) = NaN;

    % ---------------- Peak detection parameters ----------------
    % MinPeakDistance: equivalent to 10 cpm (6s period) — prevents within-cycle
    % double detection without preventing genuine slow wave absences
    minPeakDist_samp    = round(6 * fs);
    maxSlowWaveRate_cpm = 8;

    % Slow wave rate window parameters
    rateWinSec          = 60;
    minStretchSec       = 30;
    minPeaksInStretch   = 3;

    % Rate time axis: 1-second steps
    rateStepSamp = round(fs);
    rateT_idx    = 1 : rateStepSamp : N;
    rateT_idx    = rateT_idx(rateT_idx >= 1 & rateT_idx <= N);
    nRateT       = numel(rateT_idx);
    slowWaveRateTime = t(rateT_idx);

    nSwChan = numel(swChans);
    slowWaveRateSeries     = NaN(nRateT, nSwChan);
    avgSlowWave            = NaN(1, nSwChan);
    slowWavePeakLocs       = cell(1, nSwChan);
    sw_implausibleFraction = NaN(1, nSwChan);

    % Combined clean mask: valid AND not in edge buffer
    cleanMask        = ~invalidMask & ~edgeMask;
    minStretchSamp   = round(minStretchSec * fs);
    sigDurSec        = (N - 1) / fs;

    % ---------------- Peak detection and rate per stomach channel ----------------
    for ci = 1:nSwChan
        i = swChans(ci);

        % Smooth the lowpass+detrended signal for peak detection
        % (matches the stored slowWaveTimeSeries built above).
        yForPeaks = smoothdata(filteredSignal(:,i), 1, 'gaussian', windowlen);

        [~, locsRaw] = findpeaks(yForPeaks, ...
            'MinPeakDistance', minPeakDist_samp);

        % Reject peaks in invalid or edge regions
        validPeakMask = ~invalidMask(locsRaw) & ~edgeMask(locsRaw);
        locs = locsRaw(validPeakMask);
        slowWavePeakLocs{ci} = locs;

        % Flag implausible inter-peak rates (do not discard)
        if numel(locs) >= 2
            peakIntervals_sec = diff(locs) / fs;
            peakRates_cpm     = 60 ./ peakIntervals_sec;
            implausible       = peakRates_cpm > maxSlowWaveRate_cpm;
            sw_implausibleFraction(ci) = mean(implausible);
            if sw_implausibleFraction(ci) > 0.05
                warning('%s ch%d: %.1f%% of slow wave intervals exceed %.0f cpm — check detection.', ...
                    conditionLabel, i, sw_implausibleFraction(ci)*100, maxSlowWaveRate_cpm);
            end
        end

        % Moving rate: longest continuous clean stretch in each centered window.
        % Behavior consistent with HR_BR_HRVAnalysis_new.m:
        %   - Centered window of size rateWinSec around metrics_t(ti)
        %   - NaN where the window cannot fit fully in the signal
        %     (edges: t < rateWinSec/2 or t > T - rateWinSec/2)
        %   - NaN where the center sample is itself inside invalidMask
        %     (don't report a rate at a blanked timepoint).
        % Per-window inter-peak plausibility gate intentionally removed.
        % Global sw_implausibleFraction (computed above) still reports the
        % fraction of intervals exceeding maxSlowWaveRate_cpm for QC.
        halfWinSamp = round(rateWinSec * fs / 2);
        halfWinSec  = rateWinSec / 2;

        for ti = 1:nRateT
            ctrSamp = rateT_idx(ti);
            tc      = slowWaveRateTime(ti);

            % Edge guard: skip if centered window does not fit fully.
            if tc < halfWinSec || tc > sigDurSec - halfWinSec
                continue;
            end

            % Skip windows whose center is inside a blanked region.
            if invalidMask(ctrSamp)
                continue;
            end

            winStartSamp = max(1, ctrSamp - halfWinSamp);
            winEndSamp   = min(N, ctrSamp + halfWinSamp);

            winClean = cleanMask(winStartSamp:winEndSamp);
            [stretchLen, stretchStart] = longestContiguousRun(winClean);

            if stretchLen < minStretchSamp
                continue;
            end

            absStart = winStartSamp + stretchStart - 1;
            absEnd   = min(N, absStart + stretchLen - 1);

            peaksInStretch = locs(locs >= absStart & locs <= absEnd);

            if numel(peaksInStretch) >= minPeaksInStretch
                stretchDurSec = stretchLen / fs;
                slowWaveRateSeries(ti, ci) = numel(peaksInStretch) / (stretchDurSec / 60);
            end
        end

        avgSlowWave(ci) = mean(slowWaveRateSeries(:,ci), 'omitnan');
        fprintf('%s average slow wave rate ch%d: %.2f cpm\n', conditionLabel, i, avgSlowWave(ci));
    end

    % ---------------- Diagnostic figure (channel 1 only) ----------------
    if figtrue
        invalidColor = [0.6 0.6 0.6];   % gray   — blanked/invalid
        edgeColor    = [1.0 0.95 0.6];  % yellow — edge buffer only

        fig = figure('Name', sprintf('%s — Slow waves ch%d', conditionLabel, plotChan), ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.75]);

        ax = gobjects(3,1);
        for p = 1:3
            ax(p) = subplot(3,1,p, 'Parent', fig);
        end

        % Panel 1: Raw signal channel 1
        rawCh = signal(:, plotChan);
        plot(ax(1), t, rawCh, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
        ylabel(ax(1), 'Raw (a.u.)');
        title(ax(1), sprintf('%s — Raw signal (ch%d)', conditionLabel, plotChan));
        applyShading(ax(1), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);
        addShadingLegend(ax(1), invalidColor, edgeColor);

        % Panel 2: Smoothed signal + peaks for channel 1
        locs1 = slowWavePeakLocs{plotChan};
        plot(ax(2), t, slowWaveTimeSeries(:, plotChan), ...
            'Color', [0.15 0.5 0.25], 'LineWidth', 0.8);
        hold(ax(2), 'on');
        if ~isempty(locs1)
            plot(ax(2), t(locs1), slowWaveTimeSeries(locs1, plotChan), ...
                'v', 'Color', [0.8 0.2 0.1], 'MarkerFaceColor', [0.8 0.2 0.1], ...
                'MarkerSize', 6, 'LineStyle', 'none');
        end
        ylabel(ax(2), 'LP+smooth (a.u.)');
        title(ax(2), sprintf('Lowpass + Gaussian-smoothed signal + peaks  (detected: %d)', numel(locs1)));
        applyShading(ax(2), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 3: Slow wave rate series for channel 1
        plot(ax(3), slowWaveRateTime, slowWaveRateSeries(:, plotChan), ...
            'Color', [0.4 0.2 0.7], 'LineWidth', 1);
        ylabel(ax(3), 'Rate (cpm)');
        xlabel(ax(3), 'Time (s)');
        title(ax(3), sprintf('Slow wave rate ch%d  (mean: %.2f cpm)', ...
            plotChan, avgSlowWave(plotChan)));
        applyShading(ax(3), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        linkaxes(ax, 'x');
        xlim(ax(1), [t(1) t(end)]);

        saveFigure(fig, folderpath, sprintf('%s_slowwave', conditionLabel));
    end

    % ---------------- Save ----------------
    swFile = fullfile(folderpath, sprintf('%s_slowWaves.mat', conditionLabel));
    save(swFile, ...
        'slowWaveTimeSeries', 'slowWaveRateSeries', 'slowWaveRateTime', ...
        'avgSlowWave', 'slowWavePeakLocs', 'sw_implausibleFraction', ...
        'invalidMask', 'edgeMask', 'blankIdx', 'edgeBufferSec', ...
        'rateWinSec', 'minStretchSec', 'minPeaksInStretch', ...
        't', 'fs', 'window', 'windowlen');

    % ---------------- Output struct ----------------
    out                        = struct();
    out.slowWaveTimeSeries     = slowWaveTimeSeries;
    out.slowWaveRateSeries     = slowWaveRateSeries;
    out.slowWaveRateTime       = slowWaveRateTime;
    out.avgSlowWave            = avgSlowWave;
    out.slowWavePeakLocs       = slowWavePeakLocs;
    out.sw_implausibleFraction = sw_implausibleFraction;
    out.invalidMask            = invalidMask;
    out.edgeMask               = edgeMask;
    out.edgeOnlyMask           = edgeOnlyMask;
    out.blankIdx               = blankIdx;
    out.t                      = t;
    out.fs                     = fs;
    out.rateWinSec             = rateWinSec;
    out.minStretchSec          = minStretchSec;
    out.minPeaksInStretch      = minPeaksInStretch;
    out.window                 = window;
    out.windowlen              = windowlen;
end


% =========================================================================
function [runLen, runStart] = longestContiguousRun(binaryVec)
% Find length and 1-based start index of the longest contiguous run of true values.
    if ~any(binaryVec)
        runLen = 0; runStart = 0; return;
    end
    d       = diff([0; binaryVec(:); 0]);
    starts  = find(d ==  1);
    ends    = find(d == -1) - 1;
    lengths = ends - starts + 1;
    [runLen, idx] = max(lengths);
    runStart = starts(idx);
end


% =========================================================================
function applyShading(ax, t, invalidMask, edgeOnlyMask, invalidColor, edgeColor)
% Shade invalid regions (gray) and edge-buffer-only regions (yellow) on ax.
    shadeRegions(ax, t, edgeOnlyMask, edgeColor);
    shadeRegions(ax, t, invalidMask,  invalidColor);
end


% =========================================================================
function shadeRegions(ax, t, mask, color)
% Shade contiguous true regions of mask on axes ax using patch.
    if ~any(mask)
        return;
    end
    yl = ylim(ax);
    if yl(1) == yl(2)
        yl = [-1 1];
    end
    d      = diff([0; mask(:); 0]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;
    starts = min(starts, numel(t));
    ends   = min(ends,   numel(t));
    hold(ax, 'on');
    for k = 1:numel(starts)
        x0 = t(starts(k));
        x1 = t(ends(k));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], color, ...
            'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end


% =========================================================================
function addShadingLegend(ax, invalidColor, edgeColor)
% Add legend entries for shading colors on the given axes.
    hold(ax, 'on');
    patch(ax, NaN, NaN, invalidColor, 'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
        'DisplayName', 'Blanked (invalid)');
    patch(ax, NaN, NaN, edgeColor,    'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
        'DisplayName', 'Edge buffer');
    legend(ax, 'show', 'Location', 'northeast');
end


% =========================================================================
function saveFigure(fig, folderpath, baseName)
% Save figure as both .fig and .png.
    figPath = fullfile(folderpath, [baseName '.fig']);
    pngPath = fullfile(folderpath, [baseName '.png']);
    savefig(fig, figPath);
    exportgraphics(fig, pngPath, 'Resolution', 150);
    close(fig);
end