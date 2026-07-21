function results = HR_BR_HRVAnalysis_new(signal, fs, cutoff, order, folderpath, condition, chanidx, blankIdx, edgeBufferSec, winSec, stepSec, figtrue, hrBrWinSec)
% HR_BR_HRVAnalysis
%
% Analyze heart rate, breath rate, and HRV on a single signal matrix.
% (Previously supported a 'stim_rec' struct mode — that branch was removed;
% callers are expected to slice the signal themselves before calling.)
%
% Inputs:
%   signal        - [N x C] signal matrix
%   fs            - sampling rate in Hz (e.g. 25000)
%   cutoff        - lowpass filter cutoff in Hz (e.g. 8 Hz)
%   order         - lowpass filter order (e.g. 4)
%                   Filter settling time = order/cutoff seconds.
%                   edgeBufferSec must exceed this value.
%   folderpath    - output folder for saved .mat files and figures
%   condition     - string label for this condition
%   chanidx       - channel index to analyze
%   blankIdx      - Nx2 matrix of [start end] sample indices to blank
%   edgeBufferSec - seconds to mask around blank boundaries and signal edges
%                   default = 0.75s (covers order=4, cutoff=8Hz settling of 0.5s)
%   winSec        - window size (seconds) for HRV metrics and per-window beat
%                   count series. User-facing analysis window. Default 60.
%   stepSec       - step size in seconds for all moving metrics (default: 1)
%   figtrue       - logical, generate and save diagnostic figures (default: false)
%   hrBrWinSec    - window size (seconds) for heart rate and breath rate
%                   moving series. Default 60. Independent of winSec so HRV
%                   can be examined at a shorter user-chosen window while HR/BR
%                   stay on a stable 60s window.
%
% Diagnostic figures (saved as .fig and .png):
%   Figure 1 — Detection overview (6 panels):
%     1. Raw signal with invalid (gray) and edge-buffer-only (light yellow) regions shaded
%     2. Filtered+detrended signal with heart peaks marked
%     3. Signal at heartbeat times with breath troughs marked
%     4. Heart rate time series (bpm)
%     5. Heart beats per user-provided winSec window (raw count)
%     6. Breath rate time series
%   Figure 2 — HRV metrics (7 panels):
%     1. Filtered+detrended signal with heart peaks (reference)
%     2. Signal at heartbeat times with breath troughs (reference)
%     3. HRV (std RR) time series             [centered winSec window]
%     4. RMSSD time series                    [centered winSec window]
%     5. pNN5 time series                     [centered winSec window]
%     6. SD1 and SD2 overlaid                 [centered winSec window]
%     7. Sample entropy time series           [centered fixed 60s window]
%
% Moving metrics (all CENTERED windows of size W, same stepSec):
%   Heart rate   : peaks per valid minute in centered hrBrWinSec window (bpm)
%   Breath rate  : continuous-stretch based (>=5s stretch, >=3 peaks required);
%                  centered hrBrWinSec window. Per-window ISI plausibility
%                  gate removed (see movingCardiacMetrics for rationale);
%                  global br_implausibleFraction still flags QC issues.
%   Heart count  : raw heart beats within centered user-provided winSec window
%                  (beats / winSec). Saved as heartCountSeries.
%   HRV (std RR) : centered winSec window; requires >= minRR valid RR intervals
%   RMSSD, pNN5, SD1, SD2 : centered winSec window; same guard as HRV
%   Sample entropy : centered FIXED 60s window (independent of winSec).
%                    Requires >= 40 RR intervals in the window for stability.
%
% At each time t, a metric is NaN if:
%   (a) The window centered at t cannot fit fully inside the signal, i.e.
%       t < W/2 or t > T - W/2 where W is that metric's window size.
%   (b) The sample at t is inside the invalidMask (don't report a metric at
%       a time the underlying signal is blanked).
%
% minRR scales with winSec at 10% of expected RR count at ~400 bpm rat HR.

    if nargin < 8  || isempty(blankIdx),       blankIdx      = [];   end
    if nargin < 9  || isempty(edgeBufferSec),  edgeBufferSec = 0.75; end
    if nargin < 10 || isempty(winSec),         winSec        = 60;   end
    if nargin < 11 || isempty(stepSec),        stepSec       = 1;    end
    if nargin < 12 || isempty(figtrue),        figtrue       = false;end
    if nargin < 13 || isempty(hrBrWinSec),     hrBrWinSec    = 60;   end

    results = analyzeOneSignal(signal, fs, cutoff, order, folderpath, ...
        condition, chanidx, blankIdx, edgeBufferSec, winSec, stepSec, figtrue, hrBrWinSec);
end


% =========================================================================
function out = analyzeOneSignal(signal, fs, cutoff, order, folderpath, condition, ...
    chanidx, blankIdx, edgeBufferSec, winSec, stepSec, figtrue, hrBrWinSec)

    % --- Input validation ---
    assert(fs >= 100 && fs <= 30000, ...
        'fs = %.1f Hz seems implausible. Check units.', fs);

    settlingTimeSec = order / cutoff;
    if edgeBufferSec < settlingTimeSec
        warning('HR_BR_HRVAnalysis: edgeBufferSec (%.2fs) < filter settling time (%.2fs). Consider increasing it.', ...
            edgeBufferSec, settlingTimeSec);
    end

    condition      = string(condition);
    condition      = condition(1);
    folderpath     = char(string(folderpath));
    conditionLabel = char(condition);

    validateattributes(signal, {'numeric'}, {'2d','nonempty'});
    if chanidx > size(signal,2)
        error('chanidx exceeds number of channels.');
    end
    if isempty(blankIdx)
        blankIdx = zeros(0,2);
    end

    x = signal(:, chanidx);
    N = length(x);
    t = (0:N-1)' / fs;

    % --- minRR: 10% of expected RR intervals at rat HR ~400 bpm ---
    expectedRR = 400 * (winSec / 60);
    minRR      = max(3, round(0.1 * expectedRR));

    % --- Breath rate plausibility bounds ---
    minBreathRate_bpm = 40;
    maxBreathRate_bpm = 170;
    minBreathRateHz   = minBreathRate_bpm / 60;
    maxBreathRateHz   = maxBreathRate_bpm / 60;

    % ---------------- Valid / invalid masks ----------------
    invalidMask = isnan(x);

    if ~isempty(blankIdx)
        blankIdx      = round(blankIdx);
        blankIdx(:,1) = max(blankIdx(:,1), 1);
        blankIdx(:,2) = min(blankIdx(:,2), N);
        for k = 1:size(blankIdx,1)
            invalidMask(blankIdx(k,1):blankIdx(k,2)) = true;
        end
    end
    
    fprintf('N = %d samples (%.1f s at fs=%.1f)\n', N, N/fs, fs);
    fprintf('invalidMask covers %d samples (%.1f s)\n', sum(invalidMask), sum(invalidMask)/fs);
    if ~isempty(blankIdx)
        fprintf('blankIdx: %d segments, first=[%d %d], last=[%d %d]\n', ...
            size(blankIdx,1), blankIdx(1,1), blankIdx(1,2), ...
            blankIdx(end,1), blankIdx(end,2));
        fprintf('blankIdx max index = %d vs N = %d\n', max(blankIdx(:,2)), N);
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

    % ---------------- Signal preparation ----------------
    xFill = fillmissing(x, 'linear', 'EndValues', 'nearest');

    % Filter commented out — using raw signal directly
    % highFreqCutoff = cutoff / (fs/2);
    % [z,p,k_]       = butter(order, highFreqCutoff, 'low');
    % [sos,g]        = zp2sos(z, p, k_);
    % yFilt          = filtfilt(sos, g, xFill);
    yFilt = xFill;

    % Detrend on fully populated signal (no NaNs) for clean peak detection.
    % NaNs are reinserted afterward for display/bookkeeping only.
    heartBeatSeriesClean = detrend(yFilt);

    % ---------------- Heart peaks ----------------
    minPeakDist_samp = round(0.1 * fs);   % max ~600 bpm

    [~, heartlocsRaw] = findpeaks(heartBeatSeriesClean, ...
        'MinPeakDistance', minPeakDist_samp);

    % Reject peaks in invalid or edge regions
    validHeartPeakMask = ~invalidMask(heartlocsRaw) & ~edgeMask(heartlocsRaw);
    heartlocs          = heartlocsRaw(validHeartPeakMask);

    % Store with NaNs for display only
    heartBeatSeries = heartBeatSeriesClean;
    heartBeatSeries(invalidMask) = NaN;

    heartPeakTrain = zeros(N,1);
    heartPeakTrain(heartlocs)   = 1;
    heartPeakTrain(invalidMask) = 0;

    % ---------------- RR intervals ----------------
    [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs);

    if ~isempty(RR_intervals)
        minRR_sec              = 0.1;
        maxRR_sec              = 0.5;
        implausibleRR          = RR_intervals < minRR_sec | RR_intervals > maxRR_sec;
        RR_implausibleFraction = mean(implausibleRR);
        if RR_implausibleFraction > 0.05
            warning('%s: %.1f%% of RR intervals outside plausible range — check peak detection.', ...
                conditionLabel, RR_implausibleFraction * 100);
        end
    else
        RR_implausibleFraction = NaN;
    end

    % ---------------- Breath peaks ----------------
    % Breath detected as troughs of raw signal at heartbeat locations,
    % reflecting respiratory modulation of the vagal/gastric signal amplitude.
    br_implausibleFraction = NaN;
    br_locs_true           = [];
    envelopeAtHeartlocs    = [];

    if ~isempty(heartlocs)
        candidateVals       = xFill(heartlocs);
        envelopeAtHeartlocs = candidateVals;

        if ~all(candidateVals == 0)
            meanRR_samp       = mean(diff(heartlocs));
            minBreathSepSec   = 1 / maxBreathRateHz;
            minBreathSepBeats = max(2, round(minBreathSepSec * fs / meanRR_samp));

            % Detect troughs (minima) by inverting signal
            [~, br_locs] = findpeaks(-candidateVals, 'MinPeakDistance', minBreathSepBeats);
            br_locs_true = heartlocs(br_locs);
        end

        br_locs_true = br_locs_true(~invalidMask(br_locs_true) & ~edgeMask(br_locs_true));

        % Flag implausible breath intervals (do not discard)
        if numel(br_locs_true) >= 2
            br_intervals_sec       = diff(br_locs_true) / fs;
            br_rate_hz             = 1 ./ br_intervals_sec;
            implausibleBreath      = br_rate_hz < minBreathRateHz | br_rate_hz > maxBreathRateHz;
            br_implausibleFraction = mean(implausibleBreath);
            if br_implausibleFraction > 0.05
                warning('%s: %.1f%% of breath intervals outside plausible range — check breath detection.', ...
                    conditionLabel, br_implausibleFraction * 100);
            end
        end
    end

    % ---------------- Moving metrics ----------------
    brMinStretchSec     = 5;    % breath: 3 events at 40 bpm min -> 4.5s, rounded up
    brMinPeaksInStretch = 3;

    % HR stretch: at ~400 bpm, 3 beats occur in ~0.45s.
    % Use 1s minimum stretch — enough for 3+ beats and avoids
    % reporting rate from a single isolated valid sample cluster.
    hrMinStretchSec     = 1;
    hrMinPeaksInStretch = 3;

    [metrics_t, heartRateSeries, heartCountSeries, hrv_series, rmssd_series, pnn5_series, ...
        sd1_series, sd2_series, sampEn_series, nRR_used, breathRateSeries] = ...
        movingCardiacMetrics(heartPeakTrain, RR_intervals, RR_times, ...
            br_locs_true, invalidMask, edgeMask, fs, winSec, hrBrWinSec, stepSec, minRR, ...
            brMinStretchSec, brMinPeaksInStretch, minBreathRate_bpm, maxBreathRate_bpm, ...
            hrMinStretchSec, hrMinPeaksInStretch);

    avgHeartRate       = mean(heartRateSeries,  'omitnan');
    avgBreathRate      = mean(breathRateSeries, 'omitnan');
    avgHeartCount      = mean(heartCountSeries, 'omitnan');

    fprintf('%s Heart Rate (centered %ds window):  %.2f bpm\n', ...
        conditionLabel, hrBrWinSec, avgHeartRate);
    fprintf('%s Heart Beats per %ds window (centered): %.2f beats\n', ...
        conditionLabel, winSec, avgHeartCount);
    fprintf('%s Breath Rate (centered %ds window): %.2f bpm\n', ...
        conditionLabel, hrBrWinSec, avgBreathRate);

    % ---------------- Scalar HRV summary metrics ----------------
    if isempty(RR_intervals) || numel(RR_intervals) < 2
        hrv    = NaN;   rmssd  = NaN;   pnn5   = NaN;
        sd1    = NaN;   sd2    = NaN;
        sampEn = NaN;   appxEn = NaN;
    else
        hrv    = std(RR_intervals, 'omitnan');
        diffRR = diff(RR_intervals);
        rmssd  = sqrt(mean(diffRR.^2, 'omitnan'));

        absdiffRR    = abs(diff(RR_intervals .* 1000));
        count_over_5 = sum(absdiffRR > 5, 'omitnan');
        pnn5         = (count_over_5 / numel(absdiffRR)) * 100;

        sd1 = sqrt(0.5) * std(diffRR, 'omitnan');
        sd2 = sqrt(max(0, 2*hrv^2 - sd1^2));

        sampEn = sampleEntropy(RR_intervals);
        appxEn = approximateEntropy(RR_intervals);

        fprintf('%s HRV: %.4f s  RMSSD: %.4f s  pNN5: %.2f%%  SD1: %.4f  SD2: %.4f\n', ...
            conditionLabel, hrv, rmssd, pnn5, sd1, sd2);
    end

    % ---------------- Gap-aware DFA (fractal scaling) ----------------
    if isempty(RR_intervals) || numel(RR_intervals) < 2
        dfaOut = struct('alpha1', NaN, 'alpha2', NaN, 'alphaFull', NaN, ...
            'R2_1', NaN, 'R2_2', NaN, 'nCross', NaN, 'nWindows', [], ...
            'excludedScales', []);
    else
        dfaOut = dfaRR_gapAware(RR_intervals, RR_times, pipeline_params());
    end

    % ---------------- Diagnostic figures ----------------
    if figtrue
        invalidColor = [0.6 0.6 0.6];   % gray   — blanked/invalid
        edgeColor    = [1.0 0.95 0.6];  % yellow — edge buffer only

        % -- Figure 1: Detection overview --
        fig1 = figure('Name', sprintf('%s — Detection', conditionLabel), ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.92]);

        ax1 = gobjects(6,1);
        for p = 1:6
            ax1(p) = subplot(6,1,p, 'Parent', fig1);
        end

        % Panel 1: Raw signal
        plot(ax1(1), t, x, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
        ylabel(ax1(1), 'Raw signal');
        title(ax1(1), sprintf('%s — Raw signal (chan %d)', conditionLabel, chanidx));
        applyShading(ax1(1), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 2: Detrended signal + heart peaks
        plot(ax1(2), t, heartBeatSeries, 'Color', [0.15 0.45 0.75], 'LineWidth', 0.7);
        hold(ax1(2), 'on');
        if ~isempty(heartlocs)
            plot(ax1(2), t(heartlocs), heartBeatSeries(heartlocs), ...
                'v', 'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], ...
                'MarkerSize', 3, 'LineStyle', 'none');
        end
        ylabel(ax1(2), 'Detrended (a.u.)');
        title(ax1(2), sprintf('Heart peaks detected: %d', numel(heartlocs)));
        applyShading(ax1(2), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 3: Signal at heartbeat times + breath troughs
        if ~isempty(heartlocs) && ~isempty(envelopeAtHeartlocs)
            tHeartlocs = t(heartlocs);
            plot(ax1(3), tHeartlocs, envelopeAtHeartlocs, ...
                '.-', 'Color', [0.2 0.65 0.35], 'MarkerSize', 4, 'LineWidth', 0.6);
            hold(ax1(3), 'on');
            if ~isempty(br_locs_true)
                [~, brEnvIdx] = ismember(br_locs_true, heartlocs);
                brEnvIdx      = brEnvIdx(brEnvIdx > 0);
                plot(ax1(3), t(br_locs_true), envelopeAtHeartlocs(brEnvIdx), ...
                    '^', 'Color', [0.85 0.5 0.1], 'MarkerFaceColor', [0.85 0.5 0.1], ...
                    'MarkerSize', 5, 'LineStyle', 'none');
            end
        else
            text(0.5, 0.5, 'No heartlocs — envelope unavailable', ...
                'Parent', ax1(3), 'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
        ylabel(ax1(3), 'Signal at beats (a.u.)');
        title(ax1(3), sprintf('Breath troughs detected: %d', numel(br_locs_true)));
        applyShading(ax1(3), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 4: Heart rate series (centered hrBrWinSec window, bpm)
        plot(ax1(4), metrics_t, heartRateSeries, 'Color', [0.85 0.15 0.15], 'LineWidth', 1);
        ylabel(ax1(4), 'HR (bpm)');
        title(ax1(4), sprintf('Heart rate  (centered %ds window, mean: %.1f bpm)', ...
            hrBrWinSec, avgHeartRate));
        applyShading(ax1(4), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 5: Heart beats per user-provided winSec window (centered, raw count)
        plot(ax1(5), metrics_t, heartCountSeries, 'Color', [0.5 0.1 0.45], 'LineWidth', 1);
        ylabel(ax1(5), sprintf('Beats / %ds', winSec));
        title(ax1(5), sprintf('Heart beats per %ds window  (centered, mean: %.1f beats)', ...
            winSec, avgHeartCount));
        applyShading(ax1(5), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 6: Breath rate series (centered hrBrWinSec window)
        plot(ax1(6), metrics_t, breathRateSeries, 'Color', [0.85 0.5 0.1], 'LineWidth', 1);
        ylabel(ax1(6), 'BR (bpm)');
        xlabel(ax1(6), 'Time (s)');
        title(ax1(6), sprintf('Breath rate  (centered %ds window, mean: %.1f bpm)', ...
            hrBrWinSec, avgBreathRate));
        applyShading(ax1(6), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        linkaxes(ax1, 'x');
        xlim(ax1(1), [t(1) t(end)]);
        addShadingLegend(ax1(1), invalidColor, edgeColor);
        saveFigure(fig1, folderpath, sprintf('%s_detection', conditionLabel));

        % -- Figure 2: HRV metrics --
        fig2 = figure('Name', sprintf('%s — HRV metrics', conditionLabel), ...
            'NumberTitle', 'off', 'Color', 'w', ...
            'Units', 'normalized', 'Position', [0.05 0.02 0.9 0.95]);

        ax2 = gobjects(7,1);
        for p = 1:7
            ax2(p) = subplot(7,1,p, 'Parent', fig2);
        end

        % Panel 1: Detrended signal + heart peaks (reference)
        plot(ax2(1), t, heartBeatSeries, 'Color', [0.15 0.45 0.75], 'LineWidth', 0.7);
        hold(ax2(1), 'on');
        if ~isempty(heartlocs)
            plot(ax2(1), t(heartlocs), heartBeatSeries(heartlocs), ...
                'v', 'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], ...
                'MarkerSize', 3, 'LineStyle', 'none');
        end
        ylabel(ax2(1), 'Detrended (a.u.)');
        title(ax2(1), 'Heart signal + peaks (reference)');
        applyShading(ax2(1), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 2: Signal at heartbeat times + breath troughs (reference)
        if ~isempty(heartlocs) && ~isempty(envelopeAtHeartlocs)
            tHeartlocs = t(heartlocs);
            plot(ax2(2), tHeartlocs, envelopeAtHeartlocs, ...
                '.-', 'Color', [0.2 0.65 0.35], 'MarkerSize', 4, 'LineWidth', 0.6);
            hold(ax2(2), 'on');
            if ~isempty(br_locs_true)
                [~, brEnvIdx] = ismember(br_locs_true, heartlocs);
                brEnvIdx      = brEnvIdx(brEnvIdx > 0);
                plot(ax2(2), t(br_locs_true), envelopeAtHeartlocs(brEnvIdx), ...
                    '^', 'Color', [0.85 0.5 0.1], 'MarkerFaceColor', [0.85 0.5 0.1], ...
                    'MarkerSize', 5, 'LineStyle', 'none');
            end
        else
            text(0.5, 0.5, 'No heartlocs — envelope unavailable', ...
                'Parent', ax2(2), 'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
        ylabel(ax2(2), 'Signal at beats (a.u.)');
        title(ax2(2), 'Breath troughs (reference)');
        applyShading(ax2(2), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 3: HRV (std RR) — centered winSec window
        plot(ax2(3), metrics_t, hrv_series, 'Color', [0.4 0.2 0.7], 'LineWidth', 1);
        ylabel(ax2(3), 'HRV (s)');
        title(ax2(3), sprintf('HRV std(RR)  (centered %ds window, global: %.4f s)', ...
            winSec, hrv));
        applyShading(ax2(3), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 4: RMSSD — centered winSec window
        plot(ax2(4), metrics_t, rmssd_series, 'Color', [0.2 0.6 0.7], 'LineWidth', 1);
        ylabel(ax2(4), 'RMSSD (s)');
        title(ax2(4), sprintf('RMSSD  (centered %ds window, global: %.4f s)', winSec, rmssd));
        applyShading(ax2(4), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 5: pNN5 — centered winSec window
        plot(ax2(5), metrics_t, pnn5_series, 'Color', [0.7 0.4 0.1], 'LineWidth', 1);
        ylabel(ax2(5), 'pNN5 (%)');
        title(ax2(5), sprintf('pNN5  (centered %ds window, global: %.2f%%)', winSec, pnn5));
        applyShading(ax2(5), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 6: SD1 and SD2 — centered winSec window
        plot(ax2(6), metrics_t, sd1_series, 'Color', [0.2 0.5 0.8], 'LineWidth', 1);
        hold(ax2(6), 'on');
        plot(ax2(6), metrics_t, sd2_series, 'Color', [0.8 0.2 0.3], 'LineWidth', 1);
        legend(ax2(6), 'SD1', 'SD2', 'Location', 'northeast');
        ylabel(ax2(6), 'SD (s)');
        title(ax2(6), sprintf('Poincare  (centered %ds window) SD1: %.4f s   SD2: %.4f s', ...
            winSec, sd1, sd2));
        applyShading(ax2(6), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        % Panel 7: Sample entropy — centered fixed 60s window
        sampEn_global = sampEn;  % scalar global SampEn over all RR
        avgSampEn     = mean(sampEn_series, 'omitnan');
        plot(ax2(7), metrics_t, sampEn_series, 'Color', [0.45 0.25 0.55], 'LineWidth', 1);
        ylabel(ax2(7), 'SampEn');
        xlabel(ax2(7), 'Time (s)');
        title(ax2(7), sprintf( ...
            'Sample Entropy  (centered 60s window, mean: %.3f  global: %.3f)', ...
            avgSampEn, sampEn_global));
        applyShading(ax2(7), t, invalidMask, edgeOnlyMask, invalidColor, edgeColor);

        linkaxes(ax2, 'x');
        xlim(ax2(1), [t(1) t(end)]);
        addShadingLegend(ax2(1), invalidColor, edgeColor);
        saveFigure(fig2, folderpath, sprintf('%s_HRVmetrics', conditionLabel));
    end

    % ---------------- Save outputs ----------------
    hrbrFile = fullfile(folderpath, sprintf('%s_HRBR.mat', conditionLabel));
    hrvFile  = fullfile(folderpath, sprintf('%s_HRVMeasures.mat', conditionLabel));

    save(hrbrFile, ...
        'chanidx', 'heartBeatSeries', 'heartlocs', ...
        'heartRateSeries', 'avgHeartRate', ...
        'heartCountSeries', 'avgHeartCount', ...
        'breathRateSeries', 'avgBreathRate', ...
        'br_locs_true', 'br_implausibleFraction', ...
        'invalidMask', 'edgeMask', 'blankIdx', 'edgeBufferSec', ...
        'metrics_t', 't', 'winSec', 'hrBrWinSec', 'stepSec', ...
        'brMinStretchSec', 'brMinPeaksInStretch', ...
        'minBreathRate_bpm', 'maxBreathRate_bpm', ...
        'hrMinStretchSec', 'hrMinPeaksInStretch');

    sampEnWinSec = 60;  % must match value used inside movingCardiacMetrics
    dfa_alpha1         = dfaOut.alpha1;
    dfa_alpha2         = dfaOut.alpha2;
    dfa_alphaFull      = dfaOut.alphaFull;
    dfa_R2_1           = dfaOut.R2_1;
    dfa_R2_2           = dfaOut.R2_2;
    dfa_nCross         = dfaOut.nCross;
    dfa_nWindows       = dfaOut.nWindows;
    dfa_excludedScales = dfaOut.excludedScales;
    save(hrvFile, ...
        'hrv', 'rmssd', 'pnn5', 'sd1', 'sd2', ...
        'sampEn', 'appxEn', ...
        'RR_intervals', 'RR_times', 'RR_implausibleFraction', ...
        'heartlocs', 'invalidMask', 'edgeMask', 'blankIdx', 'edgeBufferSec', ...
        'metrics_t', 'hrv_series', 'rmssd_series', 'pnn5_series', ...
        'sd1_series', 'sd2_series', 'sampEn_series', 'sampEnWinSec', 'nRR_used', ...
        't', 'winSec', 'hrBrWinSec', 'stepSec', 'minRR', ...
        'dfa_alpha1', 'dfa_alpha2', 'dfa_alphaFull', 'dfa_R2_1', 'dfa_R2_2', ...
        'dfa_nCross', 'dfa_nWindows', 'dfa_excludedScales');

    % ---------------- Output struct ----------------
    out                        = struct();
    out.chanidx                = chanidx;
    out.heartBeatSeries        = heartBeatSeries;
    out.heartlocs              = heartlocs;
    out.metrics_t              = metrics_t;
    out.heartRateSeries        = heartRateSeries;
    out.avgHeartRate           = avgHeartRate;
    out.heartCountSeries       = heartCountSeries;
    out.avgHeartCount          = avgHeartCount;
    out.breathRateSeries       = breathRateSeries;
    out.avgBreathRate          = avgBreathRate;
    out.br_locs_true           = br_locs_true;
    out.br_implausibleFraction = br_implausibleFraction;
    out.RR_intervals           = RR_intervals;
    out.RR_times               = RR_times;
    out.RR_implausibleFraction = RR_implausibleFraction;
    out.hrv                    = hrv;
    out.rmssd                  = rmssd;
    out.pnn5                   = pnn5;
    out.sd1                    = sd1;
    out.sd2                    = sd2;
    out.dfa_alpha1             = dfaOut.alpha1;
    out.dfa_alpha2             = dfaOut.alpha2;
    out.dfa_alphaFull          = dfaOut.alphaFull;
    out.dfa_R2_1               = dfaOut.R2_1;
    out.dfa_R2_2               = dfaOut.R2_2;
    out.dfa_nCross             = dfaOut.nCross;
    out.dfa_nWindows           = dfaOut.nWindows;
    out.dfa_excludedScales     = dfaOut.excludedScales;
    out.hrv_series             = hrv_series;
    out.rmssd_series           = rmssd_series;
    out.pnn5_series            = pnn5_series;
    out.sd1_series             = sd1_series;
    out.sd2_series             = sd2_series;
    out.sampEn_series          = sampEn_series;
    out.sampEnWinSec           = sampEnWinSec;
    out.nRR_used               = nRR_used;
    out.sampEn                 = sampEn;
    out.appxEn                 = appxEn;
    out.invalidMask            = invalidMask;
    out.edgeMask               = edgeMask;
    out.edgeOnlyMask           = edgeOnlyMask;
    out.blankIdx               = blankIdx;
    out.t                      = t;
    out.winSec                 = winSec;
    out.stepSec                = stepSec;
    out.minRR                  = minRR;
    out.brMinStretchSec        = brMinStretchSec;
    out.brMinPeaksInStretch    = brMinPeaksInStretch;
    out.minBreathRate_bpm      = minBreathRate_bpm;
    out.maxBreathRate_bpm      = maxBreathRate_bpm;
    out.hrMinStretchSec        = hrMinStretchSec;
    out.hrMinPeaksInStretch    = hrMinPeaksInStretch;
end


% =========================================================================
function [metrics_t, heartRateSeries, heartCountSeries, hrv_series, rmssd_series, pnn5_series, ...
    sd1_series, sd2_series, sampEn_series, nRR_used, breathRateSeries] = ...
    movingCardiacMetrics(heartPeakTrain, RR_intervals, RR_times, ...
        br_locs_true, invalidMask, edgeMask, fs, winSec, hrBrWinSec, stepSec, minRR, ...
        brMinStretchSec, brMinPeaksInStretch, minBreathRate_bpm, maxBreathRate_bpm, ...
        hrMinStretchSec, hrMinPeaksInStretch)
% Compute moving cardiac metrics on CENTERED windows.
%
% Window sizes:
%   hrBrWinSec      : heart rate (bpm), breath rate (bpm). Default 60s.
%   winSec          : HRV (std RR), RMSSD, pNN5, SD1, SD2, per-window beat count.
%                     User-chosen analysis window.
%   sampEnWinSec=60 : sample entropy series, fixed 60s window
%                     (hardcoded; SampEn needs a stable, long enough RR
%                     embedding to be meaningful — keeping it independent of
%                     the user's HRV winSec).
%
% At each time t = metrics_t(i):
%   - Each metric is NaN if its centered window cannot fit fully inside the
%     signal: t < W/2 or t > T - W/2, where W is that metric's window.
%   - Each metric is NaN if the center sample t is inside invalidMask
%     (don't report metrics at a time the signal itself is blanked).
%
% Rate logic (HR, BR) still uses the longest contiguous clean stretch within
% the centered window so partial blanks inside a window don't fabricate beats
% from no signal. This matches the design in pipeline_summary §3.1.2 / §3.1.8;
% only the window centering changes.

    N         = numel(heartPeakTrain);
    sigDurSec = (N - 1) / fs;
    metrics_t = (0 : stepSec : sigDurSec)';
    nT        = numel(metrics_t);

    heartRateSeries  = NaN(nT, 1);
    heartCountSeries = NaN(nT, 1);
    hrv_series       = NaN(nT, 1);
    rmssd_series     = NaN(nT, 1);
    pnn5_series      = NaN(nT, 1);
    sd1_series       = NaN(nT, 1);
    sd2_series       = NaN(nT, 1);
    sampEn_series    = NaN(nT, 1);
    nRR_used         = zeros(nT, 1);
    breathRateSeries = NaN(nT, 1);

    brMinStretchSamp = round(brMinStretchSec * fs);
    hrMinStretchSamp = round(hrMinStretchSec * fs);
    cleanMask        = ~invalidMask & ~edgeMask;

    halfHrBr     = hrBrWinSec / 2;
    halfWin      = winSec / 2;
    sampEnWinSec = 60;           % fixed by design — see header docstring
    halfSampEn   = sampEnWinSec / 2;
    % Min RR count for a usable SampEn estimate. 10% of expected RR at rat HR
    % ~400 bpm over 60s = 40. Below this, SampEn estimates are unstable.
    minRR_sampEn = max(3, round(0.1 * 400 * (sampEnWinSec / 60)));

    for i = 1:nT
        tc       = metrics_t(i);
        idxCtr   = max(1, min(N, round(tc * fs) + 1));

        % Skip everything when the center itself is in a blanked region.
        if invalidMask(idxCtr)
            continue;
        end

        % ---- HR / BR : centered hrBrWinSec window ----
        if tc >= halfHrBr && tc <= sigDurSec - halfHrBr
            t0   = tc - halfHrBr;
            t1   = tc + halfHrBr;
            idx0 = max(1, round(t0 * fs) + 1);
            idx1 = min(N, round(t1 * fs) + 1);

            winClean = cleanMask(idx0:idx1);
            [stretchLen, stretchStart] = longestContiguousRun(winClean);

            % Heart rate (bpm) from longest clean stretch in the window.
            if stretchLen >= hrMinStretchSamp
                hrAbsStart       = idx0 + stretchStart - 1;
                hrAbsEnd         = min(N, hrAbsStart + stretchLen - 1);
                hrPeaksInStretch = sum(heartPeakTrain(hrAbsStart:hrAbsEnd));
                if hrPeaksInStretch >= hrMinPeaksInStretch
                    heartRateSeries(i) = hrPeaksInStretch / ((stretchLen / fs) / 60);
                end
            end

            % Breath rate (bpm) from longest clean stretch.
            % Per-window inter-peak plausibility gate intentionally removed
            % (symmetrical to slow wave fix). One stray missed breath (rate
            % slightly outside [40, 170] bpm) would otherwise nuke a whole
            % window even though >98% of intervals in the window are fine.
            % Global br_implausibleFraction (computed above) still reports
            % the fraction of intervals exceeding the plausibility band for QC.
            if ~isempty(br_locs_true) && stretchLen >= brMinStretchSamp
                absStart       = idx0 + stretchStart - 1;
                absEnd         = min(N, absStart + stretchLen - 1);
                peaksInStretch = br_locs_true(br_locs_true >= absStart & br_locs_true <= absEnd);

                if numel(peaksInStretch) >= brMinPeaksInStretch
                    breathRateSeries(i) = numel(peaksInStretch) / ((stretchLen / fs) / 60);
                end
            end
        end

        % ---- Per-window beat count + HRV : centered winSec window ----
        if tc >= halfWin && tc <= sigDurSec - halfWin
            t0w   = tc - halfWin;
            t1w   = tc + halfWin;
            idx0w = max(1, round(t0w * fs) + 1);
            idx1w = min(N, round(t1w * fs) + 1);

            % Heart beats per user-provided window: literal count over the
            % full centered window. heartPeakTrain zeros out invalid samples
            % so blanked portions contribute zero.
            heartCountSeries(i) = sum(heartPeakTrain(idx0w:idx1w));

            % HRV metrics: RR intervals whose start time is within the window.
            keep        = RR_times >= t0w & RR_times <= t1w;
            nRR_used(i) = sum(keep);
            if nRR_used(i) >= minRR
                RR_win  = RR_intervals(keep);
                diffRR  = diff(RR_win);
                hrv_val = std(RR_win, 'omitnan');

                hrv_series(i)   = hrv_val;
                rmssd_series(i) = sqrt(mean(diffRR.^2, 'omitnan'));

                absdiffRR      = abs(diffRR * 1000);
                pnn5_series(i) = (sum(absdiffRR > 5, 'omitnan') / numel(absdiffRR)) * 100;

                sd1_val        = sqrt(0.5) * std(diffRR, 'omitnan');
                sd1_series(i)  = sd1_val;
                sd2_series(i)  = sqrt(max(0, 2 * hrv_val^2 - sd1_val^2));
            end
        end

        % ---- Sample entropy : centered FIXED 60s window ----
        % Independent of winSec because SampEn needs a long, stable RR
        % embedding to give meaningful estimates.
        if tc >= halfSampEn && tc <= sigDurSec - halfSampEn
            t0se  = tc - halfSampEn;
            t1se  = tc + halfSampEn;
            keepSE = RR_times >= t0se & RR_times <= t1se;
            nSE    = sum(keepSE);
            if nSE >= minRR_sampEn
                sampEn_series(i) = sampleEntropyFast( ...
                    RR_intervals(keepSE), 2, []);
            end
        end
    end
end


% =========================================================================
function se = sampleEntropyFast(rr, m, r)
% Vectorized SampEn (Chebyshev distance) — equivalent to sampleEntropy.m but
% O(N^2) memory via pdist2 in one shot instead of a double for-loop.
%   rr : 1D series (e.g. RR intervals)
%   m  : embedding dimension (default 2)
%   r  : tolerance (default 0.2 * std(rr))
%
% Returns NaN when no matches are found (consistent with sampleEntropy.m).

    if nargin < 2 || isempty(m), m = 2; end
    rr = rr(:)';
    N  = numel(rr);
    if nargin < 3 || isempty(r), r = 0.2 * std(rr); end
    if N < m + 2 || ~isfinite(r) || r <= 0
        se = NaN;
        return;
    end

    Nm = N - m;
    % Embedding matrices: Nm rows, each a length-m or length-(m+1) window
    Xm  = zeros(Nm, m);
    Xm1 = zeros(Nm, m + 1);
    for k = 1:m
        Xm(:, k) = rr(k : k + Nm - 1).';
    end
    for k = 1:(m + 1)
        Xm1(:, k) = rr(k : k + Nm - 1).';
    end

    Dm  = pdist2(Xm,  Xm,  'chebychev');
    Dm1 = pdist2(Xm1, Xm1, 'chebychev');

    mask = triu(true(Nm), 1);   % i < j pairs only — matches sampleEntropy.m
    B = sum(Dm(mask)  < r);
    A = sum(Dm1(mask) < r);

    if B == 0 || A == 0
        se = NaN;
    else
        se = -log(A / B);
    end
end


% =========================================================================
function [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs)
    if numel(heartlocs) < 2
        RR_intervals = []; RR_times = []; return;
    end
    csumInvalid  = [0; cumsum(invalidMask(:))];
    RR_intervals = [];
    RR_times     = [];
    for i = 1:numel(heartlocs)-1
        s1 = heartlocs(i);
        s2 = heartlocs(i+1);
        if s2 <= s1 + 1
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1)     = s1 / fs;         %#ok<AGROW>
            continue;
        end
        if csumInvalid(s2) - csumInvalid(s1+1) == 0
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1)     = s1 / fs;         %#ok<AGROW>
        end
    end
end


% =========================================================================
function [runLen, runStart] = longestContiguousRun(binaryVec)
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
    shadeRegions(ax, t, edgeOnlyMask, edgeColor);
    shadeRegions(ax, t, invalidMask,  invalidColor);
end


% =========================================================================
function shadeRegions(ax, t, mask, color)
    if ~any(mask), return; end
    yl = ylim(ax);
    if yl(1) == yl(2), yl = [-1 1]; end
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
    hold(ax, 'on');
    patch(ax, NaN, NaN, invalidColor, 'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
        'DisplayName', 'Blanked (invalid)');
    patch(ax, NaN, NaN, edgeColor,    'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
        'DisplayName', 'Edge buffer');
    legend(ax, 'show', 'Location', 'northeast');
end


% =========================================================================
function saveFigure(fig, folderpath, baseName)
    figPath = fullfile(folderpath, [baseName '.fig']);
    pngPath = fullfile(folderpath, [baseName '.png']);
    savefig(fig, figPath);
    exportgraphics(fig, pngPath, 'Resolution', 150);
    close(fig);
end