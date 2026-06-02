function D = step3_detect(D, P, plotMode)
% STEP3_DETECT  Threshold-crossing spike detection on the time-varying sigma.
%
%   D = step3_detect(D, P, plotMode)
%
% For each neural channel, detects events where the filtered signal crosses
% the LOCAL threshold P.threshSigma * sigma(t), enforcing a P.refractoryMs
% dead time and restricting detection to valid samples (D.validMask). The
% time-varying threshold is handled by normalising the signal by sigma(t)
% and finding peaks of the normalised trace above 1.
%
% Polarity (P.detectPolarity): 'neg' detects negative peaks (default for
% extracellular spikes), 'pos' positive, 'both' either.
%
% This is candidate detection only -- waveform extraction, alignment, and
% amplitude/width screening come in step 4. Peaks above P.maxThreshSigma are
% FLAGGED as likely artifacts (not removed here).
%
% Requires D.filtered (step 1), D.sigma and D.validMask (step 2).
%
% Adds to D:
%   D.spikes        1 x nNeural struct array, per channel:
%       .channel .label
%       .centers          sample indices of detected peaks
%       .times            spike times (s)
%       .peakAmp_uv       single-sided peak amplitude at detection (uV)
%       .threshAtSpike_uv local threshold at each spike (uV)
%       .artifactMask     true where peakAmp exceeds maxThreshSigma
%       .nSpikes .rate_hz .validSec
%   D.detectInfo    parameters used
%
% plotMode (single-file): per channel, a representative trace window with
% detected events marked, and the spike-amplitude histogram with the
% threshold marked (to reveal a hard cut-off = systematically missed
% low-amplitude events).

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    for f = {'filtered','sigma','validMask'}
        if ~isfield(D, f{1})
            error('step3_detect:missing', 'Run steps 1-2 first (missing D.%s).', f{1});
        end
    end

    fs   = D.fs;
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;
    pol  = lower(string(P.detectPolarity));
    refr = max(1, round(P.refractoryMs * 1e-3 * fs));

    spikes = repmat(emptySpikeStruct(), 1, nCh);

    fprintf('[step3] Detection: %.1f sigma, polarity=%s, refractory=%.2f ms (%d samp).\n', ...
        P.threshSigma, pol, P.refractoryMs, refr);

    for k = 1:nCh
        xf    = D.filtered(:, k);
        sig   = D.sigma(:, k);
        valid = D.validMask(:, k);
        thr   = P.threshSigma * sig;

        switch pol
            case "neg";  s = -xf;
            case "pos";  s =  xf;
            case "both"; s = abs(xf);
            otherwise; error('step3_detect:pol', 'Unknown DetectionPolarity.');
        end

        z = s ./ thr;                       % >1 means it crossed local threshold
        z(~valid | ~isfinite(z)) = -Inf;

        [~, locs] = findpeaks(z, 'MinPeakHeight', 1, 'MinPeakDistance', refr);
        locs = locs(:);

        peakAmp_uv       = abs(xf(locs)) * 1e6;
        threshAtSpike_uv = thr(locs) * 1e6;
        artifactMask     = abs(xf(locs)) > (P.maxThreshSigma * sig(locs));

        validSec = nnz(valid) / fs;
        nsp      = numel(locs);
        rate     = nsp / max(validSec, eps);

        spikes(k).channel          = ch(k);
        spikes(k).label            = lab{k};
        spikes(k).centers          = locs;
        spikes(k).times            = (locs - 1) / fs;
        spikes(k).peakAmp_uv       = peakAmp_uv;
        spikes(k).threshAtSpike_uv = threshAtSpike_uv;
        spikes(k).artifactMask     = artifactMask;
        spikes(k).nSpikes          = nsp;
        spikes(k).rate_hz          = rate;
        spikes(k).validSec         = validSec;

        fprintf(['        ch %d (%s): %d events | %.2f spk/s over %.1f s valid | ' ...
                 'median amp %.1f uV | %d flagged artifact(s)\n'], ...
            ch(k), lab{k}, nsp, rate, validSec, ...
            median(peakAmp_uv, 'omitnan'), nnz(artifactMask));
    end

    D.spikes = spikes;
    D.detectInfo = struct('threshSigma', P.threshSigma, 'polarity', char(pol), ...
                          'refractoryMs', P.refractoryMs, 'maxThreshSigma', P.maxThreshSigma);

    if plotMode
        for k = 1:nCh
            plot_detection(D, k, pol, P.threshSigma, P.maxThreshSigma, P.plotMaxPoints);
        end
    end
end

% ========================================================================
function r = emptySpikeStruct()
    r = struct('channel', NaN, 'label', '', 'centers', [], 'times', [], ...
        'peakAmp_uv', [], 'threshAtSpike_uv', [], 'artifactMask', [], ...
        'nSpikes', 0, 'rate_hz', NaN, 'validSec', NaN);
end

% ========================================================================
function plot_detection(D, k, pol, threshSigma, maxThreshSigma, maxPts) %#ok<INUSD>
    fs  = D.fs;
    t   = D.t;
    xf  = D.filtered(:, k);
    sig = D.sigma(:, k);
    thr = threshSigma * sig;
    sp  = D.spikes(k);
    lab = sp.label;

    fig = figure('Color', 'w', 'Name', sprintf('Step 3 — %s', lab), ...
        'Position', [160 120 1150 760]);
    tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 3 detection  |  %s  |  %d events, %.2f spk/s', ...
        lab, sp.nSpikes, sp.rate_hz), 'Interpreter', 'none');

    % ---- representative window: bin nearest the median spike count ----
    [winLo, winHi] = pick_window(sp.times, t, D.validMask(:, k), fs, 5);
    i0 = max(1, floor(winLo * fs) + 1);
    i1 = min(numel(t), floor(winHi * fs) + 1);

    nexttile;
    plot(t(i0:i1), xf(i0:i1) * 1e6, 'Color', [0.45 0.45 0.45], 'LineWidth', 0.5, ...
        'DisplayName', 'filtered'); hold on;
    if pol == "pos" || pol == "both"
        plot(t(i0:i1),  thr(i0:i1) * 1e6, 'r', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    if pol == "neg" || pol == "both"
        plot(t(i0:i1), -thr(i0:i1) * 1e6, 'r', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    inWin = sp.centers >= i0 & sp.centers <= i1;
    spc = sp.centers(inWin);
    if ~isempty(spc)
        plot(t(spc), xf(spc) * 1e6, 'r.', 'MarkerSize', 10, 'DisplayName', 'spikes');
    end
    % overlay R-peaks to check whether large regular events are cardiac-locked
    if isfield(D, 'rpeakTimes') && ~isempty(D.rpeakTimes)
        rp = D.rpeakTimes(D.rpeakTimes >= t(i0) & D.rpeakTimes <= t(i1));
        yl = ylim;
        for rr = 1:numel(rp)
            plot([rp(rr) rp(rr)], yl, ':', 'Color', [0 0.45 1], ...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        if ~isempty(rp)
            plot(nan, nan, ':', 'Color', [0 0.45 1], 'LineWidth', 0.8, ...
                'DisplayName', 'R-peaks');
        end
    end
    grid on; xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    title(sprintf('Representative %.0f s window with detected events (%d shown)', ...
        winHi - winLo, numel(spc)));
    legend('show', 'Location', 'northeast');
    xlim([t(i0) t(i1)]);

    % ---- amplitude histogram with threshold line ----
    nexttile;
    amp = sp.peakAmp_uv;
    if isempty(amp)
        text(0.4, 0.5, 'No events detected', 'Units', 'normalized'); axis off;
    else
        hiEdge = pct(amp, 0.99);
        if ~isfinite(hiEdge) || hiEdge <= 0; hiEdge = max(amp); end
        edges = linspace(0, hiEdge, 60);
        histogram(amp, edges); hold on;
        medThr = median(sp.threshAtSpike_uv, 'omitnan');
        xline(medThr, 'r--', 'LineWidth', 1.5, ...
            'Label', sprintf('median %.1f\\sigma thr (%.1f \\muV)', threshSigma, medThr));
        grid on; xlabel('Peak amplitude (\muV)'); ylabel('Count');
        title('Spike-amplitude histogram (hard wall at threshold = low-amplitude events missed)');
        xlim([0 hiEdge]);
    end
end

% ========================================================================
function [lo, hi] = pick_window(times, t, valid, fs, winSec)
% Choose a winSec window that is (a) essentially fully valid -- not sitting on
% blanked/dead/edge samples -- and (b) representative in spike count (closest
% to the median of non-empty bins, not the busiest stretch).
    T = t(end);
    N = numel(valid);
    if isempty(times)
        % fall back to the middle of the valid region
        vi = find(valid, 1, 'first'); ve = find(valid, 1, 'last');
        if isempty(vi); lo = 0; hi = min(winSec, T); return; end
        c = ((vi + ve) / 2 - 1) / fs;
        lo = max(0, c - winSec/2); hi = min(T, lo + winSec); return;
    end
    edges = 0:winSec:T;
    if edges(end) < T; edges(end+1) = T; end
    counts = histcounts(times, edges);
    nb = numel(counts);

    % valid fraction inside each candidate window (avoid blanked/dead regions)
    vfrac = zeros(1, nb);
    for bi = 1:nb
        i0 = max(1, floor(edges(bi)   * fs) + 1);
        i1 = min(N, floor(edges(bi+1) * fs) + 1);
        vfrac(bi) = mean(valid(i0:i1));
    end

    target = median(counts(counts > 0));
    cand = find(counts > 0 & vfrac >= 0.95);          % representative AND clean
    if isempty(cand); cand = find(counts > 0 & vfrac >= 0.80); end
    if isempty(cand); cand = find(counts > 0); end     % last resort
    [~, j] = min(abs(counts(cand) - target));
    bi = cand(j);
    lo = edges(bi); hi = min(T, edges(bi) + winSec);
end

% ========================================================================
function v = pct(x, q)
% q-th quantile via sort (no Statistics Toolbox dependency).
    x = x(isfinite(x));
    if isempty(x); v = NaN; return; end
    xs = sort(x);
    v = xs(max(1, min(numel(xs), round(q * numel(xs)))));
end
