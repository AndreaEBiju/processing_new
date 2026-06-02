function D = step2_noise_sigma(D, P, plotMode)
% STEP2_NOISE_SIGMA  Validity mask + time-varying robust noise estimate.
%
%   D = step2_noise_sigma(D, P, plotMode)
%
% Builds, per neural channel:
%   * a validity mask that excludes
%       - NaN-blanked samples,
%       - flat / zero runs (recording-error dropouts) lasting >= P.zeroRunMinSec,
%       - samples inside D.removedSegmentIdx,
%       - a P.edgeBufferMs pad around all of the above (filter ringing);
%   * a per-sample robust noise sigma track, sigma = median(|x|)/0.6745,
%     computed in P.sigmaWindowSec sliding windows over VALID samples only
%     and interpolated to every sample. This makes the detector adapt to a
%     non-stationary noise floor instead of using one global threshold.
%
% Requires D.filtered from step1_bandpass.
%
% Adds to D:
%   D.validMask   samples x nNeural   true where data is usable
%   D.sigma       samples x nNeural   per-sample noise sigma (volts)
%   D.sigmaWin    struct              per-window centres and sigma per channel
%   D.noiseInfo   struct              parameters + per-channel median sigma
%
% plotMode (single-file): per channel, shows the filtered signal with the
% +/- detection-threshold envelope, and the sigma track vs time, with all
% invalid regions shaded so dropouts and drift are obvious.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'filtered')
        error('step2_noise_sigma:noFiltered', 'Run step1_bandpass first.');
    end

    fs   = D.fs;
    N    = size(D.filtered, 1);
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;

    edgePad = round(P.edgeBufferMs * 1e-3 * fs);
    win     = max(2, round(P.sigmaWindowSec * fs));
    step    = max(1, round(win * P.sigmaStepFrac));

    % global mask from removedSegmentIdx (shared by all channels)
    removedMask = false(N, 1);
    if isfield(D, 'removedSegmentIdx') && ~isempty(D.removedSegmentIdx)
        segs = round(D.removedSegmentIdx);
        for i = 1:size(segs, 1)
            s1 = max(1, segs(i, 1));
            s2 = min(N, segs(i, 2));
            if s2 >= s1; removedMask(s1:s2) = true; end
        end
    end

    D.validMask = false(N, nCh);
    D.sigma     = nan(N, nCh);
    sigWinAll   = cell(1, nCh);
    cWinAll     = cell(1, nCh);
    medSigma    = nan(1, nCh);

    fprintf('[step2] Noise sigma: %.1f s windows (%.0f%% overlap), MAD estimator.\n', ...
        P.sigmaWindowSec, 100 * (1 - P.sigmaStepFrac));

    for k = 1:nCh
        xraw = D.y(:, ch(k));
        xf   = D.filtered(:, k);

        % --- invalid core: NaN, flat/dead runs, removed segments ---
        flatMask    = flat_run_mask(xraw, fs, P.zeroRunMinSec, P.flatRangeFrac);
        invalidCore = isnan(xf) | flatMask | removedMask;

        % --- dilate by edge pad to drop filter ringing at boundaries ---
        if edgePad > 0
            invalid = movmax(double(invalidCore), [edgePad edgePad]) > 0;
        else
            invalid = invalidCore;
        end
        valid = ~invalid;
        D.validMask(:, k) = valid;

        % --- sliding-window MAD sigma over valid samples ---
        nWin   = max(1, floor((N - win) / step) + 1);
        sigWin = nan(nWin, 1);
        cWin   = nan(nWin, 1);
        minValid = max(100, round(P.sigmaMinValidFrac * win));
        for w = 1:nWin
            i0 = (w - 1) * step + 1;
            i1 = min(N, i0 + win - 1);
            seg = xf(i0:i1);
            vv  = valid(i0:i1);
            s   = seg(vv & isfinite(seg));
            if numel(s) >= minValid
                sigWin(w) = median(abs(s)) / 0.6745;
            end
            cWin(w) = (i0 + i1) / 2;
        end

        % --- build per-sample sigma track ---
        good = isfinite(sigWin) & sigWin > 0;
        if nnz(good) >= 2
            sig = interp1(cWin(good), sigWin(good), (1:N)', 'linear');
            sig = fillmissing(sig, 'nearest');          % extend to ends
        elseif nnz(good) == 1
            sig = repmat(sigWin(good), N, 1);
        else
            % fallback: single global MAD over all valid samples
            s = xf(valid & isfinite(xf));
            gsig = median(abs(s)) / 0.6745;
            if ~isfinite(gsig) || gsig <= 0; gsig = std(s, 'omitnan'); end
            sig = repmat(gsig, N, 1);
            warning('step2:globalSigma', ...
                'ch %d (%s): too few valid windows; used a global sigma.', ch(k), lab{k});
        end
        sig(~isfinite(sig) | sig <= 0) = NaN;
        sig = fillmissing(sig, 'nearest');

        D.sigma(:, k) = sig;
        sigWinAll{k}  = sigWin;
        cWinAll{k}    = cWin;
        medSigma(k)   = median(sig(valid), 'omitnan');

        validSec = nnz(valid) / fs;
        fprintf(['        ch %d (%s): median sigma = %.2f uV | valid %.1f s ' ...
                 '(%.0f%% of record) | flagged flat: %.1f s\n'], ...
            ch(k), lab{k}, medSigma(k) * 1e6, validSec, ...
            100 * nnz(valid) / N, nnz(flatMask) / fs);
    end

    D.sigmaWin  = struct('centers', {cWinAll}, 'sigma', {sigWinAll}, ...
                         'windowSec', P.sigmaWindowSec, 'stepFrac', P.sigmaStepFrac);
    D.noiseInfo = struct('windowSec', P.sigmaWindowSec, 'edgeBufferMs', P.edgeBufferMs, ...
                         'zeroRunMinSec', P.zeroRunMinSec, 'threshSigma', P.threshSigma, ...
                         'medianSigma_uV', medSigma * 1e6);

    if plotMode
        for k = 1:nCh
            plot_sigma(D.t, D.filtered(:, k), D.sigma(:, k), D.validMask(:, k), ...
                fs, lab{k}, P.threshSigma, P.plotMaxPoints);
        end
    end
end

% ========================================================================
function mask = flat_run_mask(x, fs, minRunSec, rangeFrac)
% Flag dead / flat runs lasting >= minRunSec. Operates on the RAW signal.
% A region is "dead" when its local peak-to-peak range collapses far below
% the recording's typical local range -- this catches near-zero dropouts
% from recording errors, not just literal zeros or stuck-constant samples.
    n = numel(x);
    w = max(3, round(0.2 * fs));                 % 200 ms local window
    rng = movmax(x, w, 'omitnan') - movmin(x, w, 'omitnan');
    % Typical *live* local range. Use a high percentile, NOT the median:
    % the dead region can be the majority of the record, so the median would
    % fall inside it and drive the threshold to ~0. (Computed via sort to
    % avoid a Statistics-Toolbox dependency.)
    r = rng(isfinite(rng));
    if isempty(r)
        g = 1;
    else
        rs = sort(r);
        g = rs(max(1, round(0.90 * numel(rs))));   % 90th percentile
    end
    if ~isfinite(g) || g <= 0; g = max(r); end
    deadThresh = rangeFrac * g;

    flatpt = (rng <= deadThresh) | (x == 0);     % low-range OR exact zero
    flatpt(isnan(x)) = false;                    % NaNs handled separately

    minRun = max(1, round(minRunSec * fs));
    mask = false(n, 1);
    d = diff([0; flatpt(:); 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    for i = 1:numel(starts)
        if ends(i) - starts(i) + 1 >= minRun
            mask(starts(i):ends(i)) = true;
        end
    end
end

% ========================================================================
function plot_sigma(t, xf, sigma, valid, fs, label, threshSigma, maxPts)
    n = numel(t);
    stride = max(1, floor(n / maxPts));
    idx = 1:stride:n;

    fig = figure('Color', 'w', 'Name', sprintf('Step 2 — %s', label), ...
        'Position', [140 120 1100 720]);
    tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 2 noise sigma  |  %s', label), 'Interpreter', 'none');

    thr = threshSigma * sigma;

    % ---- signal with +/- threshold envelope ----
    ax1 = nexttile;
    plot(t(idx), xf(idx) * 1e6, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.4); hold on;
    plot(t(idx),  thr(idx) * 1e6, 'r', 'LineWidth', 0.8);
    plot(t(idx), -thr(idx) * 1e6, 'r', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    shade_runs(ax1, t, ~valid);
    grid on; xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    legend({'filtered', sprintf('\\pm%.1f\\sigma', threshSigma)}, ...
        'Location', 'northeast');
    title('Filtered signal with detection-threshold envelope');
    xlim([t(1) t(end)]);

    % ---- sigma track over time ----
    ax2 = nexttile;
    plot(t(idx), sigma(idx) * 1e6, 'b', 'LineWidth', 1.0); hold on;
    shade_runs(ax2, t, ~valid);
    grid on; xlabel('Time (s)'); ylabel('\sigma (\muV)');
    title('Robust noise \sigma vs time (flat = stationary; drift = non-stationary)');
    xlim([t(1) t(end)]);
end

% ========================================================================
function shade_runs(ax, t, mask)
% Shade contiguous true-runs of mask as translucent patches, without
% disturbing the y-limits set by the real data.
    if ~any(mask); return; end
    yl = ylim(ax);
    d = diff([0; mask(:); 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    maxRuns = 300;
    for i = 1:min(numel(starts), maxRuns)
        x0 = t(starts(i)); x1 = t(min(ends(i), numel(t)));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.4 0.4], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end
