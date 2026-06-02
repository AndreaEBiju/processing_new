function R = step5d_threshold_sweep(D, P, plotMode)
% STEP5D_THRESHOLD_SWEEP  Quick "are there sub-/supra-threshold populations?" check.
%
%   R = step5d_threshold_sweep(D, P, plotMode)
%
% Run AFTER step0 -> step1a -> step1 -> step2 on ONE file (needs D.filtered
% and D.validMask). This is an exploratory single-file diagnostic, not part of
% the production flow.
%
% Idea (the design we agreed on):
%   * Detection threshold is swept over k = 2..8 sigma on the blanked+filtered
%     signal. Amplitude is ONLY the detection axis.
%   * A phase-randomized SURROGATE (same power spectrum, no spikes) is detected
%     identically. excess = real_rate - surrogate_rate tells you how low you can
%     meaningfully go: where excess -> 0 is the empirical noise floor.
%   * Populations are tested on the TIMESCALE axis, not amplitude: per spike we
%     compute trough-to-peak duration (ms) and spectral centroid (Hz). The
%     Silverman unimodality p is computed on those features at each threshold,
%     for BOTH real and surrogate detections (surrogate = the noise smear floor).
%
% A real second population must be (a) multimodal in a timescale feature,
% (b) NOT multimodal in the surrogate, and (c) present where excess > 0.
%
% Returns R (struct) with the swept arrays; draws one figure per channel.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'filtered'); error('step5d:noFilt','Run step1_bandpass first.'); end

    fs   = D.fs;
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;
    refr = max(1, round(P.refractoryMs * 1e-3 * fs));
    pol  = lower(string(P.detectPolarity)); if pol == ""; pol = "neg"; end

    % ---- sweep + sizes (self-contained defaults) ----
    kSweep    = [2 2.5 3 3.5 4 4.5 5 6 8];
    kShow     = [3 4.5 6];                 % thresholds shown in the feature panels
    winS      = P.sigmaWindowSec;          % sigma(t) window (s)
    stepFrac  = P.sigmaStepFrac;
    minVFrac  = P.sigmaMinValidFrac;
    prePost   = [P.wfPreMs P.wfPostMs];    % waveform window (ms)
    bandHz    = [P.bandpassLow min(P.bandpassHigh, fs/2-1)];
    maxFeat   = 3000;                      % cap spikes used for per-spike features
    modMaxN   = 1500;                      % subsample for Silverman
    nBoot     = 80;                        % bootstrap reps (quick check)
    grid      = 256;

    R = struct('kSweep', kSweep, 'channel', num2cell(ch));

    fprintf('[step5d] Threshold sweep + phase-randomized surrogate (quick check).\n');

    for c = 1:nCh
        xf = D.filtered(:, c);
        if isfield(D, 'validMask'); valid = D.validMask(:, c);
        else;                       valid = isfinite(xf); end
        valid = valid & isfinite(xf);
        validSec = nnz(valid) / fs;

        % sigma(t) on the real signal, and a phase-randomized surrogate + its sigma(t)
        sigtR = local_sigma_t(xf,  valid, fs, winS, stepFrac, minVFrac);
        surr  = phase_randomize(xf, valid);
        sigtS = local_sigma_t(surr, valid, fs, winS, stepFrac, minVFrac);

        nK = numel(kSweep);
        [rateReal, rateSurr, nReal, nSurr] = deal(nan(1, nK));
        [pTTP_R, pTTP_S, pCEN_R, pCEN_S]   = deal(nan(1, nK));
        featR = cell(1, nK); featS = cell(1, nK);   % per-k feature structs (subsampled)

        for ik = 1:nK
            k = kSweep(ik);
            locsR = detect_k(xf,  valid, sigtR, k, refr, pol);
            locsS = detect_k(surr, valid, sigtS, k, refr, pol);
            nReal(ik) = numel(locsR);  nSurr(ik) = numel(locsS);
            rateReal(ik) = nReal(ik) / max(validSec, eps);
            rateSurr(ik) = nSurr(ik) / max(validSec, eps);

            fR = wf_features(xf,  sub_locs(locsR, maxFeat), valid, fs, prePost, bandHz, pol);
            fS = wf_features(surr, sub_locs(locsS, maxFeat), valid, fs, prePost, bandHz, pol);
            featR{ik} = fR; featS{ik} = fS;

            pTTP_R(ik) = silverman_p(fR.ttp_ms,    modMaxN, nBoot, grid);
            pTTP_S(ik) = silverman_p(fS.ttp_ms,    modMaxN, nBoot, grid);
            pCEN_R(ik) = silverman_p(fR.cen_hz,    modMaxN, nBoot, grid);
            pCEN_S(ik) = silverman_p(fS.cen_hz,    modMaxN, nBoot, grid);

            fprintf(['   ch %d (%s) k=%.1f: real %.2f spk/s | surr %.2f | excess %.2f | ' ...
                     'p(ttp) real %.3f surr %.3f\n'], ch(c), lab{c}, k, ...
                     rateReal(ik), rateSurr(ik), rateReal(ik)-rateSurr(ik), pTTP_R(ik), pTTP_S(ik));
        end

        % empirical floor: lowest k whose excess is still clearly positive
        excess = rateReal - rateSurr;
        frac   = excess ./ max(rateReal, eps);          % fraction of detections that beat noise
        kFloor = NaN; good = find(frac >= 0.5);          % >=50% excess
        if ~isempty(good); kFloor = kSweep(good(1)); end

        R(c).label    = lab{c};
        R(c).rateReal = rateReal;  R(c).rateSurr = rateSurr;  R(c).excess = excess;
        R(c).nReal    = nReal;     R(c).nSurr   = nSurr;      R(c).fracExcess = frac;
        R(c).pTTP_real = pTTP_R;   R(c).pTTP_surr = pTTP_S;
        R(c).pCEN_real = pCEN_R;   R(c).pCEN_surr = pCEN_S;
        R(c).kFloor    = kFloor;
        R(c).featReal  = featR;    R(c).featSurr = featS;

        fprintf('   ch %d (%s): empirical noise floor ~ %.2g sigma (excess>=50%%).\n', ...
            ch(c), lab{c}, kFloor);

        if plotMode
            plot_sweep(R(c), kSweep, kShow, lab{c});
        end
    end
end

% ======================================================================
% detection at threshold k*sigma(t) (same logic as step3, local copy)
function locs = detect_k(x, valid, sigt, k, refr, pol)
    thr = k * sigt;
    switch pol
        case "pos";  s =  x;
        case "both"; s =  abs(x);
        otherwise;   s = -x;          % 'neg' default
    end
    z = s ./ thr;
    z(~valid | ~isfinite(z)) = -Inf;
    [~, locs] = findpeaks(z, 'MinPeakHeight', 1, 'MinPeakDistance', refr);
    locs = locs(:);
end

function l = sub_locs(locs, maxN)
    if numel(locs) <= maxN; l = locs; return; end
    idx = randperm(numel(locs), maxN);
    l = sort(locs(idx));
end

% ======================================================================
% per-spike TIMESCALE features: trough-to-peak (ms) + spectral centroid (Hz)
function F = wf_features(x, locs, valid, fs, prePost, bandHz, pol)
    pre  = round(prePost(1) * 1e-3 * fs);
    post = round(prePost(2) * 1e-3 * fs);
    N = numel(x); L = pre + post + 1;
    ttp = []; cen = []; vpp = [];
    f = ((0:L-1) * fs / L)';                    % FFT bin frequencies (column)
    fb = f >= bandHz(1) & f <= bandHz(2);
    for i = 1:numel(locs)
        p = locs(i); a = p - pre; b = p + post;
        if a < 1 || b > N; continue; end
        if ~all(valid(a:b)); continue; end
        w = x(a:b);
        if any(~isfinite(w)); continue; end
        % trough-to-peak: from the extremum to the following opposite peak
        if pol == "pos"; [~, ie] = max(w); else; [~, ie] = min(w); end
        seg = w(ie:end);
        if pol == "pos"; [~, io] = min(seg); else; [~, io] = max(seg); end
        ttp(end+1,1) = (io - 1) / fs * 1000;    %#ok<AGROW>  ms
        vpp(end+1,1) = (max(w) - min(w)) * 1e6; %#ok<AGROW>  uV
        % spectral centroid of the waveform within the spike band
        Wm = abs(fft(w - mean(w)));
        num = sum(f(fb) .* Wm(fb)); den = sum(Wm(fb));
        cen(end+1,1) = num / max(den, eps); %#ok<AGROW>  Hz
    end
    F = struct('ttp_ms', ttp, 'cen_hz', cen, 'vpp_uv', vpp, 'n', numel(ttp));
end

% ======================================================================
% phase-randomized surrogate: same power spectrum, random phase (kills spikes)
function surr = phase_randomize(x, valid)
    x0 = x; x0(~valid) = NaN;
    x0 = fillmissing(x0, 'linear', 'EndValues', 'nearest');
    x0(~isfinite(x0)) = 0;
    n = numel(x0);
    X = fft(x0);
    ph = angle(fft(randn(n, 1)));               % real-noise phase => conj-symmetric
    surr = real(ifft(abs(X) .* exp(1i * ph)));
end

% ======================================================================
% sliding-window MAD sigma(t), interpolated to every sample (valid only)
function sigt = local_sigma_t(sig, valid, fs, winS, stepFrac, minVFrac)
    n = numel(sig);
    w = max(1, round(winS * fs));
    stp = max(1, round(w * stepFrac));
    if w >= n
        x = sig(valid); s = madsig(x);
        sigt = repmat(max(s, eps), n, 1); return;
    end
    ctr = []; sg = [];
    for a = 1:stp:(n - w + 1)
        b = a + w - 1;
        vv = valid(a:b);
        if mean(vv) < minVFrac; continue; end
        x = sig(a:b); x = x(vv);
        s = madsig(x);
        if s > 0 && isfinite(s); ctr(end+1) = (a + b) / 2; sg(end+1) = s; end %#ok<AGROW>
    end
    if numel(sg) < 2                              % 0 or 1 window -> constant sigma
        if isempty(sg); x = sig(valid); s = madsig(x); else; s = sg(1); end
        sigt = repmat(max(s, eps), n, 1); return;
    end
    sigt = interp1(ctr(:), sg(:), (1:n)', 'linear', 'extrap');
    sigt = max(sigt, eps);
end

function s = madsig(x)
    x = x(isfinite(x));
    if isempty(x); s = eps; return; end
    s = median(abs(x - median(x))) / 0.6745;
end

% ======================================================================
% Silverman smoothed-bootstrap unimodality p (compact, no toolbox deps)
% p small  -> reject unimodality (multimodal);  p large -> unimodal
function p = silverman_p(x, maxN, nBoot, grid)
    x = x(isfinite(x)); n = numel(x);
    if n < 30; p = NaN; return; end
    if n > maxN; x = x(randperm(n, maxN)); n = maxN; end
    x = x(:);
    sd = std(x); if sd <= 0; p = NaN; return; end
    gx = linspace(min(x) - 3*sd, max(x) + 3*sd, grid)';

    % critical bandwidth: smallest h giving a unimodal KDE (bisection)
    hHi = sd;                                  % wide -> unimodal
    while count_modes(x, hHi, gx) > 1; hHi = hHi * 1.5; end
    hLo = hHi / 50;
    for it = 1:40
        hMid = 0.5 * (hLo + hHi);
        if count_modes(x, hMid, gx) <= 1; hHi = hMid; else; hLo = hMid; end
    end
    hCrit = hHi;

    % smoothed bootstrap at hCrit (variance-corrected), count multimodal reps
    v = var(x); cnt = 0;
    for b = 1:nBoot
        ys = x(randi(n, n, 1));
        yb = ys + hCrit * randn(n, 1);
        yb = mean(x) + (yb - mean(x)) / sqrt(1 + hCrit^2 / v);
        if count_modes(yb, hCrit, gx) > 1; cnt = cnt + 1; end
    end
    p = cnt / nBoot;
end

function m = count_modes(x, h, gx)
    % Gaussian KDE on grid gx, count interior local maxima (vectorized)
    d = sum(exp(-0.5 * ((gx - x(:)') / h).^2), 2);   % grid x n -> grid x 1
    dd = diff(d);
    m = sum(dd(1:end-1) > 0 & dd(2:end) <= 0);
    if m < 1; m = 1; end
end

% ======================================================================
function plot_sweep(Rc, kSweep, kShow, label)
    fig = figure('Color','w','Name',sprintf('Step 5d sweep — %s', label), ...
        'Position',[120 90 1250 880]);
    tl = tiledlayout(fig, 3, 3, 'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('Threshold sweep + noise surrogate  |  %s  (floor ~ %.2g\\sigma)', ...
        label, Rc.kFloor), 'Interpreter','tex');

    % --- row 1: detection rate vs k (real / surrogate / excess) ---
    ax = nexttile([1 3]);
    semilogy(kSweep, max(Rc.rateReal, 1e-3), '-o', 'LineWidth', 1.6, ...
        'Color', [0.1 0.3 0.8], 'DisplayName', 'real'); hold on;
    semilogy(kSweep, max(Rc.rateSurr, 1e-3), '--s', 'LineWidth', 1.4, ...
        'Color', [0.8 0.2 0.2], 'DisplayName', 'surrogate (noise)');
    semilogy(kSweep, max(Rc.excess, 1e-3), '-^', 'LineWidth', 1.6, ...
        'Color', [0 0.5 0], 'DisplayName', 'excess (real-surr)');
    xline(4.5, ':', '4.5\sigma', 'Color', [0.4 0.4 0.4], 'HandleVisibility','off');
    if isfinite(Rc.kFloor); xline(Rc.kFloor, '-', sprintf('floor %.2g\\sigma', Rc.kFloor), ...
            'Color', [0 0.5 0], 'HandleVisibility','off'); end
    grid on; xlabel('threshold k (\sigma)'); ylabel('rate (spk/s, log)');
    title('Detection rate vs threshold: where excess \rightarrow 0 is the noise floor');
    legend('Location','northeast');

    % --- rows 2-3: timescale-feature distributions at kShow ---
    ik = arrayfun(@(kk) find(kSweep == kk, 1), kShow);
    for j = 1:numel(kShow)
        % trough-to-peak (row 2)
        nexttile;
        fR = Rc.featReal{ik(j)}; fS = Rc.featSurr{ik(j)};
        plot_feat(fR.ttp_ms, fS.ttp_ms, 'trough-to-peak (ms)', ...
            Rc.pTTP_real(ik(j)), Rc.pTTP_surr(ik(j)), kShow(j));
    end
    for j = 1:numel(kShow)
        % spectral centroid (row 3)
        nexttile;
        fR = Rc.featReal{ik(j)}; fS = Rc.featSurr{ik(j)};
        plot_feat(fR.cen_hz, fS.cen_hz, 'spectral centroid (Hz)', ...
            Rc.pCEN_real(ik(j)), Rc.pCEN_surr(ik(j)), kShow(j));
    end
end

function plot_feat(xr, xs, xlab, pR, pS, k)
    xr = xr(isfinite(xr)); xs = xs(isfinite(xs));
    if isempty(xr)
        text(0.5,0.5,'no spikes','Units','normalized','HorizontalAlignment','center');
        axis off; title(sprintf('k=%.1f\\sigma', k)); return;
    end
    lo = min([xr; xs]); hi = max([xr; xs]);
    if hi <= lo; hi = lo + 1; end
    edges = linspace(lo, hi, 30);
    histogram(xr, edges, 'Normalization','pdf', 'FaceColor',[0.2 0.4 0.8], ...
        'FaceAlpha',0.55, 'EdgeColor','none'); hold on;
    histogram(xs, edges, 'Normalization','pdf', 'DisplayStyle','stairs', ...
        'EdgeColor',[0.8 0.2 0.2], 'LineWidth',1.4);
    grid on; xlabel(xlab); ylabel('pdf');
    mod_r = ternary(pR < 0.05, 'MULTI', 'uni');
    title(sprintf('k=%.1f\\sigma:  real p=%.3f (%s)   surr p=%.3f', k, pR, mod_r, pS));
    if k == 3
        legend({'real','surrogate'}, 'Location','northeast', 'FontSize', 7);
    end
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end
