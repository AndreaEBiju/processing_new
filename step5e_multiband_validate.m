function R = step5e_multiband_validate(D, P, plotMode)
% STEP5E_MULTIBAND_VALIDATE  Colleague's 3-band detection + cardiac validation.
%
%   R = step5e_multiband_validate(D, P, plotMode)
%
% Replicates the multi-band ENG analysis (bands A/B/C) from the colleague's
% notebook ON YOUR DATA, then runs the test that notebook is missing: a
% peri-R-peak cross-correlogram for each band, computed on BOTH the unblanked
% signal (with_ecg, positive control) and the ±blank signal (no_ecg, the test).
%
% Pipeline per channel (matches the colleague's logic):
%   * valid mask from removedSegmentIdx + NaNs
%   * per-channel QRS realignment (local |extremum| within +/-30 ms of each R)
%   * blank +/- P.cardiacRemoveWinMs around the realigned QRS  (the no_ecg mask)
%   * per band: bandpass (NaN-safe, no blank-edge ringing) -> sigma(t) MAD ->
%     detect at P.threshSigma * sigma(t), neg polarity, P.refractoryMs refractory
%   * per-spike wide-band (100-5000) Vpp + FWHM, rate, ISI
%
% VALIDATION (the point of this script):
%   peri-R cross-correlogram of each band's spikes vs the R-peaks, +/-150 ms,
%   2 ms bins. Read it as:
%     - bump near 0 (within a few ms)            -> QRS / ringing (cardiac)
%     - bump at +20..100 ms (T-wave window)      -> T-wave (cardiac, survives a
%                                                   QRS-only blank!)
%     - flat across the whole cardiac cycle      -> NOT cardiac (real population)
%   A Poisson z of the largest cardiac-window bin vs the far-lag baseline gives
%   a scalar verdict per band/variant.
%
% Needs D.y (raw), D.fs, D.neuralChannels, D.channelLabels, D.rpeakSamples,
% and (optional) D.removedSegmentIdx. Respiration peri-event is added if
% D.breathSamples is present; otherwise skipped.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D,'y') || ~isfield(D,'rpeakSamples')
        error('step5e:missing','Needs D.y (raw) and D.rpeakSamples.');
    end

    fs   = D.fs;
    chs  = D.neuralChannels;
    nCh  = numel(chs);
    lab  = D.channelLabels;
    N    = size(D.y,1);

    % ---- params (colleague's values; clamp HF to Nyquist) ----
    nyq   = fs/2 - 1;
    bands = {'A',[1500 min(5000,nyq)]; 'B',[500 1500]; 'C',[100 500]};
    nB    = size(bands,1);
    featBand = [100 min(5000,nyq)];          % wide-band for Vpp/FWHM
    threshK  = P.threshSigma;                % 4.5
    refrMs   = P.refractoryMs;               % 1.0
    blankMs  = P.cardiacRemoveWinMs;         % 15
    realignMs= 30;
    wfHalfMs = 3.0;                          % WAVEFORM_WIDE_MS
    ord      = 4;
    % sigma(t)
    winS = P.sigmaWindowSec; stepFrac = P.sigmaStepFrac; minVF = P.sigmaMinValidFrac;
    % peri-event
    periHalfMs = 150; periBinMs = 2;
    twaveMs    = [20 100];                    % rat repolarization / T-wave window

    rpAll = round(D.rpeakSamples(:)); rpAll = rpAll(rpAll>=1 & rpAll<=N);
    refr  = max(1, round(refrMs*1e-3*fs));

    fprintf('[step5e] Multi-band (A/B/C) + peri-R validation | %d R-peaks | blank +/-%g ms\n', ...
        numel(rpAll), blankMs);

    R = struct('channel', num2cell(chs));

    for c = 1:nCh
        raw = double(D.y(:, chs(c)));
        valid0 = isfinite(raw);
        if isfield(D,'removedSegmentIdx') && ~isempty(D.removedSegmentIdx)
            seg = round(D.removedSegmentIdx);
            for i = 1:size(seg,1)
                a = max(1,seg(i,1)); b = min(N,seg(i,2)); if b>=a; valid0(a:b)=false; end
            end
        end
        % per-channel QRS realignment + blank mask
        rp = realign_qrs(raw, rpAll, fs, realignMs);
        validBlank = blank_mask(valid0, rp, fs, blankMs);

        % wide-band feature signals (one per variant)
        featW.with = bandpass_nan(raw, valid0,     fs, featBand(1), featBand(2), ord);
        featW.no   = bandpass_nan(raw, validBlank, fs, featBand(1), featBand(2), ord);

        variants = {'with_ecg', valid0; 'no_ecg', validBlank};
        Rc = struct(); Rc.label = lab{c};

        for vi = 1:2
            vname = variants{vi,1}; vmask = variants{vi,2};
            if vi==1; fwide = featW.with; else; fwide = featW.no; end
            validSec = nnz(vmask)/fs;

            for bi = 1:nB
                bn = bands{bi,1}; brg = bands{bi,2};
                bsig = bandpass_nan(raw, vmask, fs, brg(1), brg(2), ord);
                sigt = sigma_t(bsig, vmask, fs, winS, stepFrac, minVF);
                locs = detect_k(bsig, vmask, sigt, threshK, refr);

                % per-spike features on the wide-band signal
                [vpp, fwhm, wfm, twf] = wf_features(fwide, locs, vmask, fs, wfHalfMs);
                isi = diff(sort(locs))/fs*1000;          % ms
                M = struct();
                M.band = bn; M.range = brg; M.n = numel(locs);
                M.rate = numel(locs)/max(validSec,eps);
                M.sigma_uv = median(sigt(vmask),'omitnan')*1e6;
                M.vpp_med = median(vpp,'omitnan'); M.fwhm_med = median(fwhm,'omitnan');
                M.isi_med = median(isi,'omitnan');
                M.isi_cv  = std(isi,'omitnan')/max(mean(isi,'omitnan'),eps);
                M.meanWave = wfm; M.t_wf = twf; M.locs = locs;

                % ---- peri-R cross-correlogram (the validation) ----
                [lags, cc, base, pz, pkLagMs] = peri_event(locs, rp, fs, periHalfMs, periBinMs, twaveMs);
                M.periLags = lags; M.periCC = cc; M.periBase = base;
                M.cardiacZ = pz; M.cardiacPeakLagMs = pkLagMs;
                M.verdict = ternary(pz >= 4, 'CARDIAC-LOCKED', 'flat (not cardiac)');

                Rc.(vname).(bn) = M;

                fprintf(['   ch%d %-8s %s: N=%5d rate=%5.2f Hz | Vpp=%4.1f FWHM=%4.2fms ' ...
                         'ISImed=%5.1fms | peri-R z=%.1f @%+.0fms -> %s\n'], ...
                    chs(c), vname, bn, M.n, M.rate, M.vpp_med, M.fwhm_med, M.isi_med, ...
                    pz, pkLagMs, M.verdict);
            end
        end

        % ---- A vs B: one population or two? modality on cardiac-cleaned A union B ----
        % Union the no_ecg A and B detections (dedupe co-detections within the
        % refractory), measure each event's wide-band timescale (trough-to-peak +
        % spectral centroid), remove the cardiac-locked fraction (events at
        % over-represented R-peak phases), then Silverman-test the survivors.
        % Unimodal => one continuum (A/B are filter slices); bimodal => candidate
        % two types.
        Aloc = Rc.no_ecg.A.locs; Bloc = Rc.no_ecg.B.locs;
        locsAB = dedupe(sort([Aloc(:); Bloc(:)]), refr);
        capAB = 6000;
        if numel(locsAB) > capAB
            locsAB = locsAB(round(linspace(1, numel(locsAB), capAB)));
        end
        F = wf_timescale(featW.no, locsAB, validBlank, fs, wfHalfMs, featBand);
        used   = locsAB(F.idx);
        locked = cardiac_locked_events(used, rp, fs, periBinMs);
        keep   = ~locked;
        AB = struct();
        AB.nUnion = numel(locsAB); AB.nFeat = numel(F.idx); AB.nClean = nnz(keep);
        AB.ttp_all = F.ttp; AB.cen_all = F.cen;
        AB.ttp = F.ttp(keep); AB.cen = F.cen(keep);
        AB.pTTP = silverman_p(AB.ttp, 1500, 120, 256);
        AB.pCEN = silverman_p(AB.cen, 1500, 120, 256);
        AB.verdict = ternary(AB.pCEN < 0.05 || AB.pTTP < 0.05, ...
            'MULTIMODAL (candidate 2 types)', 'unimodal (one continuum)');
        Rc.AB = AB;
        fprintf(['   ch%d A&B union: %d events, %d after cardiac removal | ' ...
                 'modality p(ttp)=%.3f p(centroid)=%.3f -> %s\n'], ...
            chs(c), AB.nUnion, AB.nClean, AB.pTTP, AB.pCEN, AB.verdict);

        % optional respiration peri-event (band C, no_ecg) if breaths provided
        if isfield(D,'breathSamples') && ~isempty(D.breathSamples)
            bsig = bandpass_nan(raw, validBlank, fs, bands{3,2}(1), bands{3,2}(2), ord);
            sigt = sigma_t(bsig, validBlank, fs, winS, stepFrac, minVF);
            locsC = detect_k(bsig, validBlank, sigt, threshK, refr);
            br = round(D.breathSamples(:)); br = br(br>=1 & br<=N);
            [rl, rcc, rb, rz, rpk] = peri_event(locsC, br, fs, 2000, 25, [0 0]); %#ok<ASGLU>
            Rc.respC = struct('lags',rl,'cc',rcc,'base',rb,'z',rz);
            fprintf('   ch%d band C peri-respiration z=%.1f\n', chs(c), rz);
        end

        R(c).res = Rc;
        if plotMode
            plot_waveforms(Rc, lab{c});
            plot_periR(Rc, lab{c}, twaveMs);
            plot_AB(Rc, lab{c});
        end
    end
end

% ======================================================================
function rs2 = realign_qrs(raw, rs, fs, winMs)
    w = round(winMs/1000*fs); N = numel(raw); rs2 = rs;
    for i = 1:numel(rs)
        a = max(1, rs(i)-w); b = min(N, rs(i)+w);
        [~, o] = max(abs(raw(a:b))); rs2(i) = a + o - 1;
    end
end

function v = blank_mask(valid, rs, fs, blankMs)
    w = round(blankMs/1000*fs); v = valid; N = numel(valid);
    for i = 1:numel(rs)
        a = max(1, rs(i)-w); b = min(N, rs(i)+w); v(a:b) = false;
    end
end

% NaN-safe bandpass: interpolate across invalid (no blank-edge ringing), filter,
% then re-zero invalid. zp2sos for stability at low corners / high fs.
function y = bandpass_nan(x, valid, fs, lo, hi, ord)
    hi = min(hi, fs/2 - 1); lo = max(lo, 0.1);
    [z,p,k] = butter(ord, [lo hi]/(fs/2), 'bandpass');
    [sos,g] = zp2sos(z,p,k);
    xi = x; xi(~valid) = NaN;
    xi = fillmissing(xi, 'linear', 'EndValues','nearest');
    xi(~isfinite(xi)) = 0;
    y = filtfilt(sos, g, xi);
    y(~valid) = 0;
end

function locs = detect_k(x, valid, sigt, k, refr)
    thr = k * sigt;
    z = (-x) ./ thr;                     % negative-going
    z(~valid | ~isfinite(z)) = -Inf;
    [~, locs] = findpeaks(z, 'MinPeakHeight', 1, 'MinPeakDistance', refr);
    locs = locs(:);
end

function sigt = sigma_t(sig, valid, fs, winS, stepFrac, minVF)
    n = numel(sig); w = max(1, round(winS*fs)); stp = max(1, round(w*stepFrac));
    if w >= n
        s = madsig(sig(valid)); sigt = repmat(max(s,eps), n, 1); return;
    end
    ctr = []; sg = [];
    for a = 1:stp:(n-w+1)
        b = a+w-1; vv = valid(a:b);
        if mean(vv) < minVF; continue; end
        x = sig(a:b); s = madsig(x(vv));
        if s>0 && isfinite(s); ctr(end+1)=(a+b)/2; sg(end+1)=s; end %#ok<AGROW>
    end
    if numel(sg) < 2
        if isempty(sg); s = madsig(sig(valid)); else; s = sg(1); end
        sigt = repmat(max(s,eps), n, 1); return;
    end
    sigt = interp1(ctr(:), sg(:), (1:n)', 'linear', 'extrap');
    sigt = max(sigt, eps);
end

function s = madsig(x)
    x = x(isfinite(x)); if isempty(x); s = eps; return; end
    s = median(abs(x - median(x))) / 0.6745;
end

function [vpp, fwhm, meanWave, twf] = wf_features(x, locs, valid, fs, halfMs)
    h = round(halfMs*1e-3*fs); L = 2*h+1; Nn = numel(x);
    twf = ((-h:h)/fs)*1000;
    vpp = []; fwhm = []; W = [];
    cap = 4000; if numel(locs) > cap; locs = locs(round(linspace(1,numel(locs),cap))); end
    for i = 1:numel(locs)
        p = locs(i); a = p-h; b = p+h;
        if a<1 || b>Nn; continue; end
        if ~all(valid(a:b)); continue; end
        w = x(a:b); if any(~isfinite(w)); continue; end
        vpp(end+1,1) = (max(w)-min(w))*1e6;       %#ok<AGROW>
        fwhm(end+1,1) = fwhm_ms(w, fs);           %#ok<AGROW>
        W(end+1,:) = w(:)';                        %#ok<AGROW>
    end
    if isempty(W); meanWave = zeros(1,L); else; meanWave = mean(W,1)*1e6; end
end

function ms = fwhm_ms(w, fs)
    [tr, it] = min(w); half = tr/2;             % half-min level (negative spike)
    li = it; while li>1 && w(li)<=half; li=li-1; end
    ri = it; while ri<numel(w) && w(ri)<=half; ri=ri+1; end
    ms = (ri-li)/fs*1000;
end

% cross-correlogram of spikes vs events; returns counts, far-lag baseline,
% Poisson z of the largest bin in [-10, twave_hi] ms, and that peak's lag.
function [lags, cc, base, pz, pkLagMs] = peri_event(spk, ev, fs, halfMs, binMs, twaveMs)
    half = halfMs/1000; bin = binMs/1000;
    edges = -half:bin:half; lags = (edges(1:end-1)+bin/2)*1000;   % ms
    st = sort(spk(:))/fs; ev = sort(ev(:))/fs;
    cc = zeros(1, numel(lags));
    for i = 1:numel(ev)
        d = st - ev(i); d = d(d>=-half & d<half);
        if ~isempty(d); cc = cc + histcounts(d, edges); end
    end
    far = abs(lags) > 100;                       % baseline from far lags
    base = mean(cc(far)); if base<=0; base = mean(cc)+eps; end
    win = lags >= -10 & lags <= twaveMs(2);      % cardiac-relevant window
    [pk, j] = max(cc(win)); wl = lags(win);
    pz = (pk - base)/sqrt(max(base,1));
    pkLagMs = wl(j);
end

% ======================================================================
function plot_waveforms(Rc, label)
    figure('Color','w','Name',sprintf('Step 5e waveforms — %s', label), ...
        'Position',[120 360 1100 320]);
    tl = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('%s — band waveforms (no\\_ecg, blanked)', label));
    cols = {'A',[0.12 0.47 0.71]; 'B',[0.17 0.63 0.17]; 'C',[1 0.5 0.05]};
    for bi = 1:3
        bn = cols{bi,1}; M = Rc.no_ecg.(bn);
        nexttile;
        if M.n>0; plot(M.t_wf, M.meanWave, 'Color', cols{bi,2}, 'LineWidth', 2); end
        hold on; xline(0,'k:'); grid on; xlabel('ms'); ylabel('\muV');
        title(sprintf('%s (%d-%d Hz)  N=%d, %.1f Hz', bn, M.range(1), M.range(2), M.n, M.rate));
    end
end

function plot_periR(Rc, label, twaveMs)
    figure('Color','w','Name',sprintf('Step 5e peri-R — %s', label), ...
        'Position',[120 80 1200 640]);
    tl = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf(['%s — peri-R cross-correlogram (bump@0=QRS, bump@%g-%gms=T-wave, ' ...
        'flat=not cardiac)'], label, twaveMs(1), twaveMs(2)));
    variants = {'with_ecg','no_ecg'}; bandsN = {'A','B','C'};
    for vi = 1:2
        for bi = 1:3
            M = Rc.(variants{vi}).(bandsN{bi});
            nexttile;
            bar(M.periLags, M.periCC, 1, 'FaceColor',[0.4 0.55 0.8], 'EdgeColor','none'); hold on;
            yl = ylim;
            patch([twaveMs(1) twaveMs(2) twaveMs(2) twaveMs(1)], [0 0 yl(2) yl(2)], ...
                [1 0.8 0.4], 'FaceAlpha',0.25, 'EdgeColor','none');
            yline(M.periBase, 'r--', 'LineWidth',1);
            xline(0,'k:');
            grid on; xlabel('lag from R (ms)'); ylabel('count');
            title(sprintf('%s | %s:  z=%.1f  %s', variants{vi}, bandsN{bi}, M.cardiacZ, M.verdict), ...
                'Interpreter','none');
            xlim([-150 150]);
        end
    end
end

function locs = dedupe(locs, refr)
% merge events closer than refr samples (A/B co-detections of one spike)
    if isempty(locs); return; end
    locs = sort(locs(:));
    locs = locs([true; diff(locs) >= refr]);
end

function F = wf_timescale(x, locs, valid, fs, halfMs, bandHz)
% per-event wide-band trough-to-peak (ms) + spectral centroid (Hz); returns the
% indices of locs actually used (others skipped for bounds/validity).
    h = round(halfMs*1e-3*fs); Nn = numel(x); L = 2*h+1;
    f = ((0:L-1)*fs/L)'; fb = f>=bandHz(1) & f<=bandHz(2);
    ttp = []; cen = []; idx = [];
    for i = 1:numel(locs)
        p = locs(i); a = p-h; b = p+h;
        if a<1 || b>Nn; continue; end
        if ~all(valid(a:b)); continue; end
        w = x(a:b); if any(~isfinite(w)); continue; end
        [~, ie] = min(w); seg = w(ie:end); [~, io] = max(seg);
        ttp(end+1,1) = (io-1)/fs*1000;                 %#ok<AGROW>
        Wm = abs(fft(w - mean(w)));
        num = sum(f(fb).*Wm(fb)); den = sum(Wm(fb));
        cen(end+1,1) = num/max(den,eps);               %#ok<AGROW>
        idx(end+1,1) = i;                              %#ok<AGROW>
    end
    F = struct('ttp', ttp, 'cen', cen, 'idx', idx);
end

function lag = nearest_lag(st, ev)
% signed lag (s) from each spike time st to its nearest event time ev
    lag = nan(size(st));
    if isempty(ev); return; end
    ev = sort(ev(:));
    k = interp1(ev, 1:numel(ev), st, 'nearest', 'extrap');
    k = min(max(round(k),1), numel(ev));
    lag = st - ev(k);
    for d = [-1 1]
        j = min(max(k+d,1), numel(ev));
        alt = st - ev(j);
        better = abs(alt) < abs(lag);
        lag(better) = alt(better);
    end
end

function locked = cardiac_locked_events(spk, ev, fs, binMs)
% classify each spike as cardiac-locked if its nearest-R lag falls in a bin
% that is significantly over-represented (vs the median phase level).
    st = sort(spk(:))/fs; evs = sort(ev(:))/fs;
    lag = nearest_lag(st, evs);
    half = 0.150; bin = binMs/1000;
    edges = -half:bin:half;
    inr = abs(lag) < half;
    counts = histcounts(lag(inr), edges);
    nz = counts(counts>0);
    baseLvl = median(nz); if isempty(baseLvl) || baseLvl<=0; baseLvl = mean(counts)+eps; end
    sigbin = counts > baseLvl + 4*sqrt(max(baseLvl,1));
    bi = discretize(lag, edges);
    locked = false(size(st));
    ok = inr & ~isnan(bi);
    locked(ok) = sigbin(bi(ok));
end

% --- Silverman smoothed-bootstrap unimodality p (compact, no toolbox deps) ---
function p = silverman_p(x, maxN, nBoot, grid)
    x = x(isfinite(x)); n = numel(x);
    if n < 30; p = NaN; return; end
    if n > maxN; x = x(randperm(n, maxN)); n = maxN; end
    x = x(:); sd = std(x); if sd<=0; p = NaN; return; end
    gx = linspace(min(x)-3*sd, max(x)+3*sd, grid)';
    hHi = sd; guard = 0;
    while count_modes(x, hHi, gx) > 1 && guard < 60; hHi = hHi*1.5; guard = guard+1; end
    hLo = hHi/50;
    for it = 1:40
        hMid = 0.5*(hLo+hHi);
        if count_modes(x, hMid, gx) <= 1; hHi = hMid; else; hLo = hMid; end
    end
    hCrit = hHi; v = var(x); cnt = 0;
    for b = 1:nBoot
        ys = x(randi(n, n, 1)); yb = ys + hCrit*randn(n, 1);
        yb = mean(x) + (yb - mean(x))/sqrt(1 + hCrit^2/v);
        if count_modes(yb, hCrit, gx) > 1; cnt = cnt + 1; end
    end
    p = cnt / nBoot;
end

function m = count_modes(x, h, gx)
    d = sum(exp(-0.5 * ((gx - x(:)')/h).^2), 2);
    dd = diff(d);
    m = sum(dd(1:end-1) > 0 & dd(2:end) <= 0);
    if m < 1; m = 1; end
end

function plot_AB(Rc, label)
    if ~isfield(Rc, 'AB'); return; end
    AB = Rc.AB;
    figure('Color','w','Name',sprintf('Step 5e A&B modality — %s', label), ...
        'Position',[140 140 1020 420]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('%s — A\\cupB cardiac-cleaned (%d \\rightarrow %d events): %s', ...
        label, AB.nFeat, AB.nClean, AB.verdict), 'Interpreter','tex');
    nexttile; feat_hist(AB.ttp_all, AB.ttp, 'trough-to-peak (ms)', AB.pTTP);
    nexttile; feat_hist(AB.cen_all, AB.cen, 'spectral centroid (Hz)', AB.pCEN);
end

function feat_hist(xall, xkeep, xlab, p)
    xall = xall(isfinite(xall)); xkeep = xkeep(isfinite(xkeep));
    if isempty(xkeep)
        text(0.5,0.5,'no events','Units','normalized','HorizontalAlignment','center');
        axis off; return;
    end
    lo = min(xall); hi = max(xall); if hi<=lo; hi = lo+1; end
    e = linspace(lo, hi, 30);
    histogram(xall,  e, 'Normalization','pdf', 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none'); hold on;
    histogram(xkeep, e, 'Normalization','pdf', 'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.6, 'EdgeColor','none');
    grid on; xlabel(xlab); ylabel('pdf');
    md = ternary(p<0.05, 'MULTI', 'uni');
    title(sprintf('p=%.3f (%s)', p, md));
    legend({'all A\cupB','cardiac-cleaned'}, 'Location','best', 'FontSize',7);
end

function out = ternary(c,a,b); if c; out=a; else; out=b; end; end
