function D = step6_spike_report(D, P, plotMode)
% STEP6_SPIKE_REPORT  Single-population spike-train report (per channel).
%
%   D = step6_spike_report(D, P, plotMode)
%
% Treats all accepted spikes (step 4) as ONE C-like population (no
% clustering) and computes the standard spike-train characterisation, per
% neural channel, storing every scalar for the bulk stage in D.metrics(k):
%
%   firing rate (spk/s, blanking-corrected) ; mean waveform + percentile band ;
%   peak-to-peak amplitude over time ; ISI: log-ISI histogram, time x ISI
%   heatmap, return map ; CV / CV2 / LV regularity ; Fano-vs-window curve ;
%   autocorrelogram ; firing-rate power spectrum ; amplitude-vs-ISI ;
%   data-driven burst detection (log-ISI void parameter).
%
% Requires D.spikes(k).alignedTimes / .alignedCenters / .waveforms /
% .Vpp_uv / .width_ms (step 4) and D.validMask (step 2).

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D,'spikes') || ~isfield(D.spikes,'alignedTimes')
        error('step6_spike_report:missing','Run step4_waveforms first.');
    end

    fs = D.fs; N = size(D.filtered,1); Tend = (N-1)/fs;
    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    cond = ''; if isfield(D,'condition'); cond = D.condition; end

    D.metrics = repmat(emptyMetrics(), 1, nCh);

    for k = 1:nCh
        valid = D.validMask(:,k);
        validDur = nnz(valid)/fs;
        % time-order spikes and carry features along
        [st, ord] = sort(D.spikes(k).alignedTimes(:));
        cen = D.spikes(k).alignedCenters(:); cen = cen(ord);
        vpp = D.spikes(k).Vpp_uv(:); vpp = vpp(ord);
        W = D.spikes(k).waveforms; if ~isempty(W); W = W(ord,:); end
        tms = D.spikes(k).wf_t_ms;
        isi = diff(st);                        % seconds
        isims = isi*1000;

        % gap-clean mask: an ISI is "clean" only if the interval between its
        % two spikes contains no invalid (blanked / cardiac-censored) sample.
        cvv = [0; cumsum(double(valid))];
        if numel(cen) >= 2
            span = double(cen(2:end) - cen(1:end-1));
            vbet = cvv(cen(2:end)+1) - cvv(cen(1:end-1)+1);
            isClean = vbet == span;
        else
            isClean = false(0,1);
        end

        M = emptyMetrics();
        M.label = lab{k}; M.condition = cond; M.nSpikes = numel(st);
        M.validDur_s = validDur; M.meanRate_hz = numel(st)/max(validDur,eps);
        M.medianVpp_uv = median(vpp,'omitnan'); M.medianFWHM_ms = median(D.spikes(k).width_ms,'omitnan');
        M.nISItotal = numel(isi); M.nISIclean = nnz(isClean);
        M.fracISIclean = M.nISIclean / max(numel(isi),1);

        % --- firing-rate trace (blanking-corrected) ---
        [M.fr_t, M.fr_hz, M.fr_validFrac] = firing_rate(cen, valid, N, fs, P.frBinSec);

        % --- regularity (gap-clean ISIs only) ---
        if nnz(isClean) >= 2
            d = isims;
            both = isClean(1:end-1) & isClean(2:end);     % consecutive clean ISI pairs
            if any(both)
                cv2v = 2*abs(d(2:end)-d(1:end-1))./(d(2:end)+d(1:end-1));
                lvv  = 3*(d(1:end-1)-d(2:end)).^2 ./ (d(1:end-1)+d(2:end)).^2;
                M.CV2 = mean(cv2v(both),'omitnan');
                M.LV  = mean(lvv(both),'omitnan');
            end
            M.CV  = std(d(isClean),'omitnan')/mean(d(isClean),'omitnan');
        end
        M.refracViolFrac = mean(isims(isClean) < P.refractoryMs);

        % --- Fano vs window ---
        [M.fano_T, M.fano_F] = fano_curve(cen, valid, N, fs, P);
        M.fanoCanon = interp1(M.fano_T, M.fano_F, P.fanoCanonSec, 'linear', NaN);
        gd = isfinite(M.fano_F) & M.fano_F>0;
        if nnz(gd) >= 3; pp = polyfit(log10(M.fano_T(gd)), log10(M.fano_F(gd)), 1); M.fanoSlope = pp(1); end

        % --- autocorrelogram (gap-clean pairs only) ---
        [M.acg_lag, M.acg] = autocorrelogram(st, cen, cvv, P.acgMaxLagSec, P.acgBinSec);

        % --- firing-rate power spectrum (largest valid segment) ---
        [M.psd_f, M.psd_p] = rate_psd(cen, valid, N, fs, P);

        % --- bursts (data-driven log-ISI, gap-clean) ---
        M.burst = detect_bursts(st, isims, isClean, validDur, P);

        % --- rolling CV2 (gap-clean) ---
        [M.cv2_t, M.cv2_roll] = rolling_cv2(st, isims, isClean, Tend, P.cv2WinSec);

        % --- waveform band ---
        if ~isempty(W)
            M.meanWaveform = mean(W,1,'omitnan');
            pb = prctile(W, P.wfBandPct, 1); M.bandLo = pb(1,:); M.bandHi = pb(2,:);
        end
        M.wf_t_ms = tms;

        D.metrics(k) = M;

        fprintf(['[step6] ch %d (%s): %d spk | rate %.2f /s | CV2 %.2f | LV %.2f | ' ...
                 'Fano@%.2gs %.2f (slope %.2f) | ISIclean %d/%d (%.0f%%) | bursts: %s\n'], ...
            ch(k), lab{k}, M.nSpikes, M.meanRate_hz, M.CV2, M.LV, P.fanoCanonSec, M.fanoCanon, ...
            M.fanoSlope, M.nISIclean, M.nISItotal, 100*M.fracISIclean, burstLine(M.burst));
        if M.fracISIclean < 0.2
            warning('step6:fewISI', ['ch %d (%s): only %.0f%% of ISIs are gap-clean -- ' ...
                'ISI-family metrics (CV2/LV/ISI-hist/bursts/ACG) are data-limited here.'], ...
                ch(k), lab{k}, 100*M.fracISIclean);
        end

        if plotMode
            plot_report(D, k, st, isims, isClean, vpp, valid);
        end
    end
end

% ======================================================================
function M = emptyMetrics()
    M = struct('label','','condition','','nSpikes',0,'validDur_s',NaN,'meanRate_hz',NaN, ...
        'medianVpp_uv',NaN,'medianFWHM_ms',NaN,'CV',NaN,'CV2',NaN,'LV',NaN,'refracViolFrac',NaN, ...
        'nISItotal',0,'nISIclean',0,'fracISIclean',NaN, ...
        'fanoCanon',NaN,'fanoSlope',NaN,'fano_T',[],'fano_F',[],'fr_t',[],'fr_hz',[],'fr_validFrac',[], ...
        'acg_lag',[],'acg',[],'psd_f',[],'psd_p',[],'cv2_t',[],'cv2_roll',[], ...
        'meanWaveform',[],'bandLo',[],'bandHi',[],'wf_t_ms',[],'burst',emptyBurst());
end
function B = emptyBurst()
    B = struct('hasBursts',false,'thrMs',NaN,'void',NaN,'nBursts',0,'rate_per_min',0, ...
        'meanDur_s',NaN,'meanSpikes',NaN,'intraRate_hz',NaN,'fracInBurst',0,'onsets',[],'offsets',[]);
end
function s = burstLine(B)
    if B.hasBursts; s = sprintf('%d (%.1f/min, %.0f%% spikes)', B.nBursts, B.rate_per_min, 100*B.fracInBurst);
    else; s = 'none'; end
end

% ======================================================================
function [t, fr, validFrac] = firing_rate(cen, valid, N, fs, binSec)
    binN = max(1, round(binSec*fs)); nB = ceil(N/binN);
    edges = (0:nB)*binN; t = nan(nB,1); fr = nan(nB,1); validFrac = zeros(nB,1);
    cnt = histcounts(cen, edges+0.5);
    cv = [0; cumsum(double(valid))];
    for b = 1:nB
        i0 = edges(b)+1; i1 = min(N, edges(b+1));
        vsec = (cv(i1+1)-cv(i0))/fs;
        validFrac(b) = vsec / ((i1-i0+1)/fs); % fraction of this bin's duration that was valid
        if vsec > 0; fr(b) = cnt(b)/vsec; end
        t(b) = ((i0+i1)/2-1)/fs;
    end
end

% ======================================================================
function [Ts, F] = fano_curve(cen, valid, N, fs, P)
    Ts = logspace(log10(P.fanoMinWinSec), log10(P.fanoMaxWinSec), P.fanoNWin);
    F = nan(size(Ts));
    cv = [0; cumsum(double(valid))];
    for ti = 1:numel(Ts)
        wN = max(1, round(Ts(ti)*fs)); nW = floor(N/wN);
        if nW < 5; continue; end
        edges = (0:nW)*wN;
        cnt = histcounts(cen, edges+0.5);
        validCount = cv(edges(2:end)+1) - cv(edges(1:end-1)+1);
        keep = validCount(:).' == wN;          % fully-valid windows only
        c = cnt(keep);
        if numel(c) >= 5 && mean(c) > 0; F(ti) = var(c)/mean(c); end
    end
end

% ======================================================================
function [lags, acg] = autocorrelogram(st, cen, cvv, maxLag, binW)
    edges = -maxLag:binW:maxLag; lags = edges(1:end-1)+binW/2; acg = zeros(1,numel(lags));
    n = numel(st); if n < 3; return; end
    diffs = zeros(1, n*10); m = 0;
    for i = 1:n
        j = i+1;
        while j <= n && (st(j)-st(i)) <= maxLag
            % only count pairs whose interval is gap-clean (no censored sample between)
            if (cvv(cen(j)+1) - cvv(cen(i)+1)) == (cen(j) - cen(i))
                m = m+1; if m > numel(diffs); diffs(end*2) = 0; end
                diffs(m) = st(j)-st(i);
            end
            j = j+1;
        end
    end
    diffs = diffs(1:m);
    acg = histcounts([diffs -diffs], edges);
end

% ======================================================================
function [f, pxx] = rate_psd(cen, valid, N, fs, P)
    f = []; pxx = [];
    % Blanking-corrected firing rate over the WHOLE record, then its spectrum.
    % Robust to the per-beat cardiac censoring (which leaves no long contiguous
    % valid run); the corrected rate removes cardiac, so this shows the
    % remaining (respiratory / gastric) rhythms.
    binN = max(1, round(P.psdRateBinSec*fs)); fsr = fs/binN;
    nB = floor(N/binN); if nB < 64; return; end
    edges = (0:nB)*binN;
    cnt = histcounts(cen, edges+0.5);
    cvv = [0; cumsum(double(valid))];
    vsec = (cvv(edges(2:end)+1) - cvv(edges(1:end-1)+1)) / fs;   % valid sec / bin
    rate = cnt(:) ./ max(vsec(:), eps);
    rate(vsec(:) < 0.2*P.psdRateBinSec) = NaN;                  % bin too censored
    rate = fillmissing(rate, 'linear', 'EndValues', 'nearest');
    if numel(rate) < 64 || all(~isfinite(rate)); return; end
    rate = detrend(rate);
    win = min(256, floor(numel(rate)/4)); if win < 16; return; end
    try
        [pxx, f] = pwelch(rate, hamming(win), round(win/2), [], fsr);
    catch
        f = []; pxx = [];
    end
end

% ======================================================================
function B = detect_bursts(st, isims, isClean, validDur, P)
    B = emptyBurst();
    if nnz(isClean) < 20; return; end
    li = log10(isims(isClean & isims>0));
    nb = 50; edges = linspace(min(li), max(li), nb+1); ctr = edges(1:end-1)+diff(edges)/2;
    h = histcounts(li, edges); h = movmean(h, 3);
    [pk, loc] = findpeaks(h);
    if isempty(pk); return; end
    % intra-burst mode = tallest peak at short ISI (< burstMaxThreshMs)
    isShort = ctr(loc) < log10(P.burstMaxThreshMs);
    if ~any(isShort); return; end
    sLoc = loc(isShort); sPk = pk(isShort);
    [~, si] = max(sPk); p1 = sLoc(si);
    % inter-burst mode = tallest peak to the right of the short peak
    rLoc = loc(loc > p1); rPk = pk(loc > p1);
    if isempty(rLoc); return; end
    [~, ri] = max(rPk); p2 = rLoc(ri);
    [valley, vrel] = min(h(p1:p2)); vloc = p1+vrel-1;
    void = 1 - valley/sqrt(h(p1)*h(p2)+eps);
    B.void = void;
    if void < P.burstVoidThresh; return; end                   % not bimodal enough
    thrMs = 10^ctr(vloc); B.thrMs = thrMs;

    inb = (isims < thrMs) & isClean;           % intra-burst intervals (gap-clean only)
    d = diff([false; inb(:); false]);
    bs = find(d==1); be = find(d==-1)-1;       % runs of in-burst intervals
    onsets = []; offsets = []; spk = []; durs = [];
    for r = 1:numel(bs)
        nSpk = (be(r)-bs(r)+1) + 1;            % spikes = intervals + 1
        if nSpk >= P.burstMinSpikes
            i0 = bs(r); i1 = be(r)+1;          % spike indices
            onsets(end+1)=st(i0); offsets(end+1)=st(i1); %#ok<AGROW>
            spk(end+1)=nSpk; durs(end+1)=st(i1)-st(i0);   %#ok<AGROW>
        end
    end
    if isempty(onsets); return; end
    B.hasBursts = true; B.nBursts = numel(onsets); B.onsets = onsets; B.offsets = offsets;
    B.rate_per_min = B.nBursts / (validDur/60);
    B.meanDur_s = mean(durs); B.meanSpikes = mean(spk);
    B.intraRate_hz = mean(spk./max(durs,eps));
    B.fracInBurst = sum(spk) / numel(st);
end

% ======================================================================
function [t, cv2] = rolling_cv2(st, isims, isClean, Tend, winSec)
    edges = 0:winSec:Tend; if edges(end)<Tend; edges(end+1)=Tend; end
    t = edges(1:end-1)+winSec/2; cv2 = nan(1,numel(t));
    tm = st(2:end);                            % time of each ISI (later spike)
    for b = 1:numel(t)
        sel = tm >= edges(b) & tm < edges(b+1) & isClean;
        dd = isims(sel);
        if numel(dd) >= 3; cv2(b) = mean(2*abs(diff(dd))./(dd(1:end-1)+dd(2:end)),'omitnan'); end
    end
end

% ======================================================================
function plot_report(D, k, st, isims, isClean, vpp, valid)
    M = D.metrics(k); fs = D.fs; N = numel(valid); Tend = (N-1)/fs;
    P = pipeline_params();
    fig = figure('Color','w','Name',sprintf('Step 6 report — %s',M.label), ...
        'Position',[60 40 1500 950]);
    tl = tiledlayout(fig,3,4,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('Spike-train report  |  %s%s  |  %d spikes, %.2f spk/s', M.label, ...
        ternary(isempty(M.condition),'',[' | ' M.condition]), M.nSpikes, M.meanRate_hz), 'Interpreter','none');

    % 1 firing rate + burst shading
    ax = nexttile; plot(M.fr_t, M.fr_hz, 'b'); hold on;
    if M.burst.hasBursts
        yl = ylim;
        for i=1:numel(M.burst.onsets)
            patch(ax,[M.burst.onsets(i) M.burst.offsets(i) M.burst.offsets(i) M.burst.onsets(i)], ...
                [yl(1) yl(1) yl(2) yl(2)], [1 0.6 0.2],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
        end
        ylim(yl);
    end
    grid on; xlabel('Time (s)'); ylabel('Rate (spk/s)');
    if M.burst.hasBursts; title('Firing rate (bursts shaded)'); else; title('Firing rate'); end

    % 2 waveform mean+band
    nexttile;
    if ~isempty(M.meanWaveform)
        t = M.wf_t_ms;
        fill([t fliplr(t)], [M.bandLo fliplr(M.bandHi)]*1e6, [0.8 0.85 0.95],'EdgeColor','none','FaceAlpha',0.8); hold on;
        plot(t, M.meanWaveform*1e6,'b','LineWidth',2);
        xlabel('ms'); ylabel('\muV'); title(sprintf('Mean waveform (%d–%d%% band)',P.wfBandPct)); grid on; xlim([t(1) t(end)]);
    else; axis off; end

    % 3 Vpp over time
    nexttile;
    [vt, vmed, vlo, vhi] = binstat(st, vpp, Tend, P.vppBinSec);
    fill([vt fliplr(vt)],[vlo fliplr(vhi)],[0.95 0.85 0.8],'EdgeColor','none','FaceAlpha',0.8); hold on;
    plot(vt, vmed,'r','LineWidth',1.5); grid on; xlabel('Time (s)'); ylabel('V_{pp} (\muV)');
    title('Peak-to-peak over time (median, IQR)');

    % 4 log-ISI histogram (gap-clean) + burst threshold
    nexttile;
    liClean = log10(isims(isClean & isims>0));
    histogram(liClean, 50, 'FaceColor',[0.5 0.6 0.85]); hold on;
    if M.burst.hasBursts; xline(log10(M.burst.thrMs),'r--','LineWidth',1.5); end
    xlabel('log_{10} ISI (ms)'); ylabel('count'); grid on;
    title(sprintf('log-ISI gap-clean (void=%.2f)', M.burst.void));

    % 5 time x ISI heatmap (gap-clean)
    nexttile;
    tmAll = st(2:end); yAll = log10(isims(:));
    tm = tmAll(isClean); y = yAll(isClean);
    tE = 0:P.isiHeatTBinSec:Tend; if numel(tE)<2; tE=[0 Tend]; end
    if numel(y) < 2
        yE = linspace(0,4,40); H = zeros(numel(tE)-1, numel(yE)-1);
    else
        yE = linspace(min(y),max(y),40); H = histcounts2(tm, y, tE, yE);
    end
    imagesc(tE(1:end-1), yE(1:end-1), H.'); set(gca,'YDir','normal');
    xlabel('Time (s)'); ylabel('log_{10} ISI (ms)'); title('Time × ISI density'); colorbar;

    % 6 ISI return map (consecutive gap-clean pairs)
    nexttile;
    both = isClean(1:end-1) & isClean(2:end);
    if any(both)
        xn = isims(1:end-1); xn1 = isims(2:end);
        loglog(xn(both), xn1(both),'.','Color',[0.3 0.3 0.6],'MarkerSize',3);
        xlabel('ISI_n (ms)'); ylabel('ISI_{n+1} (ms)'); title('Return map (clean pairs)'); grid on;
    else; axis off; text(0.2,0.5,'no clean pairs','Units','normalized'); end

    % 7 Fano vs window
    nexttile;
    loglog(M.fano_T, M.fano_F,'o-','Color',[0.2 0.4 0.8]); hold on; yline(1,'k--');
    xlabel('window T (s)'); ylabel('Fano'); grid on;
    title(sprintf('Fano vs window (slope %.2f)', M.fanoSlope));

    % 8 autocorrelogram
    nexttile;
    bar(M.acg_lag*1000, M.acg, 1, 'FaceColor',[0.4 0.5 0.7],'EdgeColor','none');
    xlabel('lag (ms)'); ylabel('count'); title('Autocorrelogram'); grid on;

    % 9 firing-rate power spectrum
    nexttile;
    if ~isempty(M.psd_f)
        plot(M.psd_f, 10*log10(M.psd_p+eps),'b'); hold on;
        for fb = [1 2 6.7]; xline(fb,':','Color',[0.6 0.6 0.6]); end
        xlim([0 P.psdMaxHz]); xlabel('Hz'); ylabel('dB'); title('Firing-rate power spectrum'); grid on;
    else; axis off; text(0.2,0.5,'PSD n/a','Units','normalized'); end

    % 10 amplitude vs preceding ISI (gap-clean)
    nexttile;
    if any(isClean)
        vlater = vpp([false; isClean(:)]);
        semilogx(isims(isClean), vlater,'.','Color',[0.5 0.3 0.3],'MarkerSize',3);
        xlabel('preceding ISI (ms)'); ylabel('V_{pp} (\muV)'); title('Amplitude vs ISI (clean)'); grid on;
    else; axis off; end

    % 11 rolling CV2
    nexttile;
    plot(M.cv2_t, M.cv2_roll,'m','LineWidth',1.2); hold on; yline(1,'k--');
    xlabel('Time (s)'); ylabel('CV2'); ylim([0 2]); grid on; title('Rolling CV2 (regularity)');

    % 12 scalar summary
    nexttile; axis off;
    txt = {
        sprintf('mean rate: %.2f spk/s', M.meanRate_hz)
        sprintf('CV: %.2f   CV2: %.2f   LV: %.2f', M.CV, M.CV2, M.LV)
        sprintf('Fano@%.2gs: %.2f   slope: %.2f', P.fanoCanonSec, M.fanoCanon, M.fanoSlope)
        sprintf('refractory viol: %.3f', M.refracViolFrac)
        sprintf('median Vpp: %.1f uV   FWHM: %.2f ms', M.medianVpp_uv, M.medianFWHM_ms)
        '---'
        sprintf('bursts: %s', burstLine(M.burst))
    };
    if M.burst.hasBursts
        txt{end+1} = sprintf('  thr %.0f ms, dur %.0f ms, %.1f spk/burst', M.burst.thrMs*1, M.burst.meanDur_s*1000, M.burst.meanSpikes);
    end
    text(0.02, 0.98, txt, 'Units','normalized','VerticalAlignment','top','FontName','Consolas','FontSize',9);
    title('Summary');
end

function [t, med, lo, hi] = binstat(st, v, Tend, binSec)
    edges = 0:binSec:Tend; if edges(end)<Tend; edges(end+1)=Tend; end
    t = edges(1:end-1)+binSec/2; med = nan(1,numel(t)); lo = med; hi = med;
    for b=1:numel(t)
        sel = st>=edges(b) & st<edges(b+1);
        if any(sel); q = prctile(v(sel),[25 50 75]); lo(b)=q(1); med(b)=q(2); hi(b)=q(3); end
    end
    ok = isfinite(med); t=t(ok); med=med(ok); lo=lo(ok); hi=hi(ok);
end
function out = ternary(c,a,b); if c; out=a; else; out=b; end; end
