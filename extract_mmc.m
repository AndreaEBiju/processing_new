function mmc = extract_mmc(blankFile, hrbrFile, opts)
% EXTRACT_MMC  Gastric/intestinal motor ("MMC-like") spike-burst activity from a
% 3-channel time series stored in a _blankmotion.mat file. Cardiac artifacts are
% removed with R-peaks detected from the paired _HRBR ECG (blank-before-filter),
% the signal is bandpassed to the spike/burst band, bursts are detected by
% adaptive median/sigma thresholding (one event per burst), and a per-file cache
% <base>_mmc.mat is written next to the input file.
%
%   mmc = extract_mmc(blankFile, hrbrFile, opts)
%
% REQUIRED opts:
%   .gastricVar   name of the N x 3 array in blankFile
%   .ecgVar       name of the ECG trace in hrbrFile
% OPTIONAL opts (defaults in brackets):
%   .fsVar        fs variable name in blankFile         [auto: fs/Fs/sampleRate/samplerate]
%   .ecgFsVar     fs variable name in hrbrFile           [auto, else = gastric fs]
%   .band         bandpass band, Hz                      [2 50]
%   .W .S         rate window / step, s                  [10  1]
%   .k            threshold multiplier (k*sigma)         [4]
%   .sigmaWin     moving robust-sigma window, s          [1]
%   .refractory   merge gap that groups spikes->burst, s [0.5]
%   .cardiacBlankMs  half-width blanked around each R, ms[25]
%   .minValidFrac min non-NaN fraction to score a window [0.5]
%   .delayW .delayStep .delayMaxLag  cross-channel delay [30 5 3] s
%   .save         write <base>_mmc.mat                   [true]
%   .rpeakTimes   (override) supplied R-peak times (s); skips ECG detection []

    if nargin < 3; opts = struct(); end
    req = {'gastricVar','ecgVar'};
    for r = req; assert(isfield(opts,r{1}) && ~isempty(opts.(r{1})), ...
            'extract_mmc: opts.%s is required.', r{1}); end
    g = @(f,d) getdef(opts,f,d);
    band   = g('band',[2 50]);   W = g('W',10);  S = g('S',1);
    k      = g('k',4);           sigmaWin = g('sigmaWin',1);
    refr   = g('refractory',0.5);cardMs = g('cardiacBlankMs',25);
    minVF  = g('minValidFrac',0.5);
    dW = g('delayW',30); dS = g('delayStep',5); dMax = g('delayMaxLag',3);
    doSave = g('save',true);

    % ---- load gastric N x 3 + fs ----
    Sg = load(blankFile, opts.gastricVar);
    G = Sg.(opts.gastricVar); G = double(G);
    if size(G,1) < size(G,2); G = G.'; end           % force N x 3
    assert(size(G,2)>=3, 'extract_mmc: %s is not >=3 channels.', opts.gastricVar);
    G = G(:,1:3);
    fs = read_fs(blankFile, getdef(opts,'fsVar',''));
    N = size(G,1); t = (0:N-1).'/fs;

    % ---- R-peak times from HRBR ECG (or supplied) ----
    if isfield(opts,'rpeakTimes') && ~isempty(opts.rpeakTimes)
        rT = opts.rpeakTimes(:);
    else
        Se = load(hrbrFile, opts.ecgVar); ecg = double(Se.(opts.ecgVar)(:));
        fsE = read_fs(hrbrFile, getdef(opts,'ecgFsVar',''), fs);   % fall back to gastric fs
        rT = detect_rpeaks(ecg, fsE);
    end
    rIdx = unique(min(max(round(rT*fs),1),N));        % R-peaks in gastric samples
    half = max(1, round(cardMs/1000*fs));

    % ---- per-channel processing ----
    cond = nan(N,3);                                  % conditioned signal (full rate)
    evIdx = cell(1,3);                                % burst-peak sample indices
    [bb,aa] = butter(4, [band(1) min(band(2),0.45*fs)]/(fs/2), 'bandpass');
    for ch = 1:3
        x = G(:,ch);
        % cardiac blank-before-filter: NaN around R, interp, filtfilt, restore NaN
        bl = false(N,1);
        for ri = rIdx(:)'
            lo = max(1,ri-half); hi = min(N,ri+half); bl(lo:hi) = true;
        end
        xb = x; xb(bl) = NaN;
        xf = fillmissing(xb,'linear','EndValues','nearest');
        xf(~isfinite(xf)) = 0;
        y = filtfilt(bb,aa,xf);
        y(bl) = NaN;                                  % keep cardiac gaps as NaN
        cond(:,ch) = y;
        % adaptive median/sigma detection -> burst peaks
        evIdx{ch} = detect_bursts(y, fs, k, sigmaWin, refr);
    end

    % ---- sliding-window burst rate + per-window peak amplitude ----
    centers = (W/2 : S : (t(end)-W/2)).';
    M = numel(centers);
    rate = nan(M,3); peakAmp = nan(M,3);
    for ch = 1:3
        valid = ~isnan(cond(:,ch));
        cumv = [0; cumsum(double(valid))];
        pk = evIdx{ch}; pkt = pk/fs; pka = abs(cond(pk,ch));
        for w = 1:M
            lo = max(1, floor((centers(w)-W/2)*fs)+1);
            hi = min(N, floor((centers(w)+W/2)*fs));
            validDur = (cumv(hi+1)-cumv(lo))/fs;
            if validDur < minVF*W; continue; end
            inw = pkt>=centers(w)-W/2 & pkt<centers(w)+W/2;
            rate(w,ch) = sum(inw)/validDur;
            if any(inw); peakAmp(w,ch) = mean(pka(inw)); end
        end
    end
    avgMMCRate = zeros(1,3);
    for ch = 1:3
        valid = ~isnan(cond(:,ch));
        avgMMCRate(ch) = numel(evIdx{ch}) / max(sum(valid)/fs, eps);
    end

    % ---- inter-channel propagation delay (pairs 1-2,1-3,2-3) on the rate series
    [delay, delay_t] = xchan_delay(rate, S, dW, dS, dMax);

    % ---- store conditioned signal + event boolean at FULL fs ----
    evS = false(N,3);
    for ch = 1:3
        evS(evIdx{ch},ch) = true;
    end

    % ---- QC products (kept small; make plot_mmc self-contained) ----
    qc = struct();
    qc.srcFile = blankFile; qc.gastricVar = opts.gastricVar;
    qc.rpeakT = rIdx/fs;
    if numel(rIdx) > 2; qc.meanHR = 60/median(diff(rIdx)/fs); else; qc.meanHR = NaN; end
    qc.pctBlanked  = mean(isnan(cond),1)*100;          % 1x3, % time cardiac-blanked
    qc.rateNanFrac = mean(isnan(rate),1)*100;          % 1x3, % NaN rate windows
    qc.nBursts     = cellfun(@numel, evIdx);           % 1x3
    Lh = round(0.1*fs); seg = (-Lh:Lh).'; qc.periR_t = seg/fs;     % peri-R average +/-100ms
    pr = zeros(numel(seg),3); pc = zeros(numel(seg),3); cnt = 0;
    for ri = rIdx(:)'
        if ri-Lh<1 || ri+Lh>N; continue; end
        w = ri-Lh:ri+Lh; pr = pr + G(w,:);
        cc = cond(w,:); cc(isnan(cc)) = 0; pc = pc + cc; cnt = cnt + 1;
    end
    if cnt > 0; qc.periR_raw = pr/cnt; qc.periR_cond = pc/cnt;
    else; qc.periR_raw = nan(numel(seg),3); qc.periR_cond = nan(numel(seg),3); end
    try                                                 % Welch PSD raw vs conditioned
        Lp = 2^floor(log2(min(N, round(4*fs)))); wp = hann(Lp);
        [Pr,fp] = pwelch(fillmissing(G,'constant',0), wp, Lp/2, Lp, fs);
        Pc = pwelch(fillmissing(double(cond),'constant',0), wp, Lp/2, Lp, fs);
        s = fp <= 100; qc.psd_f = fp(s); qc.psd_raw = Pr(s,:); qc.psd_cond = Pc(s,:);
    catch
        qc.psd_f = []; qc.psd_raw = []; qc.psd_cond = [];
    end

    mmc = struct();
    mmc.fs = fs; mmc.t = t;
    mmc.signal = single(cond);    % full-rate conditioned (cardiac-blanked, bandpassed)
    mmc.eventSeries = evS;        % full-rate logical, 1 at each burst peak
    mmc.rate_t = centers; mmc.rate = rate; mmc.peakAmp = peakAmp;
    mmc.avgMMCRate = avgMMCRate;
    mmc.delay_t = delay_t; mmc.delay = delay;   % cols = pairs [1-2 1-3 2-3]
    mmc.pairs = {'1-2','1-3','2-3'};
    mmc.params = struct('band',band,'W',W,'S',S,'k',k,'sigmaWin',sigmaWin, ...
        'refractory',refr,'cardiacBlankMs',cardMs,'minValidFrac',minVF, ...
        'delayW',dW,'delayStep',dS,'delayMaxLag',dMax, ...
        'gastricVar',opts.gastricVar,'ecgVar',opts.ecgVar,'nRpeaks',numel(rIdx));
    mmc.qc = qc;

    if doSave
        [d,b] = fileparts(blankFile);
        save(fullfile(d,[b '_mmc.mat']), 'mmc');
        fprintf('  [mmc] %s  rate(ch)= %s Hz\n', b, num2str(avgMMCRate,'%.3f '));
    end
end

% ======================================================================
function v = getdef(s,f,d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end

function fs = read_fs(file, name, fallback)
    if nargin<3; fallback=[]; end
    if ~isempty(name)
        S = load(file,name); fs = double(S.(name)(1)); return;
    end
    cands = {'fs','Fs','sampleRate','samplerate','SampleRate','samplingRate'};
    info = who('-file',file);
    for c = cands
        if any(strcmp(info,c{1})); S=load(file,c{1}); fs=double(S.(c{1})(1)); return; end
    end
    assert(~isempty(fallback), 'read_fs: no fs variable found in %s; set opts.fsVar.', file);
    fs = fallback;
end

function rT = detect_rpeaks(ecg, fsE)
% Simple robust R-peak detector: bandpass 5-30 Hz, rectify, MAD threshold,
% min RR refractory. Returns peak times (s).
    ecg = ecg(:); ecg(~isfinite(ecg)) = 0;
    hi = min(30, 0.45*fsE);
    [b,a] = butter(2, [5 hi]/(fsE/2), 'bandpass');
    z = filtfilt(b,a,ecg); z = abs(z);
    thr = median(z) + 5*1.4826*median(abs(z-median(z)));
    minRR = round(0.08*fsE);                          % up to ~12 Hz HR
    [~,locs] = findpeaks(z,'MinPeakHeight',thr,'MinPeakDistance',max(1,minRR));
    rT = locs(:)/fsE;
end

function pk = detect_bursts(y, fs, k, sigmaWin, refr)
% Adaptive median/sigma threshold (robust), group supra-threshold crossings into
% bursts using a refractory merge gap, return one peak (max |y|) per burst.
    win = max(3, round(sigmaWin*fs));
    med = movmedian(y, win, 'omitnan');
    sig = movmedian(abs(y-med), win, 'omitnan') / 0.6745;   % robust moving sigma (MAD)
    sig(sig==0 | ~isfinite(sig)) = median(sig(isfinite(sig)&sig>0),'omitnan');
    above = (abs(y-med) > k*sig) & isfinite(y);
    idx = find(above);
    if isempty(idx); pk = []; return; end
    gaps = [Inf; diff(idx)];
    burstId = cumsum(gaps > refr*fs);
    nb = burstId(end); pk = zeros(nb,1);
    for b = 1:nb
        gi = idx(burstId==b);
        [~,mi] = max(abs(y(gi)-med(gi)));
        pk(b) = gi(mi);
    end
end

function [delay, delay_t] = xchan_delay(rate, S, dW, dS, dMax)
% Sliding-window lag (s) of peak cross-correlation between channel pairs, on the
% rate series. Positive lag = second channel leads. cols: [1-2 1-3 2-3].
    pairs = [1 2; 1 3; 2 3];
    wlen = max(4, round(dW/S)); step = max(1, round(dS/S)); mlag = round(dMax/S);
    M = size(rate,1);
    starts = 1:step:max(1,M-wlen+1);
    delay = nan(numel(starts),3); delay_t = nan(numel(starts),1);
    for s = 1:numel(starts)
        lo = starts(s); hi = min(M, lo+wlen-1);
        delay_t(s) = (lo+hi)/2 * S;
        for p = 1:3
            a = rate(lo:hi,pairs(p,1)); b = rate(lo:hi,pairs(p,2));
            ok = isfinite(a)&isfinite(b);
            if nnz(ok) < 4 || std(a(ok))==0 || std(b(ok))==0; continue; end
            a=a-mean(a(ok)); b=b-mean(b(ok)); a(~ok)=0; b(~ok)=0;
            [c,lags] = xcorr(a, b, mlag, 'coeff');
            [~,mi] = max(c); delay(s,p) = lags(mi)*S;
        end
    end
end
