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
%   .gastricCols  3 column indices of the stomach channels within the combined
%                 (nerve+stomach) array in blankFile, e.g. [5 6 7]
%   cardiac source -- ONE of:
%     .rpeakVar   name of a precomputed R-peak vector in hrbrFile (preferred), OR
%     .ecgVar     name of a raw ECG trace in hrbrFile (R-peaks detected here), OR
%     .rpeakTimes R-peak times (s) passed directly
% OPTIONAL opts (defaults in brackets):
%   .rpeakUnits   'seconds' | 'samples' for .rpeakVar     [seconds]
%   .rpeakFs      fs for 'samples' R-peaks                 [HRBR fs / gastric fs]
%   .rpeakFile    file holding .rpeakVar                   [hrbrFile]
%   .dataVar      name of the combined array in blankFile  [auto: largest 2-D numeric]
%   .fsVar        fs variable name in blankFile            [auto: fs/Fs/sampleRate/samplerate]
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
    assert(isfield(opts,'gastricCols') && ~isempty(opts.gastricCols), ...
        'extract_mmc: opts.gastricCols is required.');
    assert((isfield(opts,'rpeakVar')  && ~isempty(opts.rpeakVar))  || ...
           (isfield(opts,'ecgVar')    && ~isempty(opts.ecgVar))    || ...
           (isfield(opts,'rpeakTimes')&& ~isempty(opts.rpeakTimes)), ...
        'extract_mmc: give a cardiac source: opts.rpeakVar OR opts.ecgVar OR opts.rpeakTimes.');
    g = @(f,d) getdef(opts,f,d);
    band   = g('band',[2 50]);   W = g('W',10);  S = g('S',1);
    k      = g('k',3);           sigmaWin = g('sigmaWin',30);   % LONG (>> burst duration)
    spikeRefr = g('spikeRefractory', 0.05);                 % individual muscle firings
    burstRefr = g('burstRefractory', g('refractory', 0.5)); % grouped burst episodes
    cardMs = g('cardiacBlankMs',25);
    minVF  = g('minValidFrac',0.5);
    dW = g('delayW',30); dS = g('delayStep',5); dMax = g('delayMaxLag',3);
    doSave = g('save',true);

    % ---- load combined (nerve+stomach) array, select the 3 stomach columns ----
    cols = opts.gastricCols(:).';
    assert(numel(cols)==3, 'extract_mmc: opts.gastricCols must list 3 stomach channel indices.');
    dataVar = getdef(opts,'dataVar','');
    if isempty(dataVar); dataVar = biggest_matrix_var(blankFile); end
    Sg = load(blankFile, dataVar);
    D = double(Sg.(dataVar));
    if size(D,1) < size(D,2); D = D.'; end           % N x nChannels (time down columns)
    assert(max(cols) <= size(D,2), 'extract_mmc: gastricCols [%s] exceed #channels (%d) in %s.', ...
        num2str(cols), size(D,2), dataVar);
    G = D(:, cols);
    fs = read_fs(blankFile, getdef(opts,'fsVar',''));
    N = size(G,1); t = (0:N-1).'/fs;

    % ---- R-peak times: precomputed locations preferred; else detect from ECG ----
    if isfield(opts,'rpeakTimes') && ~isempty(opts.rpeakTimes)
        rT = opts.rpeakTimes(:);
    elseif isfield(opts,'rpeakVar') && ~isempty(opts.rpeakVar)
        rf = getdef(opts,'rpeakFile', hrbrFile);
        Sp = load(rf, opts.rpeakVar); pk = double(Sp.(opts.rpeakVar)(:));
        if startsWith(lower(getdef(opts,'rpeakUnits','seconds')), 'sample')
            pfs = getdef(opts,'rpeakFs', fs);       % default: indices are into the main (gastric) recording
            rT = pk / pfs;                          % sample indices -> seconds
        else
            rT = pk;                                % already in seconds
        end
    else
        Se = load(hrbrFile, opts.ecgVar); ecg = double(Se.(opts.ecgVar)(:));
        fsE = read_fs(hrbrFile, getdef(opts,'ecgFsVar',''), fs);
        rT = detect_rpeaks(ecg, fsE);
    end
    rIdx = unique(min(max(round(rT*fs),1),N));        % R-peaks in gastric samples
    half = max(1, round(cardMs/1000*fs));

    % ---- per-channel processing ----
    cond = nan(N,3);                                  % conditioned signal (full rate)
    spkIdx = cell(1,3); brsIdx = cell(1,3);           % event peaks: firings / grouped bursts
    [zf,pf,kf] = butter(4, [band(1) min(band(2),0.45*fs)]/(fs/2), 'bandpass');
    [sos,gd]   = zp2sos(zf,pf,kf);   % SOS form: numerically stable at very low normalized cutoffs (high fs)
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
        y = filtfilt(sos,gd,xf);
        y(bl) = NaN;                                  % keep cardiac gaps as NaN
        cond(:,ch) = y;
        spkIdx{ch} = detect_bursts(y, fs, k, sigmaWin, spikeRefr);  % individual firings
        brsIdx{ch} = detect_bursts(y, fs, k, sigmaWin, burstRefr);  % grouped episodes
    end

    % ---- sliding-window rate + per-window peak amplitude, for BOTH levels ----
    centers = (W/2 : S : (t(end)-W/2)).';
    [fRate,fAmp,fAvg] = event_rate(cond, spkIdx, fs, N, centers, W, minVF); % individual firings
    [bRate,bAmp,bAvg] = event_rate(cond, brsIdx, fs, N, centers, W, minVF); % grouped bursts

    % ---- inter-channel propagation delay on the (denser) firing-rate series ----
    [delay, delay_t] = xchan_delay(fRate, S, dW, dS, dMax);

    % ---- QC products (kept small; make plot_mmc self-contained) ----
    qc = struct();
    qc.srcFile = blankFile; qc.dataVar = dataVar; qc.gastricCols = cols;
    qc.rpeakT = rIdx/fs;
    if numel(rIdx) > 2; qc.meanHR = 60/median(diff(rIdx)/fs); else; qc.meanHR = NaN; end
    qc.pctBlanked  = mean(isnan(cond),1)*100;          % 1x3, % time cardiac-blanked
    qc.rateNanFrac = mean(isnan(bRate),1)*100;         % 1x3, % NaN rate windows
    qc.nFirings    = cellfun(@numel, spkIdx);          % 1x3, individual firings
    qc.nBursts     = cellfun(@numel, brsIdx);          % 1x3, grouped bursts
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
    mmc.rate_t = centers;
    mmc.firing = struct('events', ev_bool(spkIdx,N), 'rate',fRate,'peakAmp',fAmp, ...
                        'avgRate',fAvg, 'refractory',spikeRefr);   % individual muscle firings
    mmc.burst  = struct('events', ev_bool(brsIdx,N), 'rate',bRate,'peakAmp',bAmp, ...
                        'avgRate',bAvg, 'refractory',burstRefr);   % grouped burst episodes
    mmc.delay_t = delay_t; mmc.delay = delay;   % cols = pairs [1-2 1-3 2-3], on firing rate
    mmc.pairs = {'1-2','1-3','2-3'};
    mmc.params = struct('band',band,'W',W,'S',S,'k',k,'sigmaWin',sigmaWin, ...
        'spikeRefractory',spikeRefr,'burstRefractory',burstRefr,'cardiacBlankMs',cardMs, ...
        'minValidFrac',minVF,'delayW',dW,'delayStep',dS,'delayMaxLag',dMax, ...
        'dataVar',dataVar,'gastricCols',cols, ...
        'ecgVar',getdef(opts,'ecgVar',''),'rpeakVar',getdef(opts,'rpeakVar',''),'nRpeaks',numel(rIdx));
    mmc.qc = qc;

    if doSave
        [d,b] = fileparts(blankFile);
        save(fullfile(d,[b '_mmc.mat']), 'mmc');
        fprintf('  [mmc] %s  firings/s= %s | bursts/s= %s\n', b, ...
            num2str(fAvg,'%.2f '), num2str(bAvg,'%.2f '));
    end
end

% ======================================================================
function v = getdef(s,f,d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end

function name = biggest_matrix_var(file)
% Pick the largest 2-D numeric variable in the .mat file (the data matrix).
    w = whos('-file', file); best = ''; bestN = 0;
    for i = 1:numel(w)
        if numel(w(i).size)==2 && min(w(i).size)>=1 && ...
           any(strcmp(w(i).class,{'double','single','int16','int32','uint16'}))
            n = prod(w(i).size);
            if n > bestN; bestN = n; best = w(i).name; end
        end
    end
    assert(~isempty(best), 'biggest_matrix_var: no numeric matrix found in %s; set opts.dataVar.', file);
    name = best;
end

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
% Detect muscle-firing bursts with a moving median/MAD adaptive threshold at
% k*sigma. The sigma window (sigmaWin) must be LONG -- many burst-and-gap cycles
% -- so the MAD tracks the noise floor (pulled down by the quiet inter-burst
% stretches) and adapts to slow noise/gain drift, WITHOUT being inflated by the
% amplitude of the burst it currently sits inside (which a short window is, so a
% short window rides up with the bursts and misses them). Supra-threshold
% crossings are grouped into bursts by the refractory gap; one peak per burst.
    win = max(round(3*fs), round(sigmaWin*fs));             % long window (>> burst duration)
    med = movmedian(y, win, 'omitnan');
    sig = movmedian(abs(y-med), win, 'omitnan') / 0.6745;   % moving MAD (adaptive, not burst-inflated)
    sig(~isfinite(sig) | sig==0) = median(sig(isfinite(sig)&sig>0),'omitnan');
    above = (abs(y-med) > k*sig) & isfinite(y);
    idx = find(above);
    if isempty(idx); pk = []; return; end
    % gaps measured in VALID (non-blanked) time: a cardiac blank between two
    % supra-threshold events is NOT counted as an inter-burst gap, so bursts are
    % bridged across the blanks rather than split/terminated by them.
    cumValid = cumsum(isfinite(y));
    vgap = [Inf; (cumValid(idx(2:end)) - cumValid(idx(1:end-1)))/fs];
    burstId = cumsum(vgap > refr);
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

function [rate, peakAmp, avgRate] = event_rate(cond, evIdx, fs, N, centers, W, minVF)
% Sliding-window event rate and mean per-window peak amplitude, per channel.
% Rate = events / VALID duration (cardiac-blanked samples are excluded from the
% denominator, so blanks are never mistaken for quiescence). Windows whose valid
% fraction < minVF are left NaN. avgRate = whole-record events / valid seconds.
    M = numel(centers); rate = nan(M,3); peakAmp = nan(M,3); avgRate = zeros(1,3);
    for ch = 1:3
        valid = ~isnan(cond(:,ch)); cumv = [0; cumsum(double(valid))];
        pk = evIdx{ch}; pkt = pk/fs; pka = abs(cond(pk,ch));
        for w = 1:M
            lo = max(1, floor((centers(w)-W/2)*fs)+1);
            hi = min(N, floor((centers(w)+W/2)*fs));
            vd = (cumv(hi+1)-cumv(lo))/fs;
            if vd < minVF*W; continue; end
            inw = pkt>=centers(w)-W/2 & pkt<centers(w)+W/2;
            rate(w,ch) = sum(inw)/vd;
            if any(inw); peakAmp(w,ch) = mean(pka(inw)); end
        end
        avgRate(ch) = numel(evIdx{ch}) / max(sum(valid)/fs, eps);
    end
end

function ev = ev_bool(evIdx, N)
% Full-rate logical event series, 1 at each event peak sample. cols = channels.
    ev = false(N,3);
    for ch = 1:3; ev(evIdx{ch},ch) = true; end
end
