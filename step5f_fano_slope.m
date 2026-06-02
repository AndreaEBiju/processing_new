function R = step5f_fano_slope(D, P, plotMode)
% STEP5F_FANO_SLOPE  Epoch-restricted, surrogate-normalized Fano-vs-window slope.
%
%   R = step5f_fano_slope(D, P, plotMode)
%
% The Fano-factor scaling exponent (slope of log F(T) vs log T) is a measure of
% LONG-RANGE temporal correlation / fractal firing (Teich & Lowen). Two confounds
% are controlled here, per our design discussion:
%   * MOTION blanks (multi-second gaps) corrupt large-window variance -> the
%     analysis is restricted to motion-blank-free EPOCHS, and only windows that
%     fall ENTIRELY inside one epoch are counted (no windows straddle a gap).
%   * CARDIAC blanks (~30 ms) inside an epoch still perturb counts -> a rate-
%     matched Poisson SURROGATE carrying the SAME cardiac blanks is generated;
%     F_norm(T) = F_real(T) / F_surrogate(T) removes the blanking-induced part.
% A robustness variant additionally removes the cardiac-phase-locked spikes
% (nearest-R lag in the T-wave window) before refitting; if the slope survives,
% it is neural, if it collapses it was cardiac.
%
% Needs (per neural channel): D.spikes(k).alignedTimes (or .times), D.validMask,
% D.cardiacBlank (step1a), D.removedSegmentIdx (motion), D.rpeakSamples, D.fs.
%
% Returns R(k): slopeReal, slopeSurr, slopeNorm, slopeNormCardiacRemoved,
% the window vector and the three Fano curves, plus the usable window range.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D,'spikes'); error('step5f:noSpikes','Run step3/step4 first.'); end

    fs  = D.fs; N = size(D.validMask,1);
    chs = D.neuralChannels; nCh = numel(chs); lab = D.channelLabels;

    % ---- params (self-contained) ----
    minEpochSec = 5;            % ignore motion-free epochs shorter than this
    Wmin = 0.02;               % smallest window (s)
    WmaxCap = 5;               % hard cap on window (s)
    nW   = 8;                  % log-spaced window sizes
    nSurr = 20;                % surrogate repetitions
    twaveMs = [20 100];        % cardiac (T-wave) lag window for the removal variant

    rp = round(D.rpeakSamples(:)); rp = rp(rp>=1 & rp<=N);

    % motion-free, non-dead "usable" mask: valid OR only cardiac-blanked
    if isfield(D,'cardiacBlank') && ~isempty(D.cardiacBlank)
        cb = logical(D.cardiacBlank);
    else
        cb = false(N,1);
    end

    R = struct('channel', num2cell(chs));
    fprintf('[step5f] Fano slope (epoch-restricted, surrogate-normalized).\n');

    for c = 1:nCh
        valid = D.validMask(:, c);
        % motion-free epochs: fill SHORT invalid gaps (cardiac blanks + their
        % step2 edge padding + brief artifacts) so they don't break a run; only
        % long motion/dead gaps remain as epoch breaks.
        usable = fill_short_gaps(valid, round(0.25*fs));
        epochs = mask_runs(usable, round(minEpochSec*fs));
        if isempty(epochs)
            warning('step5f:noEpoch','ch %d: no motion-free epoch >= %g s.', chs(c), minEpochSec);
            R(c).slopeNorm = NaN; continue;
        end
        % usable window range: <= a fifth of the longest epoch, capped
        longEp = max(epochs(:,2)-epochs(:,1)+1)/fs;
        Wmax = min(WmaxCap, 0.2*longEp);
        if Wmax <= Wmin*2
            warning('step5f:shortEpoch','ch %d: epochs too short for a slope (Wmax=%.2fs).', chs(c), Wmax);
        end
        W = logspace(log10(Wmin), log10(max(Wmax,Wmin*4)), nW);

        % spike samples for this channel
        if isfield(D.spikes,'alignedTimes') && ~isempty(D.spikes(c).alignedTimes)
            st = round(D.spikes(c).alignedTimes(:)*fs)+1;
        else
            st = round(D.spikes(c).times(:)*fs)+1;
        end
        st = st(st>=1 & st<=N);

        % cardiac-locked removal variant: drop spikes whose nearest-R lag is in T-wave window
        lagms = nearest_lag_ms(st, rp, fs);
        keepCR = ~(lagms >= twaveMs(1) & lagms <= twaveMs(2));

        Freal = fano_curve(st,        epochs, cb, fs, W);
        FrealCR = fano_curve(st(keepCR), epochs, cb, fs, W);

        % rate-matched Poisson surrogate with the SAME cardiac blanks
        Fsurr = surrogate_curve(numel(st),       epochs, cb, fs, W, nSurr);
        FsurrCR = surrogate_curve(nnz(keepCR),   epochs, cb, fs, W, nSurr);

        Fnorm   = Freal   ./ max(Fsurr,   eps);
        FnormCR = FrealCR ./ max(FsurrCR, eps);

        R(c).label = lab{c};
        R(c).W = W; R(c).Freal = Freal; R(c).Fsurr = Fsurr; R(c).Fnorm = Fnorm; R(c).FnormCR = FnormCR;
        R(c).slopeReal = loglog_slope(W, Freal);
        R(c).slopeSurr = loglog_slope(W, Fsurr);
        R(c).slopeNorm = loglog_slope(W, Fnorm);
        R(c).slopeNormCardiacRemoved = loglog_slope(W, FnormCR);
        R(c).Wmax = Wmax; R(c).nEpoch = size(epochs,1);

        fprintf(['   ch %d (%s): slope real=%.2f surr=%.2f NORM=%.2f | cardiac-removed NORM=%.2f ' ...
                 '| %d epochs, Wmax=%.2fs\n'], chs(c), lab{c}, R(c).slopeReal, R(c).slopeSurr, ...
                 R(c).slopeNorm, R(c).slopeNormCardiacRemoved, size(epochs,1), Wmax);

        if plotMode; plot_fano(R(c)); end
    end
end

% ======================================================================
function u = fill_short_gaps(valid, maxGap)
% set short invalid runs (<= maxGap samples) to valid, so cardiac/edge slivers
% don't fragment a motion-free epoch; long motion/dead gaps are left as breaks.
    u = logical(valid(:));
    inv = ~u;
    d = diff([false; inv; false]);
    s = find(d==1); e = find(d==-1)-1;
    for i = 1:numel(s)
        if (e(i)-s(i)+1) <= maxGap; u(s(i):e(i)) = true; end
    end
end

function runs = mask_runs(m, minLen)
% start/stop sample pairs of true-runs in logical m, at least minLen long
    m = m(:); d = diff([false; m; false]);
    s = find(d==1); e = find(d==-1)-1;
    keep = (e-s+1) >= minLen;
    runs = [s(keep) e(keep)];
end

function F = fano_curve(stSamp, epochs, cb, fs, W)
% pooled Fano factor at each window size, counting only windows fully inside one
% motion-free epoch. cb = cardiac-blank logical (full length), unused here for
% counting but the spikes already exclude blanked detections.
    F = nan(1, numel(W));
    for wi = 1:numel(W)
        wsamp = max(1, round(W(wi)*fs));
        counts = [];
        for ei = 1:size(epochs,1)
            a = epochs(ei,1); b = epochs(ei,2);
            nb = floor((b-a+1)/wsamp);
            if nb < 3; continue; end
            edges = a + (0:nb)*wsamp;
            inep = stSamp(stSamp>=edges(1) & stSamp<edges(end));
            ct = histcounts(inep, edges);
            counts = [counts ct]; %#ok<AGROW>
        end
        if numel(counts) >= 5 && mean(counts) > 0
            F(wi) = var(counts) / mean(counts);
        end
    end
end

function F = surrogate_curve(nSpk, epochs, cb, fs, W, nSurr)
% mean Fano curve of a rate-matched homogeneous-Poisson surrogate that carries
% the same cardiac blanks (spikes landing in cb are removed, like the real train).
    fsum = zeros(1, numel(W)); cnt = zeros(1, numel(W));
    % usable (non-cardiac) duration within epochs, for the rate
    usableSamp = 0; fullSamp = 0;
    for ei = 1:size(epochs,1)
        a = epochs(ei,1); b = epochs(ei,2);
        fullSamp = fullSamp + (b-a+1);
        usableSamp = usableSamp + nnz(~cb(a:b));
    end
    if usableSamp <= 0; F = nan(1,numel(W)); return; end
    lambda = nSpk / (usableSamp/fs);            % spikes/s over usable time
    for r = 1:nSurr
        sg = [];
        for ei = 1:size(epochs,1)
            a = epochs(ei,1); b = epochs(ei,2); dur = (b-a+1)/fs;
            ng = round(lambda*dur);           % uniform placement -> Fano~1 null (no toolbox)
            if ng==0; continue; end
            ts = a + sort(floor(rand(ng,1)*(b-a)));   % samples in epoch
            ts = ts(~cb(ts));                          % apply cardiac blank
            sg = [sg; ts]; %#ok<AGROW>
        end
        Fr = fano_curve(sg, epochs, cb, fs, W);
        ok = isfinite(Fr); fsum(ok) = fsum(ok)+Fr(ok); cnt(ok)=cnt(ok)+1;
    end
    F = fsum ./ max(cnt,1); F(cnt==0) = NaN;
end

function s = loglog_slope(W, F)
    ok = isfinite(F) & F>0 & isfinite(W) & W>0;
    if nnz(ok) < 3; s = NaN; return; end
    p = polyfit(log10(W(ok)), log10(F(ok)), 1); s = p(1);
end

function lagms = nearest_lag_ms(stSamp, rp, fs)
    lagms = nan(size(stSamp));
    if isempty(rp); return; end
    st = stSamp/fs; ev = sort(rp(:))/fs;
    k = interp1(ev, 1:numel(ev), st, 'nearest', 'extrap');
    k = min(max(round(k),1), numel(ev));
    lag = st - ev(k);
    for d = [-1 1]
        j = min(max(k+d,1), numel(ev)); alt = st - ev(j);
        better = abs(alt) < abs(lag); lag(better) = alt(better);
    end
    lagms = lag*1000;
end

function plot_fano(Rc)
    figure('Color','w','Name',sprintf('Step 5f Fano slope — %s', Rc.label), ...
        'Position',[150 150 760 560]);
    loglog(Rc.W, Rc.Freal, '-o', 'LineWidth',1.6, 'Color',[0.1 0.3 0.8], 'DisplayName','real'); hold on;
    loglog(Rc.W, Rc.Fsurr, '--s', 'LineWidth',1.3, 'Color',[0.8 0.3 0.3], 'DisplayName','surrogate (blank null)');
    loglog(Rc.W, Rc.Fnorm, '-^', 'LineWidth',1.8, 'Color',[0 0.5 0], 'DisplayName','normalized (real/surr)');
    loglog(Rc.W, Rc.FnormCR, ':d', 'LineWidth',1.5, 'Color',[0.5 0 0.5], 'DisplayName','normalized, cardiac-removed');
    yline(1,'k:','HandleVisibility','off');
    grid on; xlabel('counting window T (s)'); ylabel('Fano factor');
    title(sprintf(['%s — Fano vs window  |  slope_{norm}=%.2f  (cardiac-removed %.2f)\n' ...
        '%d epochs, Wmax=%.2fs'], Rc.label, Rc.slopeNorm, Rc.slopeNormCardiacRemoved, ...
        Rc.nEpoch, Rc.Wmax), 'Interpreter','tex');
    legend('Location','northwest');
end
