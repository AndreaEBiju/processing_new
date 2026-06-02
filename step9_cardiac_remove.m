function D = step9_cardiac_remove(D, P, plotMode)
% STEP9_CARDIAC_REMOVE  Detect and remove cardiac-locked (heartbeat-entrained) spikes.
%
%   D = step9_cardiac_remove(D, P, plotMode)
%
% Builds a peri-R-wave histogram of the accepted spikes (offset to nearest
% R-peak). If a channel is genuinely cardiac-coupled (peak z >= cardiacZThresh
% above the off-peak baseline), it:
%   * defines an exclusion window from the WIDTH of the cardiac peak
%     (half-height, capped at +/- cardiacWinMaxMs) -- data-driven, not a fixed
%     guess;
%   * removes the accepted spikes whose R-offset falls in that window;
%   * CENSORS that time (adds the per-beat windows to D.validMask) so the
%     firing-rate denominator shrinks accordingly and the rate stays unbiased.
% Channels with no significant coupling are left untouched.
%
% Set P.cardiacRemove = false to tag/measure only (no removal).
%
% Requires D.spikes(k) accepted fields (step 4), D.rpeakSamples, D.validMask.
% Updates D.spikes(k) (cleaned train, rate, waveform) and D.validMask; stores
% D.cardiacTag(k): lag, hist, baseline, z, winLo/Hi (ms), nRemoved, fracRemoved.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D,'spikes') || ~isfield(D.spikes,'alignedTimes')
        error('step9_cardiac_remove:missing','Run step4_waveforms first.');
    end

    fs = D.fs; N = size(D.filtered,1);
    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    rt = (round(D.rpeakSamples(:))-1)/fs;          % R-peak times (s)
    rs = round(D.rpeakSamples(:));
    edges = (-P.cardiacLagMaxMs):P.cardiacLagBinMs:(P.cardiacLagMaxMs);
    L = edges(1:end-1) + P.cardiacLagBinMs/2;      % bin centres (ms)

    D.cardiacTag = repmat(struct('lag',L,'hist',[],'baseline',NaN,'z',NaN, ...
        'coupled',false,'winLo',NaN,'winHi',NaN,'nRemoved',0,'fracRemoved',0), 1, nCh);

    fprintf('[step9] Cardiac-locked spike removal (peri-R-wave), remove=%d.\n', P.cardiacRemove);

    for k = 1:nCh
        st = D.spikes(k).alignedTimes(:);
        if numel(st) < 10 || isempty(rt)
            D.cardiacTag(k).hist = zeros(1,numel(L)); continue;
        end
        % offset of each spike to its nearest R-peak (ms)
        nearestRt = interp1(rt, rt, st, 'nearest', 'extrap');
        offMs = (st - nearestRt) * 1000;
        h = histcounts(offMs, edges);

        far = abs(L) > P.cardiacBaseLagMs;
        baseline = mean(h(far)); sdv = std(h(far));
        [pkc, pki] = max(h);
        z = (pkc - baseline) / max(sdv, eps);

        T = D.cardiacTag(k);
        T.hist = h; T.baseline = baseline; T.z = z;
        T.coupled = z >= P.cardiacZThresh;

        keep = true(size(st));
        % --- determine the exclusion window ---
        if P.cardiacUniform
            % uniform fixed window on every channel (comparable across conditions)
            winLo = -P.cardiacRemoveWinMs; winHi = P.cardiacRemoveWinMs;
        elseif T.coupled
            % data-driven window from the half-height width of the cardiac peak
            halfThr = baseline + 0.5*(pkc - baseline);
            lo = pki; while lo>1      && h(lo-1) > halfThr; lo = lo-1; end
            hi = pki; while hi<numel(h)&& h(hi+1) > halfThr; hi = hi+1; end
            winLo = max(L(lo)-P.cardiacLagBinMs/2, -P.cardiacWinMaxMs);
            winHi = min(L(hi)+P.cardiacLagBinMs/2,  P.cardiacWinMaxMs);
        else
            winLo = NaN; winHi = NaN;
        end
        T.winLo = winLo; T.winHi = winHi;

        if isfinite(winLo)
            inWin = offMs >= winLo & offMs <= winHi;
            T.nRemoved = nnz(inWin); T.fracRemoved = nnz(inWin)/numel(st);
            if P.cardiacRemove
                keep = ~inWin;
                % censor the per-beat window time in the validity mask
                v = D.validMask(:,k);
                aS = round(winLo/1000*fs); bS = round(winHi/1000*fs);
                for j = 1:numel(rs)
                    a = max(1, rs(j)+aS); b = min(N, rs(j)+bS);
                    if b>=a; v(a:b) = false; end
                end
                D.validMask(:,k) = v;
            end
        end

        % apply removal to the accepted train + recompute
        if P.cardiacRemove && any(~keep)
            D.spikes(k) = subset_spikes(D.spikes(k), keep);
        end
        validSec = nnz(D.validMask(:,k))/fs;
        D.spikes(k).validSec = validSec;
        D.spikes(k).nSpikes  = numel(D.spikes(k).alignedTimes);
        D.spikes(k).rate_hz  = D.spikes(k).nSpikes / max(validSec, eps);
        D.cardiacTag(k) = T;

        if P.cardiacRemove && isfinite(T.winLo)
            fprintf(['        ch %d (%s): removed %d spikes (%.0f%%) in [%.0f, %.0f] ms ' ...
                     '(coupling z=%.1f) | rate -> %.2f /s\n'], ch(k), lab{k}, ...
                T.nRemoved, 100*T.fracRemoved, T.winLo, T.winHi, z, D.spikes(k).rate_hz);
        else
            fprintf('        ch %d (%s): no removal (coupling z=%.1f)\n', ch(k), lab{k}, z);
        end

        if plotMode; plot_periR(D, k, edges); end
    end
end

% ======================================================================
function s = subset_spikes(s, keep)
    flds = {'alignedCenters','alignedTimes','Vpp_uv','width_ms'};
    for i = 1:numel(flds)
        if isfield(s, flds{i}) && numel(s.(flds{i})) == numel(keep)
            s.(flds{i}) = s.(flds{i})(keep);
        end
    end
    if isfield(s,'waveforms') && size(s.waveforms,1) == numel(keep)
        s.waveforms = s.waveforms(keep,:);
    end
    if isfield(s,'waveforms') && ~isempty(s.waveforms)
        s.meanWaveform = mean(s.waveforms,1,'omitnan');
        s.stdWaveform  = std(s.waveforms,0,1,'omitnan');
    end
end

% ======================================================================
function plot_periR(D, k, edges)
    T = D.cardiacTag(k); lab = D.channelLabels{k};
    figure('Color','w','Name',sprintf('Step 9 — %s',lab),'Position',[200 200 720 420]);
    bar(T.lag, T.hist, 1, 'FaceColor',[0.6 0.6 0.85],'EdgeColor','none'); hold on;
    yl = ylim;
    if isfinite(T.winLo)
        patch([T.winLo T.winHi T.winHi T.winLo],[yl(1) yl(1) yl(2) yl(2)], ...
            [1 0.5 0.5],'FaceAlpha',0.25,'EdgeColor','none');
    end
    yline(T.baseline,'k--');
    grid on; xlabel('spike offset from R-peak (ms)'); ylabel('count');
    if isfinite(T.winLo)
        ttl = sprintf('%s: removed %.0f%% in [%.0f,%.0f] ms (coupling z=%.1f)', ...
            lab, 100*T.fracRemoved, T.winLo, T.winHi, T.z);
    else
        ttl = sprintf('%s: no removal (z=%.1f)', lab, T.z);
    end
    title(ttl, 'Interpreter','none'); ylim(yl);
end
