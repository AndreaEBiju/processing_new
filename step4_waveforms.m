function D = step4_waveforms(D, P, plotMode)
% STEP4_WAVEFORMS  Extract, align, and characterise detected-spike waveforms.
%
%   D = step4_waveforms(D, P, plotMode)
%
% For each event from step 3, re-aligns to the local extremum (within
% +/- P.wfAlignSearchMs), extracts a [-P.wfPreMs, +P.wfPostMs] window, and
% re-enforces the refractory period on the aligned peaks (keeping the
% strongest) so a multiphasic event yields a single waveform. Aligning on a
% consistent feature is essential so identical spikes do not split across
% feature space later.
%
% Per spike it computes peak-to-peak amplitude (Vpp) and full-width at
% half-amplitude (FWHM); these feed fiber-character description and the
% sorting features in step 5.
%
% Requires D.filtered (cardiac-cleaned) and D.spikes (step 3).
%
% Augments each D.spikes(k) with:
%   .wf_t_ms        waveform time axis (ms, 0 = aligned peak)
%   .waveforms      nSpikes x winLen aligned waveforms (V)
%   .alignedCenters .alignedTimes
%   .Vpp_uv         per-spike peak-to-peak amplitude (uV)
%   .width_ms       per-spike FWHM (ms)
%   .meanWaveform .stdWaveform
%
% plotMode (single-file): per channel, the aligned waveform overlay with
% mean +/- std, and an amplitude-vs-width scatter (a first look at whether
% distinct waveform populations exist, previewing sorting).

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'filtered') || ~isfield(D, 'spikes')
        error('step4_waveforms:missing', 'Run steps 1-3 first.');
    end

    fs   = D.fs;
    N    = size(D.filtered, 1);
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;
    pol  = lower(string(P.detectPolarity));

    npre   = round(P.wfPreMs  * 1e-3 * fs);
    npost  = round(P.wfPostMs * 1e-3 * fs);
    nalign = round(P.wfAlignSearchMs * 1e-3 * fs);
    refr   = max(1, round(P.refractoryMs * 1e-3 * fs));
    wf_t_ms = (-npre:npost) / fs * 1e3;

    fprintf('[step4] Waveforms: window [-%g, +%g] ms, align +/- %g ms.\n', ...
        P.wfPreMs, P.wfPostMs, P.wfAlignSearchMs);

    for k = 1:nCh
        xf = D.filtered(:, k);
        c0 = D.spikes(k).centers(:);

        % --- re-align each detection to the local extremum ---
        aligned = c0;
        for i = 1:numel(c0)
            lo = max(1, c0(i) - nalign);
            hi = min(N, c0(i) + nalign);
            seg = xf(lo:hi);
            switch pol
                case "pos";  [~, m] = max(seg);
                case "both"; [~, m] = max(abs(seg));
                otherwise;   [~, m] = min(seg);   % 'neg'
            end
            aligned(i) = lo + m - 1;
        end

        % --- keep in-bounds, finite-window events ---
        ok = aligned > npre & aligned <= N - npost;
        aligned = aligned(ok);

        % --- re-enforce refractory on aligned peaks (keep strongest) ---
        if ~isempty(aligned)
            strength = abs(xf(aligned));
            aligned = enforce_refractory(aligned, strength, refr);
        end

        % --- extract waveforms + features ---
        nS = numel(aligned);
        W  = nan(nS, npre + npost + 1);
        Vpp = nan(nS, 1);
        wid = nan(nS, 1);
        finiteW = true(nS, 1);
        for i = 1:nS
            w = xf(aligned(i) - npre : aligned(i) + npost);
            if any(~isfinite(w)); finiteW(i) = false; continue; end
            W(i, :)  = w(:).';
            Vpp(i)   = (max(w) - min(w)) * 1e6;
            wid(i)   = fwhm_ms(w, fs);
        end

        % --- screening: keep spike-like events, drop noise / artifact ---
        ampLow  = Vpp < P.minAmpUV;
        ampHigh = Vpp > P.maxAmpUV;
        widLow  = wid < P.minWidthMs;
        widHigh = wid > P.maxWidthMs;
        accept  = finiteW & ~ampLow & ~ampHigh & ~widLow & ~widHigh;

        nDet = nS;
        W = W(accept, :); Vpp = Vpp(accept); wid = wid(accept); aligned = aligned(accept);

        if isempty(W)
            meanW = []; stdW = [];
        else
            meanW = mean(W, 1, 'omitnan');
            stdW  = std(W, 0, 1, 'omitnan');
        end

        validSec = D.spikes(k).validSec;
        D.spikes(k).wf_t_ms        = wf_t_ms;
        D.spikes(k).waveforms      = W;
        D.spikes(k).alignedCenters = aligned;
        D.spikes(k).alignedTimes   = (aligned - 1) / fs;
        D.spikes(k).Vpp_uv         = Vpp;
        D.spikes(k).width_ms       = wid;
        D.spikes(k).meanWaveform   = meanW;
        D.spikes(k).stdWaveform    = stdW;
        D.spikes(k).nDetected      = nDet;
        D.spikes(k).nSpikes        = size(W, 1);
        D.spikes(k).rate_hz        = size(W, 1) / max(validSec, eps);
        D.spikes(k).screen = struct( ...
            'minAmpUV', P.minAmpUV, 'maxAmpUV', P.maxAmpUV, ...
            'minWidthMs', P.minWidthMs, 'maxWidthMs', P.maxWidthMs, ...
            'nAmpHigh', nnz(finiteW & ampHigh), 'nAmpLow', nnz(finiteW & ampLow), ...
            'nWidHigh', nnz(finiteW & widHigh), 'nWidLow', nnz(finiteW & widLow));

        fprintf(['        ch %d (%s): %d detected -> %d accepted (%.2f spk/s) | ' ...
                 'rej Vpp>%g:%d  FWHM>%g:%d  FWHM<%g:%d | median Vpp %.1f uV, FWHM %.2f ms\n'], ...
            ch(k), lab{k}, nDet, size(W,1), D.spikes(k).rate_hz, ...
            P.maxAmpUV, nnz(finiteW&ampHigh), P.maxWidthMs, nnz(finiteW&widHigh), ...
            P.minWidthMs, nnz(finiteW&widLow), median(Vpp,'omitnan'), median(wid,'omitnan'));
    end

    if plotMode
        for k = 1:nCh
            plot_waveforms(D, k);
        end
    end
end

% ========================================================================
function w_ms = fwhm_ms(w, fs)
% Full width at half-amplitude of |w| (about the aligned peak).
    wa = abs(w(:));
    pk = max(wa);
    if pk <= 0; w_ms = NaN; return; end
    above = find(wa >= 0.5 * pk);
    if isempty(above); w_ms = NaN; return; end
    w_ms = (above(end) - above(1) + 1) / fs * 1e3;
end

% ========================================================================
function centersOut = enforce_refractory(centers, strength, refr)
% Greedy: within any run spaced <= refr, keep the strongest only.
    [cs, ord] = sort(centers);
    ss = strength(ord);
    keep = true(size(cs));
    i = 1;
    while i <= numel(cs)
        j = i + 1;
        grp = i;
        while j <= numel(cs) && (cs(j) - cs(j-1) <= refr)
            grp(end+1) = j; %#ok<AGROW>
            j = j + 1;
        end
        [~, im] = max(ss(grp));
        loser = setdiff(grp, grp(im));
        keep(loser) = false;
        i = j;
    end
    centersOut = sort(cs(keep));
end

% ========================================================================
function plot_waveforms(D, k)
    sp  = D.spikes(k);
    lab = sp.label;
    W   = sp.waveforms;
    tms = sp.wf_t_ms;

    fig = figure('Color', 'w', 'Name', sprintf('Step 4 — %s', lab), ...
        'Position', [180 120 1150 620]);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 4 waveforms  |  %s  |  %d spikes', lab, size(W,1)), ...
        'Interpreter', 'none');

    % ---- waveform overlay + mean +/- std ----
    nexttile;
    if isempty(W)
        text(0.4, 0.5, 'No waveforms', 'Units', 'normalized'); axis off;
    else
        nShow = min(200, size(W,1));
        sel = round(linspace(1, size(W,1), nShow));
        plot(tms, W(sel, :).' * 1e6, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.3); hold on;
        plot(tms, sp.meanWaveform * 1e6, 'k', 'LineWidth', 2);
        plot(tms, (sp.meanWaveform + sp.stdWaveform) * 1e6, 'r--', 'LineWidth', 1);
        plot(tms, (sp.meanWaveform - sp.stdWaveform) * 1e6, 'r--', 'LineWidth', 1);
        grid on; xlabel('Time from peak (ms)'); ylabel('Amplitude (\muV)');
        title('Aligned waveforms (mean \pm std)');
        xlim([tms(1) tms(end)]);
    end

    % ---- amplitude vs width scatter ----
    nexttile;
    if isempty(W)
        axis off;
    else
        plot(sp.width_ms, sp.Vpp_uv, '.', 'Color', [0.2 0.4 0.8], 'MarkerSize', 6);
        hold on;
        if isfield(sp, 'screen')
            sc = sp.screen;
            rectangle('Position', [sc.minWidthMs, sc.minAmpUV, ...
                sc.maxWidthMs - sc.minWidthMs, sc.maxAmpUV - sc.minAmpUV], ...
                'EdgeColor', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 1);
        end
        grid on; xlabel('FWHM (ms)'); ylabel('V_{pp} (\muV)');
        title('Amplitude vs width (accepted; red box = screening bounds)');
    end
end
