function D = step6_population_summary(D, P, plotMode)
% STEP6_POPULATION_SUMMARY  Single-population (no clustering) per-channel summary.
%
%   D = step6_population_summary(D, P, plotMode)
%
% Treats ALL accepted spikes (step 4) as one C-like population -- no
% clustering. Per neural channel it computes and plots:
%   1. Firing-rate time trace, 1 s bins, blanking-corrected (count / valid
%      seconds per bin) so masked/dropout time does not depress the rate.
%   2. Inter-spike-interval (ISI) histogram, with the refractory line marked.
%   3. Mean waveform with a shaded percentile band (P.wfBandPct) showing the
%      spread/range -- individual waveforms are NOT drawn.
%   4. Raster: spike times wrapped into rows of P.rasterRowSec seconds (one
%      row per epoch), a readable raster for a long continuous recording.
%      (No stimulus trials exist here; if you have trial markers, rows can be
%      aligned to them instead.)
%
% Requires D.spikes(k).alignedTimes / .waveforms / .wf_t_ms (step 4) and
% D.validMask (step 2).
%
% Adds D.popsummary(k): fr_t, fr_hz, isi_ms, meanWaveform, bandLo, bandHi,
% wf_t_ms -- so the bulk stage can aggregate these later.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'spikes') || ~isfield(D.spikes, 'alignedTimes')
        error('step6_population_summary:missing', 'Run step4_waveforms first.');
    end

    fs = D.fs; N = size(D.filtered, 1);
    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    cond = ''; if isfield(D, 'condition'); cond = [' | ' D.condition]; end
    binN = max(1, round(P.frBinSec * fs));

    D.popsummary = repmat(struct('fr_t', [], 'fr_hz', [], 'isi_ms', [], ...
        'meanWaveform', [], 'bandLo', [], 'bandHi', [], 'wf_t_ms', []), 1, nCh);

    for k = 1:nCh
        st = D.spikes(k).alignedTimes(:);
        valid = D.validMask(:, k);
        W = D.spikes(k).waveforms;
        tms = D.spikes(k).wf_t_ms;

        % --- 1. firing-rate trace (blanking-corrected) ---
        nBins = ceil(N / binN);
        fr_t = nan(nBins, 1); fr = nan(nBins, 1);
        edges = (0:nBins) * binN;                 % sample edges
        counts = histcounts(round(st * fs) + 1, edges + 0.5);
        for b = 1:nBins
            i0 = edges(b) + 1; i1 = min(N, edges(b+1));
            validSec = nnz(valid(i0:i1)) / fs;
            if validSec > 0; fr(b) = counts(b) / validSec; end
            fr_t(b) = ((i0 + i1) / 2 - 1) / fs;
        end

        % --- 2. ISI ---
        isi = diff(sort(st)) * 1000;              % ms

        % --- 3. mean waveform + percentile band ---
        if isempty(W)
            mw = []; lo = []; hi = [];
        else
            mw = mean(W, 1, 'omitnan');
            pb = prctile(W, P.wfBandPct, 1);       % 2 x L
            lo = pb(1, :); hi = pb(2, :);
        end

        D.popsummary(k) = struct('fr_t', fr_t, 'fr_hz', fr, 'isi_ms', isi, ...
            'meanWaveform', mw, 'bandLo', lo, 'bandHi', hi, 'wf_t_ms', tms);

        fprintf('[step6] ch %d (%s): %d spikes | mean rate %.2f spk/s | median ISI %.1f ms\n', ...
            ch(k), lab{k}, numel(st), mean(fr, 'omitnan'), median(isi, 'omitnan'));

        if plotMode
            plot_population(D, k, st, valid, cond);
        end
    end
end

% ========================================================================
function plot_population(D, k, st, valid, cond)
    ps = D.popsummary(k); lab = D.channelLabels{k}; fs = D.fs;

    fig = figure('Color', 'w', 'Name', sprintf('Step 6 — %s', lab), ...
        'Position', [120 90 1250 820]);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Population summary  |  %s%s  |  %d spikes', lab, cond, numel(st)), ...
        'Interpreter', 'none');

    % 1. firing rate
    ax1 = nexttile;
    plot(ps.fr_t, ps.fr_hz, 'b', 'LineWidth', 1.0); hold on;
    yline(mean(ps.fr_hz, 'omitnan'), 'k--', 'LineWidth', 1);
    shade_runs(ax1, (0:numel(valid)-1)/fs, ~valid);
    grid on; xlabel('Time (s)'); ylabel('Firing rate (spk/s)');
    title(sprintf('Firing rate (%g s bins, blanking-corrected)', median(diff(ps.fr_t), 'omitnan')));
    xlim([0 (numel(valid)-1)/fs]);

    % 2. ISI histogram
    nexttile;
    if ~isempty(ps.isi_ms)
        edges = 0:1:100;
        histogram(ps.isi_ms(ps.isi_ms <= 100), edges, 'FaceColor', [0.5 0.6 0.85]);
        xline(D.spikes(k).threshAtSpike_uv*0 + 1.0, 'r--');  % refractory ~1 ms marker
        grid on; xlabel('ISI (ms)'); ylabel('Count');
        title(sprintf('ISI histogram (median %.1f ms)', median(ps.isi_ms, 'omitnan')));
    else
        axis off; text(0.3, 0.5, 'No ISI', 'Units', 'normalized');
    end

    % 3. mean waveform + band
    nexttile;
    if ~isempty(ps.meanWaveform)
        t = ps.wf_t_ms;
        fill([t fliplr(t)], [ps.bandLo fliplr(ps.bandHi)]*1e6, [0.8 0.85 0.95], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.8); hold on;
        plot(t, ps.meanWaveform*1e6, 'b', 'LineWidth', 2);
        grid on; xlabel('Time from peak (ms)'); ylabel('Amplitude (\muV)');
        title(sprintf('Mean waveform (shaded = %d–%d%% range)', P_band(D)));
        xlim([t(1) t(end)]);
    else
        axis off; text(0.3, 0.5, 'No waveforms', 'Units', 'normalized');
    end

    % 4. raster (wrapped into rows)
    nexttile;
    rowSec = getRowSec(D);
    xr = mod(st, rowSec); yr = floor(st / rowSec);
    if ~isempty(st)
        X = [xr.'; xr.'; nan(1, numel(xr))];
        Y = [yr.'; yr.' + 0.8; nan(1, numel(yr))];
        plot(X(:), Y(:), 'k', 'LineWidth', 0.2);
    end
    set(gca, 'YDir', 'reverse'); grid on;
    xlabel(sprintf('Time within epoch (s)')); ylabel(sprintf('Epoch (× %g s)', rowSec));
    title('Spike raster (recording wrapped into epochs)');
    xlim([0 rowSec]);
end

% ---- tiny helpers to fetch params stored at call (kept simple) ----
function pr = P_band(D) %#ok<INUSD>
    pp = pipeline_params(); pr = pp.wfBandPct;
end
function rs = getRowSec(D) %#ok<INUSD>
    pp = pipeline_params(); rs = pp.rasterRowSec;
end

% ========================================================================
function shade_runs(ax, t, mask)
    if ~any(mask); return; end
    yl = ylim(ax);
    d = diff([0; mask(:); 0]);
    s = find(d == 1); e = find(d == -1) - 1;
    for i = 1:min(numel(s), 300)
        x0 = t(s(i)); x1 = t(min(e(i), numel(t)));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.4 0.4], ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end
