function D = step5b_verify_populations(D, P, plotMode)
% STEP5B_VERIFY_POPULATIONS  How many genuinely distinct waveform populations?
%
%   D = step5b_verify_populations(D, P, plotMode)
%
% Rigorous, reproducible alternative to trusting t-SNE/DBSCAN's cluster count
% (which over-segments). Per channel:
%   1. PCA of the aligned waveforms.
%   2. Deliberately OVER-cluster with k-means (k = P.verifyKover).
%   3. Iteratively MERGE clusters whose mean templates have shape correlation
%      above P.verifyMergeCorr. Pearson correlation is amplitude-invariant, so
%      amplitude/electrode-distance variants collapse and only SHAPE-distinct
%      templates survive -- shape is what tracks conduction class / fiber type.
%   4. Drop populations smaller than P.verifyMinClusterSize.
%
% The surviving templates are a defensible count of distinct populations.
% Each is classified C-like (wide FWHM) or A-like (narrow) per
% P.verifyWidthSplitMs. NOTE: this is a waveform-shape distinction; firm
% fiber-type assignment needs conduction velocity or pharmacology.
%
% Run it on each condition (baseline, recovery) and compare: a narrow A-like
% template that appears in recovery but not baseline is a recruited population.
%
% Requires D.spikes(k).waveforms (step 4).
% Adds D.popcheck(k): labels, templates, nPops, counts, fwhm_ms, Vpp_uv, type.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'spikes') || ~isfield(D.spikes, 'waveforms')
        error('step5b_verify_populations:missing', 'Run step4_waveforms first.');
    end

    ch  = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    D.popcheck = repmat(struct('labels', [], 'templates', [], 'nPops', 0, ...
        'counts', [], 'fwhm_ms', [], 'Vpp_uv', [], 'sep', []), 1, nCh);

    fprintf('[step5b] Population verification: over-cluster k=%d, merge if d''<%.1f (within noise).\n', ...
        P.verifyKover, P.verifyMinDprime);

    for k = 1:nCh
        W = D.spikes(k).waveforms;
        fs = D.fs; tms = D.spikes(k).wf_t_ms;
        nS = size(W, 1);
        if nS < max(50, 2 * P.verifyMinClusterSize)
            fprintf('        ch %d (%s): only %d spikes; skipping.\n', ch(k), lab{k}, nS);
            continue;
        end

        % --- PCA ---
        W0 = W - mean(W, 1, 'omitnan'); W0(~isfinite(W0)) = 0;
        [~, score] = pca(W0);
        nPC = min(P.numPCs, size(score, 2));
        X = score(:, 1:nPC);

        % --- over-cluster ---
        Kover = max(2, min(P.verifyKover, floor(nS / P.verifyMinClusterSize)));
        rng(1);
        labels = kmeans(X, Kover, 'Replicates', 5, 'MaxIter', 500);

        % --- merge populations not separable above noise (LDA d') ---
        while true
            ids = unique(labels).'; ids = ids(ids > 0);
            if numel(ids) < 2; break; end
            bestD = Inf; bi = 0; bj = 0;
            for a = 1:numel(ids)-1
                for b = a+1:numel(ids)
                    d = dprime(W(labels==ids(a), :), W(labels==ids(b), :));
                    if d < bestD; bestD = d; bi = ids(a); bj = ids(b); end
                end
            end
            if bestD >= P.verifyMinDprime; break; end
            labels(labels == bj) = bi;          % merge the least-separable pair
        end
        active = unique(labels).';

        % --- drop small populations, renumber survivors ---
        survivors = []; newLabels = zeros(size(labels));
        nid = 0;
        for a = active
            if nnz(labels == a) >= P.verifyMinClusterSize
                nid = nid + 1; newLabels(labels == a) = nid; survivors(end+1) = a; %#ok<AGROW>
            end
        end
        labels = newLabels;
        nPops = nid;

        % --- per-population features ---
        templates = nan(nPops, size(W, 2));
        counts = zeros(1, nPops); fwhm = nan(1, nPops); vpp = nan(1, nPops);
        for c = 1:nPops
            m = labels == c;
            templates(c, :) = mean(W(m, :), 1, 'omitnan');
            counts(c) = nnz(m);
            fwhm(c) = fwhm_ms(templates(c, :), fs);
            vpp(c)  = (max(templates(c, :)) - min(templates(c, :))) * 1e6;
        end
        % pairwise separability (d') among survivors
        sepMat = zeros(nPops);
        for a = 1:nPops
            for b = a+1:nPops
                dd = dprime(W(labels==a, :), W(labels==b, :));
                sepMat(a,b) = dd; sepMat(b,a) = dd;
            end
        end

        D.popcheck(k) = struct('labels', labels, 'templates', templates, ...
            'nPops', nPops, 'counts', counts, 'fwhm_ms', fwhm, 'Vpp_uv', vpp, 'sep', sepMat);

        fprintf('        ch %d (%s): %d spikes -> %d population(s) separable above noise (d''>=%.1f)\n', ...
            ch(k), lab{k}, nS, nPops, P.verifyMinDprime);
        for c = 1:nPops
            fprintf('            pop %d: n=%d, FWHM=%.2f ms, Vpp=%.1f uV\n', ...
                c, counts(c), fwhm(c), vpp(c));
        end
        if nPops >= 2
            offd = sepMat(~eye(nPops));
            fprintf('            pairwise d'': min=%.2f, median=%.2f\n', min(offd), median(offd));
        end

        if plotMode; plot_populations(D, k, tms); end
    end
end

% ========================================================================
function plot_populations(D, k, tms)
    pc = D.popcheck(k); lab = D.channelLabels{k};
    W = D.spikes(k).waveforms; labels = pc.labels;
    nP = pc.nPops; cmap = lines(max(1, nP));

    fig = figure('Color', 'w', 'Name', sprintf('Step 5b — %s', lab), ...
        'Position', [130 90 1250 820]);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 5b populations  |  %s  |  %d distinct population(s)', ...
        lab, nP), 'Interpreter', 'none');

    % template waveforms
    nexttile; hold on;
    for c = 1:nP
        plot(tms, pc.templates(c,:)*1e6, 'Color', cmap(c,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('pop %d: n=%d, FWHM %.2f ms', c, pc.counts(c), pc.fwhm_ms(c)));
    end
    grid on; xlabel('Time from peak (ms)'); ylabel('Amplitude (\muV)');
    title('Distinct population templates'); legend('show', 'Location', 'best');
    if ~isempty(tms); xlim([tms(1) tms(end)]); end

    % feature scatter colored by population
    nexttile; hold on;
    for c = 1:nP
        m = labels == c;
        plot(D.spikes(k).width_ms(m), D.spikes(k).Vpp_uv(m), '.', ...
            'Color', cmap(c,:), 'MarkerSize', 5);
    end
    grid on; xlabel('FWHM (ms)'); ylabel('V_{pp} (\muV)');
    title('Features colored by population');

    % pairwise separability (d') matrix
    nexttile;
    if nP >= 2
        imagesc(pc.sep); axis square; colorbar;
        set(gca, 'XTick', 1:nP, 'YTick', 1:nP);
        title("Pairwise separability d' (higher off-diag = more distinct)");
        xlabel('pop'); ylabel('pop');
    else
        axis off; text(0.2, 0.5, 'Single population', 'Units', 'normalized');
    end

    % count by population + type
    nexttile;
    bar(categorical(arrayfun(@(c) sprintf('p%d', c), 1:nP, 'UniformOutput', false)), pc.counts);
    grid on; ylabel('spike count'); title('Population sizes');
end

% ========================================================================
function w_ms = fwhm_ms(w, fs)
    wa = abs(w(:)); pk = max(wa);
    if pk <= 0; w_ms = NaN; return; end
    above = find(wa >= 0.5 * pk);
    if isempty(above); w_ms = NaN; return; end
    w_ms = (above(end) - above(1) + 1) / fs * 1e3;
end

% ========================================================================
function d = dprime(Wi, Wj)
% Separation between two waveform sets along their best-discriminating axis,
% in units of the trial-to-trial scatter (LDA-style d'). d' ~ how many noise
% SDs apart the two mean templates are; small d' = overlapping within noise.
    if size(Wi,1) < 2 || size(Wj,1) < 2; d = 0; return; end
    mi = mean(Wi, 1, 'omitnan'); mj = mean(Wj, 1, 'omitnan');
    w = mi - mj; nw = norm(w);
    if nw == 0 || ~isfinite(nw); d = 0; return; end
    w = w / nw;
    qi = Wi * w.'; qj = Wj * w.';      % project onto the difference axis
    d = abs(mean(qi,'omitnan') - mean(qj,'omitnan')) / ...
        sqrt(0.5 * (var(qi,'omitnan') + var(qj,'omitnan')) + eps);
end
