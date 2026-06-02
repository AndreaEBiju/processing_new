function D = step5_sort_check(D, P, plotMode)
% STEP5_SORT_CHECK  Rigorous test for sortable clusters + detection sanity.
%
%   D = step5_sort_check(D, P, plotMode)
%
% Answers two questions the amplitude/width scatter cannot:
%   (1) Are there separable waveform clusters (sortable units)?
%       Aligned waveforms -> PCA (P.numPCs) -> t-SNE 2-D embedding ->
%       DBSCAN (density clustering, discovers the cluster count, no preset).
%       The ARBITER is the per-cluster mean waveform: clusters with visibly
%       different shapes are real units; identical shapes across "clusters"
%       mean t-SNE manufactured islands from noise (a known t-SNE failure).
%   (2) Is detection catching real spikes? An ISI distribution and a gallery
%       of individual waveforms let you judge directly.
%
% All native MATLAB (pca / tsne / dbscan / knnsearch). The feature matrix is
% small (nSpikes x numPCs), so there is no large-data problem here.
%
% Requires D.spikes(k).waveforms (step 4).
%
% Adds D.sortcheck(k): score (PCA), Y (t-SNE), labels (DBSCAN, 0 = noise),
% nClusters, eps.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'spikes') || ~isfield(D.spikes, 'waveforms')
        error('step5_sort_check:missing', 'Run step4_waveforms first.');
    end

    ch  = D.neuralChannels;
    nCh = numel(ch);
    lab = D.channelLabels;

    D.sortcheck = repmat(struct('score', [], 'Y', [], 'labels', [], ...
        'nClusters', 0, 'eps', NaN), 1, nCh);

    fprintf('[step5] Sort check: PCA(%d) -> t-SNE -> DBSCAN(minPts=%d).\n', ...
        P.numPCs, P.dbscanMinPts);

    for k = 1:nCh
        W = D.spikes(k).waveforms;
        nS = size(W, 1);
        if nS < max(50, 3*P.dbscanMinPts)
            fprintf('        ch %d (%s): only %d spikes; skipping cluster test.\n', ...
                ch(k), lab{k}, nS);
            continue;
        end

        % --- PCA features ---
        W0 = W - mean(W, 1, 'omitnan');
        W0(~isfinite(W0)) = 0;
        [~, score] = pca(W0);
        nPC = min(P.numPCs, size(score, 2));
        X = score(:, 1:nPC);

        % --- t-SNE embedding (reproducible) ---
        rng(1);
        perp = min(P.tsnePerplexity, max(5, floor(nS/4)));
        Y = tsne(X, 'Perplexity', perp, 'Standardize', true);

        % --- DBSCAN on the embedding ---
        if isempty(P.dbscanEps)
            [~, kd] = knnsearch(Y, Y, 'K', P.dbscanMinPts + 1);
            epsUse = prctile(kd(:, end), P.dbscanEpsPct);
        else
            epsUse = P.dbscanEps;
        end
        labels = dbscan(Y, epsUse, P.dbscanMinPts);   % -1 = noise
        labels(labels < 0) = 0;                        % 0 = noise/unassigned
        nClust = numel(unique(labels(labels > 0)));

        D.sortcheck(k).score     = X;
        D.sortcheck(k).Y         = Y;
        D.sortcheck(k).labels    = labels;
        D.sortcheck(k).nClusters = nClust;
        D.sortcheck(k).eps       = epsUse;

        fprintf('        ch %d (%s): %d spikes -> %d cluster(s) (%.0f%% noise), eps=%.2f\n', ...
            ch(k), lab{k}, nS, nClust, 100*mean(labels==0), epsUse);
    end

    if plotMode
        for k = 1:nCh
            if ~isempty(D.sortcheck(k).Y)
                plot_sortcheck(D, k);
            end
        end
    end
end

% ========================================================================
function plot_sortcheck(D, k)
    sp  = D.spikes(k);
    sc  = D.sortcheck(k);
    lab = sp.label;
    W   = sp.waveforms;
    tms = sp.wf_t_ms;
    labels = sc.labels;
    ids = unique(labels(labels > 0));
    cmap = lines(max(1, numel(ids)));

    fig = figure('Color', 'w', 'Name', sprintf('Step 5 — %s', lab), ...
        'Position', [120 80 1250 860]);
    tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 5 sort check  |  %s  |  %d spikes, %d cluster(s)', ...
        lab, size(W,1), sc.nClusters), 'Interpreter', 'none');

    % ---- t-SNE embedding colored by cluster ----
    nexttile;
    hold on;
    isn = labels == 0;
    if any(isn)
        plot(sc.Y(isn,1), sc.Y(isn,2), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 4);
    end
    for c = 1:numel(ids)
        m = labels == ids(c);
        plot(sc.Y(m,1), sc.Y(m,2), '.', 'Color', cmap(c,:), 'MarkerSize', 5);
    end
    grid on; xlabel('t-SNE 1'); ylabel('t-SNE 2');
    title('t-SNE embedding (grey = noise/unassigned)');

    % ---- per-cluster mean waveforms (the arbiter) ----
    nexttile;
    hold on;
    if isempty(ids)
        plot(tms, mean(W,1,'omitnan')*1e6, 'k', 'LineWidth', 2);
        title('Mean waveform (no clusters found)');
    else
        for c = 1:numel(ids)
            mw = mean(W(labels==ids(c), :), 1, 'omitnan') * 1e6;
            plot(tms, mw, 'Color', cmap(c,:), 'LineWidth', 2, ...
                'DisplayName', sprintf('cl %d (n=%d)', ids(c), nnz(labels==ids(c))));
        end
        legend('show', 'Location', 'best');
        title('Per-cluster mean waveforms (different shapes = real units)');
    end
    grid on; xlabel('Time from peak (ms)'); ylabel('Amplitude (\muV)');
    xlim([tms(1) tms(end)]);

    % ---- ISI distribution ----
    nexttile;
    isi = diff(sort(sp.alignedTimes)) * 1000;   % ms
    isi = isi(isi <= 200);
    if ~isempty(isi)
        histogram(isi, 0:1:50);
        xlabel('ISI (ms)'); ylabel('Count');
        title('Inter-spike-interval distribution');
        grid on;
    else
        axis off; text(0.3,0.5,'No ISI data','Units','normalized');
    end

    % ---- individual waveform gallery (colored by cluster) ----
    nexttile;
    hold on;
    nShow = min(60, size(W,1));
    sel = round(linspace(1, size(W,1), nShow));
    for i = sel
        if labels(i) > 0
            col = cmap(find(ids==labels(i),1), :);
        else
            col = [0.8 0.8 0.8];
        end
        plot(tms, W(i,:)*1e6, 'Color', col, 'LineWidth', 0.4);
    end
    grid on; xlabel('Time from peak (ms)'); ylabel('Amplitude (\muV)');
    title(sprintf('%d individual waveforms (colored by cluster)', nShow));
    xlim([tms(1) tms(end)]);
end
