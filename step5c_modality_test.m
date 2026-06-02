function D = step5c_modality_test(D, P, plotMode)
% STEP5C_MODALITY_TEST  Is the waveform-feature distribution multimodal or a continuum?
%
%   D = step5c_modality_test(D, P, plotMode)
%
% Settles the discrete-populations-vs-continuum question the right way: a
% formal test for MULTIMODALITY (density gaps), not a forced cluster count.
% Uses Silverman's smoothed-bootstrap test for unimodality on interpretable
% 1-D features -- FWHM (width), peak-to-peak amplitude, and the top waveform
% principal component:
%   * find the smallest KDE bandwidth h_crit that makes the feature unimodal;
%   * smoothed-bootstrap from the KDE at h_crit;
%   * p = fraction of bootstrap samples still multimodal at h_crit.
% p < alpha => significantly multimodal (a real density gap => discrete
% structure on that feature). p large => unimodal => continuum.
%
% Requires D.spikes(k).waveforms, .width_ms, .Vpp_uv (step 4).
% Adds D.modality(k): feature names, p-values, h_crit, verdict.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'spikes') || ~isfield(D.spikes, 'waveforms')
        error('step5c_modality_test:missing', 'Run step4_waveforms first.');
    end

    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    featNames = {'FWHM (ms)', 'Vpp (uV)', 'waveform PC1'};
    D.modality = repmat(struct('features', {featNames}, 'p', [], 'hcrit', [], ...
        'verdict', ''), 1, nCh);

    fprintf('[step5c] Modality test (Silverman, %d bootstraps): multimodal if p<%.2f.\n', ...
        P.modalityNBoot, P.modalityAlpha);

    for k = 1:nCh
        W = D.spikes(k).waveforms;
        nS = size(W, 1);
        if nS < 50
            fprintf('        ch %d (%s): only %d spikes; skipping.\n', ch(k), lab{k}, nS);
            continue;
        end

        % PC1 of waveforms
        W0 = W - mean(W, 1, 'omitnan'); W0(~isfinite(W0)) = 0;
        [~, score] = pca(W0);
        feats = {D.spikes(k).width_ms(:), D.spikes(k).Vpp_uv(:), score(:,1)};

        pvals = nan(1, numel(feats)); hcrit = nan(1, numel(feats));
        for fi = 1:numel(feats)
            [pvals(fi), hcrit(fi)] = silverman_unimodality(feats{fi}, ...
                P.modalityNBoot, P.modalityGrid, P.modalityMaxN);
        end

        multi = find(pvals < P.modalityAlpha);
        if isempty(multi)
            verdict = 'CONTINUUM (all features unimodal)';
        else
            verdict = sprintf('MULTIMODAL on: %s', strjoin(featNames(multi), ', '));
        end

        D.modality(k).p = pvals; D.modality(k).hcrit = hcrit; D.modality(k).verdict = verdict;

        fprintf('        ch %d (%s): %s\n', ch(k), lab{k}, verdict);
        for fi = 1:numel(feats)
            fprintf('            %-14s p = %.3f\n', featNames{fi}, pvals(fi));
        end

        if plotMode
            plot_modality(feats, featNames, pvals, hcrit, P, lab{k}, ch(k));
        end
    end
end

% ========================================================================
function plot_modality(feats, names, pvals, hcrit, P, lab, ch)
    fig = figure('Color', 'w', 'Name', sprintf('Step 5c — %s', lab), ...
        'Position', [150 130 1250 380]);
    tl = tiledlayout(fig, 1, numel(feats), 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 5c modality test  |  %s (ch %d)', lab, ch), 'Interpreter', 'none');
    for fi = 1:numel(feats)
        nexttile;
        x = feats{fi}(isfinite(feats{fi}));
        histogram(x, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9], ...
            'EdgeColor', 'none'); hold on;
        if isfinite(hcrit(fi))
            g = linspace(min(x), max(x), P.modalityGrid);
            plot(g, kde(x, g, hcrit(fi)), 'k', 'LineWidth', 1.8);
        end
        grid on; xlabel(names{fi}); ylabel('density');
        if pvals(fi) < P.modalityAlpha; tag = 'MULTIMODAL'; else; tag = 'unimodal'; end
        title(sprintf('%s: p=%.3f (%s)', names{fi}, pvals(fi), tag));
    end
end

% ========================================================================
function [p, hcrit] = silverman_unimodality(x, B, G, maxN)
% Silverman (1981) smoothed-bootstrap test of H0: distribution is unimodal.
    x = x(isfinite(x));
    if numel(x) > maxN; x = x(randperm(numel(x), maxN)); end
    n = numel(x);
    if n < 30; p = NaN; hcrit = NaN; return; end
    g = linspace(min(x), max(x), G);
    hcrit = critical_bw(x, g, 1);          % smallest bw giving 1 mode
    s2 = var(x); xbar = mean(x);
    cnt = 0;
    for b = 1:B
        y  = x(randi(n, n, 1));
        z  = randn(n, 1);
        xb = xbar + (y - xbar + hcrit * z) / sqrt(1 + hcrit^2 / s2);  % variance-preserving
        gb = linspace(min(xb), max(xb), G);
        if count_modes(kde(xb, gb, hcrit)) > 1; cnt = cnt + 1; end
    end
    p = cnt / B;                            % fraction multimodal at h_crit
end

function h = critical_bw(x, g, kmodes)
% Smallest bandwidth h such that the KDE has <= kmodes modes (binary search;
% mode count is monotone non-increasing in h for a Gaussian kernel).
    s = std(x); if s <= 0; h = eps; return; end
    hhi = s;
    while count_modes(kde(x, g, hhi)) > kmodes && hhi < s*100; hhi = hhi * 2; end
    hlo = s / 1000;
    for it = 1:60
        hm = 0.5 * (hlo + hhi);
        if count_modes(kde(x, g, hm)) > kmodes; hlo = hm; else; hhi = hm; end
    end
    h = hhi;
end

function f = kde(x, g, h)
% Gaussian KDE of column data x on grid g (row), bandwidth h.
    x = x(:); g = g(:).';
    if h <= 0; h = eps; end
    U = (g - x) / h;                        % n x G
    f = sum(exp(-0.5 * U.^2), 1) / (numel(x) * h * sqrt(2*pi));
end

function m = count_modes(f)
% Number of interior local maxima of density f.
    d = sign(diff(f));
    d(d == 0) = 1;                          % treat flats as rising
    m = sum(d(1:end-1) > 0 & d(2:end) < 0);
    if m == 0; m = 1; end                   % monotone -> single (edge) mode
end
