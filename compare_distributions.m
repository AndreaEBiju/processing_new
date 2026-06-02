function compare_distributions(files)
% COMPARE_DISTRIBUTIONS  Cross-condition waveform-feature distributions.
%
%   compare_distributions()        % GUI: add summary files across folders
%   compare_distributions(files)   % cellstr of *_summary.mat paths
%
% Overlays the per-spike FWHM and peak-to-peak amplitude distributions across
% conditions and runs a two-sample Kolmogorov-Smirnov test per feature. This
% is the well-posed way to ask "does recovery recruit different fibers?" -- a
% shift toward narrow widths / larger amplitudes indicates faster/larger
% fibers recruited, visible as a distributional shift even in a continuum
% (no clustering required).

    if nargin < 1 || isempty(files)
        files = {}; startDir = '';
        while true
            [f, p] = uigetfile('*_summary.mat', ...
                'Select summary file(s) — add more from other folders next', ...
                startDir, 'MultiSelect', 'on');
            if ~isequal(f, 0)
                if ischar(f); f = {f}; end
                files = [files, cellfun(@(x) fullfile(p, x), f, 'UniformOutput', false)]; %#ok<AGROW>
                startDir = p;
            end
            if ~strcmp(questdlg(sprintf('%d file(s) so far. Add more?', numel(files)), ...
                    'Add more?', 'Add more...', 'Done', 'Done'), 'Add more...'); break; end
        end
        files = unique(files, 'stable');
        if isempty(files); fprintf('No files selected.\n'); return; end
    end

    nC = numel(files); S = cell(1, nC); names = cell(1, nC);
    for i = 1:nC
        L = load(files{i}); S{i} = L.summary; names{i} = L.summary.condition;
    end
    cmap = lines(nC);
    labels = S{1}.labels;
    featGet = {@(s,k) s.featWidth{k}, @(s,k) s.featVpp{k}};
    featName = {'FWHM (ms)', 'V_{pp} (\muV)'};

    fprintf('\n===== Cross-condition feature distributions (KS test) =====\n');
    for li = 1:numel(labels)
        lab = labels{li};
        fig = figure('Color', 'w', 'Name', sprintf('Distributions — %s', lab), ...
            'Position', [140 130 1150 460]);
        tl = tiledlayout(fig, 1, numel(featGet), 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl, sprintf('Waveform-feature distributions  |  channel %s', lab), 'Interpreter', 'none');

        for fj = 1:numel(featGet)
            nexttile; hold on;
            data = cell(1, nC);
            for i = 1:nC
                k = find(strcmp(S{i}.labels, lab), 1);
                if isempty(k); continue; end
                x = featGet{fj}(S{i}, k); x = x(isfinite(x));
                data{i} = x;
                histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
                    'EdgeColor', cmap(i,:), 'LineWidth', 1.8, 'DisplayName', names{i});
            end
            grid on; xlabel(featName{fj}); ylabel('density');
            % KS test between the first two conditions (if present)
            ttl = featName{fj};
            if nC >= 2 && ~isempty(data{1}) && ~isempty(data{2})
                [~, pks, ks] = kstest2(data{1}, data{2});
                ttl = sprintf('%s  |  KS p=%.3g (D=%.2f)', featName{fj}, pks, ks);
                fprintf('  ch %s, %s: %s vs %s  KS D=%.3f, p=%.3g\n', ...
                    lab, featName{fj}, names{1}, names{2}, ks, pks);
            end
            title(ttl); legend('show', 'Location', 'northeast');
        end
    end
    fprintf(['\nInterpret: a recovery shift toward smaller FWHM / larger Vpp = recruitment\n' ...
             'of faster/larger fibers. KS p<0.05 = distributions differ.\n']);
end
