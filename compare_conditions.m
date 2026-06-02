function compare_conditions(files)
% COMPARE_CONDITIONS  Cross-condition comparison of saved summaries.
%
%   compare_conditions()          % GUI: pick 2+ *_summary.mat files
%   compare_conditions(files)     % files = cellstr of summary paths
%
% For each neural channel label present, overlays the noise-corrected excess
% RMS and the accepted-spike rate as time courses across conditions, bars the
% per-condition means, and fits an exponential decay to each excess time
% course (recovery -> baseline). Excess and rate moving together is the
% cross-validation that an elevation is neural rather than motion.

    if nargin < 1 || isempty(files)
        % Gather summaries across MULTIPLE folders: uigetfile multi-select is
        % limited to one folder, so add files in rounds (navigate anywhere
        % each round) until you click Done.
        files = {};
        startDir = '';
        while true
            [f, p] = uigetfile('*_summary.mat', ...
                'Select summary file(s) — add more from other folders next', ...
                startDir, 'MultiSelect', 'on');
            if ~isequal(f, 0)
                if ischar(f); f = {f}; end
                files = [files, cellfun(@(x) fullfile(p, x), f, ...
                    'UniformOutput', false)]; %#ok<AGROW>
                startDir = p;
            end
            choice = questdlg(sprintf('%d file(s) selected so far. Add more from another folder?', ...
                numel(files)), 'Add more conditions?', ...
                'Add more...', 'Done', 'Done');
            if ~strcmp(choice, 'Add more...'); break; end
        end
        files = unique(files, 'stable');
        if isempty(files); fprintf('No files selected.\n'); return; end
        fprintf('Comparing %d summary file(s):\n', numel(files));
        fprintf('  %s\n', files{:});
    end

    nC = numel(files);
    S = cell(1, nC); names = cell(1, nC);
    for i = 1:nC
        L = load(files{i});
        S{i} = L.summary;
        names{i} = L.summary.condition;
    end
    cmap = lines(nC);

    % labels present in the first condition (matched by string across conditions)
    labels = S{1}.labels;

    fprintf('\n===== Cross-condition comparison =====\n');
    for li = 1:numel(labels)
        lab = labels{li};
        fig = figure('Color', 'w', 'Name', sprintf('Compare — %s', lab), ...
            'Position', [120 90 1250 820]);
        tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tl, sprintf('Baseline vs recovery  |  channel %s', lab), 'Interpreter', 'none');

        meanExc = nan(1, nC); meanRate = nan(1, nC);
        eExc = nan(1, nC); lExc = nan(1, nC); eRate = nan(1, nC); lRate = nan(1, nC);

        % ---- excess RMS time courses (raw thin + ~30 s smoothed trend) ----
        nexttile; hold on;
        for i = 1:nC
            k = find(strcmp(S{i}.labels, lab), 1);
            if isempty(k); continue; end
            t = S{i}.t{k}; y = S{i}.excess_uv{k};
            plot(t, y, 'Color', lighten(cmap(i,:)), 'LineWidth', 0.3, 'HandleVisibility', 'off');
            plot(t, movmean(y, smooth_bins(t), 'omitnan'), 'Color', cmap(i,:), ...
                'LineWidth', 2, 'DisplayName', names{i});
            meanExc(i) = S{i}.meanExcess(k);
            [eExc(i), lExc(i)] = early_late(t, y);
        end
        grid on; xlabel('Time (s)'); ylabel('Excess RMS (\muV)');
        title('Noise-corrected excess (bold = ~30 s smoothed trend)');
        legend('show', 'Location', 'northeast');

        % ---- accepted-spike rate time courses ----
        nexttile; hold on;
        for i = 1:nC
            k = find(strcmp(S{i}.labels, lab), 1);
            if isempty(k); continue; end
            t = S{i}.t{k}; y = S{i}.rate_binned{k};
            plot(t, y, 'Color', lighten(cmap(i,:)), 'LineWidth', 0.3, 'HandleVisibility', 'off');
            plot(t, movmean(y, smooth_bins(t), 'omitnan'), 'Color', cmap(i,:), ...
                'LineWidth', 2, 'DisplayName', names{i});
            meanRate(i) = S{i}.meanRate(k);
            [eRate(i), lRate(i)] = early_late(t, y);
        end
        grid on; xlabel('Time (s)'); ylabel('Accepted rate (spk/s)');
        title('Screened spike rate (bold = smoothed; cross-validates excess)');
        legend('show', 'Location', 'northeast');

        % ---- mean excess bar ----
        nexttile;
        bar(categorical(names, names), meanExc, 'FaceColor', [0.4 0.5 0.8]);
        grid on; ylabel('Mean excess RMS (\muV)'); title('Mean excess by condition');

        % ---- mean rate bar ----
        nexttile;
        bar(categorical(names, names), meanRate, 'FaceColor', [0.8 0.5 0.4]);
        grid on; ylabel('Mean rate (spk/s)'); title('Mean rate by condition');

        % ---- console table (early/late = first vs last third) ----
        fprintf('\nChannel %s:\n', lab);
        fprintf('  %-14s  %8s %8s %8s   %8s %8s %8s\n', 'condition', ...
            'mExc', 'earlyE', 'lateE', 'mRate', 'earlyR', 'lateR');
        for i = 1:nC
            fprintf('  %-14s  %8.2f %8.2f %8.2f   %8.2f %8.2f %8.2f\n', ...
                names{i}, meanExc(i), eExc(i), lExc(i), meanRate(i), eRate(i), lRate(i));
        end
    end
    fprintf('\nNote: an excess elevation that the rate trace mirrors is neural;\n');
    fprintf('excess up without rate up suggests residual motion.\n');
end

% ========================================================================
function sb = smooth_bins(t)
% Number of bins spanning ~30 s, for the smoothed trend.
    dt = median(diff(t), 'omitnan');
    if ~isfinite(dt) || dt <= 0; dt = 1; end
    sb = max(1, round(30 / dt));
end

function c = lighten(c0)
    c = 1 - (1 - c0) * 0.35;
end

function [e, l] = early_late(t, y)
% Mean of y over the first third vs the last third of the time span.
    good = isfinite(t) & isfinite(y);
    t = t(good); y = y(good);
    if isempty(t); e = NaN; l = NaN; return; end
    t0 = t(1); t1 = t(end); th = (t1 - t0) / 3;
    e = mean(y(t <= t0 + th), 'omitnan');
    l = mean(y(t >= t1 - th), 'omitnan');
end
