function plotWindowedMetrics()
%PLOTWINDOWEDMETRICS  Box plots of windowed averages over the recovery
% period, each window normalised to its animal's baseline mean.
%
% Workflow
%   1. buildFileQueue     — same persistent queue as plotAveragedMetrics
%   2. defineGroupsUI     — same condition-grouping UI
%   3. defineMetricsUI    — same metric table, but only metrics that
%                            have both 'Series field' and 'Time field'
%                            filled in are used (the rest are skipped).
%                            Defaults provided cover HR/BR/HRV/pNN5/RMSSD/
%                            SampEn/SD1/SD2/SW rate (1/2/3).
%   4. defineWindowsUI    — comma-separated list of windows in seconds,
%                            optionally with 'full' for the whole
%                            recovery duration. Cached between runs.
%   5. For each (animal, condition):
%        baseline file -> compute baseline_mean over the whole baseline
%                          time series
%        recovery file -> for each window W, compute mean of the recovery
%                          series over t in [t0, t0 + W]   (t0 is the
%                          first sample's time so partial-recovery files
%                          still align). 'full' uses the whole series.
%        norm = (window_mean - baseline_mean) / baseline_mean
%   6. One figure per (metric, user-group). Conditions on x; windows as
%      subgroups; dotted lines pair the same animal across windows;
%      ColorBySubject shades the dots; saves to a user-picked folder if
%      one is given.
%
% Adding a new metric: edit it in the metrics table (defineMetricsUI),
% fill in Label / File suffix / Field / Series field / Time field, save,
% rerun — no code change needed.

    %% ---- queue ----
    here      = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(here, 'gemsplots_queue.mat');

    state = buildFileQueue(cacheFile);
    if isempty(state.files)
        fprintf('No files. Exiting.\n');
        return;
    end

    %% ---- groups ----
    if ~isfield(state,'groups'), state.groups = {}; end
    conditions = uniqueStable(state.condition);
    [groups, ok] = defineGroupsUI(conditions, state.groups);
    if ~ok, fprintf('Cancelled.\n'); return; end
    state.groups = groups;
    save(cacheFile, '-struct', 'state');

    %% ---- metrics ----
    metricSpecs = defineMetricsUI(cacheFile, defaultWindowedMetricSpecs());
    if isempty(metricSpecs)
        fprintf('No metrics defined. Exiting.\n');
        return;
    end
    % keep only those with series + time fields filled
    hasSeries = false(numel(metricSpecs),1);
    for k = 1:numel(metricSpecs)
        hasSeries(k) = ~isempty(strtrim(char(metricSpecs(k).seriesField))) && ...
                       ~isempty(strtrim(char(metricSpecs(k).timeField)));
    end
    if ~any(hasSeries)
        warndlg(['None of the metric specs have both "Series field" and ' ...
                 '"Time field" populated. Open the metrics UI and fill ' ...
                 'them in to enable windowed plotting.'], ...
                 'plotWindowedMetrics');
        return;
    end
    skipped = metricSpecs(~hasSeries);
    metricSpecs = metricSpecs(hasSeries);
    if ~isempty(skipped)
        fprintf('[plotWindowedMetrics] Skipping %d metric(s) without ' ...
                'series/time fields: %s\n', numel(skipped), ...
                strjoin({skipped.label}, ', '));
    end

    %% ---- windows ----
    if ~isfield(state,'windowsStr') || isempty(state.windowsStr)
        state.windowsStr = '60, 120, 300, 600, full';
    end
    [winSecs, winLabels, ok, winStr] = defineWindowsUI(state.windowsStr);
    if ~ok, fprintf('Cancelled.\n'); return; end
    state.windowsStr = winStr;
    save(cacheFile, '-struct', 'state');

    %% ---- output prefs ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save every figure (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir,0)
        saveDir = '';
    else
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
    end

    defAutoClose = '5'; if isempty(saveDir), defAutoClose = '0'; end
    ac = inputdlg({sprintf(['Auto-close each figure after N seconds ' ...
        '(0 = keep all open).\n%d figures will be produced.'], ...
        numel(metricSpecs) * numel(groups))}, ...
        'Auto-close', [1 60], {defAutoClose});
    autoCloseSec = 0;
    if ~isempty(ac)
        v = str2double(ac{1});
        if ~isnan(v) && v >= 0, autoCloseSec = v; end
    end

    %% ---- compute baseline means + windowed recovery means ----
    [baselineByFile, windowByFile] = computeWindowedValues( ...
        state, metricSpecs, winSecs);

    %% ---- render ----
    try
        pubfig_setup('Theme','light', ...
            'BaseFontSize', 24, ...
            'LineWidth',    2.4, ...
            'MarkerSize',   14);
    catch ME
        warning('pubfig_setup failed: %s', ME.message);
    end

    animalsAll = uniqueStable(state.animal);
    totalSteps = numel(metricSpecs) * numel(groups);

    wbR = waitbar(0, 'Rendering windowed plots...', ...
        'Name','plotWindowedMetrics: rendering');
    try, set(wbR,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wbR,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>

    step = 0;
    try
        for k = 1:numel(metricSpecs)
            spec = metricSpecs(k);
            for gi = 1:numel(groups)
                step = step + 1;
                if ishandle(wbR)
                    msg = sprintf('[%d / %d]  %s — group %d', ...
                        step, totalSteps, displayLabel(spec), gi);
                    if ~isempty(saveDir), msg = [msg ' (saving)']; end %#ok<AGROW>
                    waitbar(step/totalSteps, wbR, msg);
                end
                figH = renderWindowedFigure(state, animalsAll, ...
                    groups{gi}, spec, gi, ...
                    baselineByFile(:,k), windowByFile(:,k,:), ...
                    winSecs, winLabels);
                if ~isempty(saveDir) && ~isempty(figH) && ishandle(figH)
                    saveFigureAllFormats(figH, saveDir);
                end
                if autoCloseSec > 0 && ~isempty(figH) && ishandle(figH)
                    drawnow; pause(autoCloseSec);
                    if ishandle(figH), close(figH); end
                end
            end
        end
    catch ME
        if ishandle(wbR), close(wbR); end
        rethrow(ME);
    end
    if ishandle(wbR), close(wbR); end

    fprintf('Rendered %d windowed figure(s).\n', step);
end


% ==========================================================================
% =========================== windows UI ==================================
% ==========================================================================
function [winSecs, winLabels, ok, winStr] = defineWindowsUI(defaultStr)
% Returns parallel arrays of seconds (Inf for 'full') + display labels,
% plus the raw string for caching.
    a = inputdlg({sprintf(['Window sizes for averaging the recovery period.\n' ...
        'Comma-separated, in SECONDS. Use ''full'' for the whole recovery.\n' ...
        'Example: 60, 120, 300, 600, full'])}, ...
        'Define windows', [1 70], {defaultStr});
    if isempty(a), winSecs = []; winLabels = {}; ok = false; winStr = defaultStr; return; end
    winStr = strtrim(a{1});
    items = strsplit(winStr, ',');
    items = strtrim(items);
    items = items(~cellfun('isempty', items));

    winSecs   = nan(numel(items),1);
    winLabels = cell(numel(items),1);
    keep = false(numel(items),1);
    for k = 1:numel(items)
        if strcmpi(items{k}, 'full')
            winSecs(k)   = Inf;
            winLabels{k} = 'full';
            keep(k) = true;
        else
            v = str2double(items{k});
            if isnan(v) || v <= 0, continue; end
            winSecs(k) = v;
            if v < 60
                winLabels{k} = sprintf('%g s', v);
            elseif v < 3600
                winLabels{k} = sprintf('%g min', v/60);
            else
                winLabels{k} = sprintf('%g hr', v/3600);
            end
            keep(k) = true;
        end
    end
    winSecs   = winSecs(keep);
    winLabels = winLabels(keep);
    ok = ~isempty(winSecs);
    if ~ok
        warndlg('No valid windows parsed.','Bad input');
    end
end


% ==========================================================================
% =========================== compute =====================================
% ==========================================================================
function [baselineByFile, windowByFile] = computeWindowedValues(state, metricSpecs, winSecs)
% Returns
%   baselineByFile : nFiles x nMetrics       (mean of the FULL series)
%   windowByFile   : nFiles x nMetrics x nW  (mean of series over each window
%                                              starting from t0 = first sample)
%
% For the windowed slot, 'full' (Inf) uses every sample. For finite
% windows, samples with t <= t0 + W are averaged (omitnan).
%
% This is the bulk-load step. Uses a waitbar to show per-file progress.

    nFiles   = numel(state.files);
    nMetrics = numel(metricSpecs);
    nW       = numel(winSecs);
    baselineByFile = nan(nFiles, nMetrics);
    windowByFile   = nan(nFiles, nMetrics, nW);

    wb = waitbar(0, sprintf('Loading time series from %d file(s)...', nFiles), ...
        'Name','plotWindowedMetrics: loading');
    try, set(wb,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wb,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>
    cleanup = onCleanup(@() safeCloseWb(wb)); %#ok<NASGU>

    for i = 1:nFiles
        if ishandle(wb)
            [~, b, e] = fileparts(state.files{i});
            waitbar(i/nFiles, wb, sprintf('[%d / %d] %s', i, nFiles, [b e]));
        end
        for k = 1:nMetrics
            spec = metricSpecs(k);
            [y, t] = loadMetricSeries(state.files{i}, spec);
            if isempty(y) || isempty(t), continue; end
            % full-series mean -> used as baseline reference if this file
            % happens to be a baseline source
            baselineByFile(i,k) = mean(y, 'omitnan');
            % per-window means -> used if this file happens to be a
            % recovery source. t is in the file's own units (seconds).
            t0 = t(1);
            for w = 1:nW
                W = winSecs(w);
                if isinf(W)
                    sel = true(size(t));
                else
                    sel = (t - t0) <= W;
                end
                if any(sel)
                    windowByFile(i,k,w) = mean(y(sel), 'omitnan');
                end
            end
        end
    end
end


% ==========================================================================
% =========================== render ======================================
% ==========================================================================
function figH = renderWindowedFigure(state, animalsAll, conds, spec, ...
        groupIdx, baselineThisMetric, windowThisMetric, winSecs, winLabels)
% One figure for (metric, user-group). For each (condition, animal),
% find the baseline file -> baseline mean; find the recovery file ->
% windowed means; emit normalised values (rec - base)/base.
%   baselineThisMetric : nFiles x 1   (baselineByFile(:,k))
%   windowThisMetric   : nFiles x 1 x nW

    figH = gobjects(0);
    nC = numel(conds);
    if nC == 0, return; end
    nW = numel(winSecs);
    nAn = numel(animalsAll);

    % data{ci, wi}(ai) = normalised window value for that condition/window/animal
    data = cell(nC, nW);
    for ci = 1:nC
        for wi = 1:nW
            data{ci,wi} = nan(nAn,1);
        end
    end

    missing = cell(0,1);
    for ci = 1:nC
        cond = conds{ci};
        for ai = 1:nAn
            aniName = animalsAll{ai};
            iBase = find( strcmpi(state.condition, cond)     & ...
                          strcmpi(state.phase,     'baseline')& ...
                          strcmpi(state.animal,    aniName), 1);
            iRec  = find( strcmpi(state.condition, cond)     & ...
                          strcmpi(state.phase,     'recovery')& ...
                          strcmpi(state.animal,    aniName), 1);
            if isempty(iBase) || isempty(iRec)
                continue;       % expected absence — not a load failure
            end
            bm = baselineThisMetric(iBase);
            if isnan(bm) || bm == 0
                missing{end+1,1} = sprintf( ...
                    '%s / animal=%s  (baseline mean %s)', ...
                    cond, aniName, ternary(isnan(bm), 'NaN', 'zero')); %#ok<AGROW>
                continue;
            end
            for wi = 1:nW
                wm = windowThisMetric(iRec, 1, wi);
                if isnan(wm), continue; end
                data{ci,wi}(ai) = (wm - bm) / bm;
            end
        end
    end

    if ~isempty(missing)
        fprintf('[%s — group %d] %d animal-condition(s) had unusable baseline:\n', ...
            displayLabel(spec), groupIdx, numel(missing));
        for n = 1:min(20, numel(missing))
            fprintf('    %s\n', missing{n});
        end
        if numel(missing) > 20
            fprintf('    ... (+%d more)\n', numel(missing) - 20);
        end
    end

    dispLab = displayLabel(spec);
    yLab    = sprintf('%s — normalised (rec − base) / base', dispLab);
    figName = sprintf('windowed - group%d - %s', groupIdx, dispLab);
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position',[100 100 1200 720]);

    if all(cellfun(@(c) all(isnan(c(:))), data(:)))
        text(0.5, 0.5, sprintf('no data for group %d / %s', groupIdx, dispLab), ...
            'HorizontalAlignment','center','Units','normalized', ...
            'Interpreter','none');
        axis off;
        return;
    end

    boxScatterPlot(data, ...
        'GroupLabels',    conds, ...
        'SubgroupLabels', winLabels, ...
        'YLabel',         yLab, ...
        'XLabel',         sprintf('Group %d', groupIdx), ...
        'Title',          sprintf('%s — group %d (windowed, baseline-normalised)', ...
                                  dispLab, groupIdx), ...
        'ColorBySubject', true, ...
        'ConnectPaired',  true, ...
        'PairLineStyle',  ':', ...
        'MarkerSize',     110);
end


% ==========================================================================
% =========================== defaults / helpers ==========================
% ==========================================================================
function specs = defaultWindowedMetricSpecs()
% Same as plotAveragedMetrics's defaults — defineMetricsUI handles
% migration if the cache already exists.
    specs = [ ...
        ms('HR',             'bpm', 'bpm', '_HRBR.mat',        'avgHeartRate',       'heartRateSeries',    'metrics_t',         'auto', NaN); ...
        ms('Breathing rate', 'bpm', 'bpm', '_HRBR.mat',        'avgBreathRate',      'breathRateSeries',   'metrics_t',         'auto', NaN); ...
        ms('HRV',            's',   'ms',  '_HRVMeasures.mat', 'hrv',                'hrv_series',         'metrics_t',         'auto', NaN); ...
        ms('pNN5',           '%',   '%',   '_HRVMeasures.mat', 'pnn5',               'pnn5_series',        'metrics_t',         'auto', NaN); ...
        ms('RMSSD',          's',   'ms',  '_HRVMeasures.mat', 'rmssd',              'rmssd_series',       'metrics_t',         'auto', NaN); ...
        ms('Sample entropy', '',    '',    '_HRVMeasures.mat', 'sampEn',             'sampEn_series',      'metrics_t',         'auto', NaN); ...
        ms('SD1',            's',   'ms',  '_HRVMeasures.mat', 'sd1',                'sd1_series',         'metrics_t',         'auto', NaN); ...
        ms('SD2',            's',   'ms',  '_HRVMeasures.mat', 'sd2',                'sd2_series',         'metrics_t',         'auto', NaN); ...
        ms('SW rate ANT1',   'cpm', 'cpm', '_slowWaves.mat',   'slowWaveRateSeries', 'slowWaveRateSeries', 'slowWaveRateTime',  'mean', 1);   ...
        ms('SW rate ANT2',   'cpm', 'cpm', '_slowWaves.mat',   'slowWaveRateSeries', 'slowWaveRateSeries', 'slowWaveRateTime',  'mean', 2);   ...
        ms('SW rate ANT3',   'cpm', 'cpm', '_slowWaves.mat',   'slowWaveRateSeries', 'slowWaveRateSeries', 'slowWaveRateTime',  'mean', 3) ];
end

function s = ms(label, unitsIn, unitsOut, suffix, field, seriesField, timeField, aggregator, channel)
    s = struct('label',label,'unitsIn',unitsIn,'unitsOut',unitsOut, ...
               'suffix',suffix,'field',field, ...
               'seriesField',seriesField,'timeField',timeField, ...
               'aggregator',aggregator,'channel',channel);
end

function s = displayLabel(spec)
    lab = strtrim(char(getf(spec,'label')));
    un  = strtrim(char(getf(spec,'unitsOut')));
    if isempty(un), s = lab; else, s = sprintf('%s (%s)', lab, un); end
end

function out = ternary(c,a,b), if c, out = a; else, out = b; end, end

function v = getf(s, name)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = ''; end
end

function u = uniqueStable(c)
    if isempty(c), u = {}; return; end
    if ~iscell(c), c = cellstr(c); end
    [~, ia] = unique(c, 'stable');
    u = c(sort(ia));
end

function safeCloseWb(wb)
    try, if ishandle(wb), close(wb); end, catch, end
end


% ==========================================================================
% =========================== save helper =================================
% ==========================================================================
function saveFigureAllFormats(fig, outDir)
    name = char(get(fig, 'Name'));
    if isempty(name), name = sprintf('figure_%d', fig.Number); end
    safe = regexprep(name, '[^\w\-.()% ]', '_');
    safe = regexprep(safe, '\s+', '_');
    safe = regexprep(safe, '_+', '_');
    safe = strtrim(safe);
    if isempty(safe), safe = sprintf('figure_%d', fig.Number); end
    base = fullfile(outDir, safe);

    try, savefig(fig, [base '.fig']); catch ME
        warning('plotWindowedMetrics:saveFig','.fig save failed: %s', ME.message);
    end
    try, exportgraphics(fig, [base '.png'], 'Resolution', 200); catch ME
        warning('plotWindowedMetrics:savePng','.png save failed: %s', ME.message);
    end
    try, print(fig, [base '.svg'], '-dsvg', '-vector'); catch ME
        warning('plotWindowedMetrics:saveSvg','.svg save failed: %s', ME.message);
    end
end
