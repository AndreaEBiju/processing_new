function plotRecoveryDynamics()
%PLOTRECOVERYDYNAMICS  Recovery-dynamics features (section 5 of the
% Dashboard_Guide) plotted as box+scatter per condition, one figure per
% (metric, user-group), with the four features in a 2x2 tiled layout.
%
% Features extracted from each animal's RECOVERY time series y(t)
% (where t starts at the first valid sample):
%   peakExc    = max(y) - min(y)
%   slope30    = OLS slope of y(t) on t in [0, 30 s]
%   final30Mean = mean of y over the last 30 s of valid signal
%   AUCdev     = integral of |y - final30Mean| over the whole recovery
%                (trapezoidal, in metric x seconds)
%
% Each panel: conditions on x (only this user-group's conditions),
% baseline NOT plotted (these features describe the recovery on its
% own), ColorBySubject so each animal carries a consistent shade.
%
% Shares gemsplots_series_cache.mat with plotWindowedMetrics /
% plotWindowedViolins / plotTimeTraces / plotSynergyHeatmaps -- the
% features themselves are derived in-memory from the cached series
% (no separate cache needed; recomputing is microseconds once the
% series are loaded).

    %% ---- queue + groups ----
    here      = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(here, 'gemsplots_queue.mat');

    state = buildFileQueue(cacheFile);
    if isempty(state.files), fprintf('No files. Exiting.\n'); return; end

    if ~isfield(state,'groups'), state.groups = {}; end
    conditions = uniqueStable(state.condition);
    [groups, ok] = defineGroupsUI(conditions, state.groups);
    if ~ok, fprintf('Cancelled.\n'); return; end
    state.groups = groups;
    save(cacheFile, '-struct', 'state');

    %% ---- validate groups ----
    queueConds = unique(state.condition);
    allGroupConds = unique(horzcat(groups{:}));
    unknown = setdiff(allGroupConds, queueConds);
    if ~isempty(unknown)
        fprintf('[plotRecoveryDynamics] !! UNKNOWN condition(s) in groups (not in queue):\n');
        for k = 1:numel(unknown), fprintf('     %s\n', unknown{k}); end
    end

    %% ---- metrics (auto-fill series/time defaults) ----
    metricSpecs = defineMetricsUI(cacheFile, defaultWindowedMetricSpecs());
    if isempty(metricSpecs), fprintf('No metrics. Exiting.\n'); return; end

    defaults = defaultWindowedMetricSpecs();
    for k = 1:numel(metricSpecs)
        if isempty(strtrim(char(metricSpecs(k).seriesField))) || ...
           isempty(strtrim(char(metricSpecs(k).timeField)))
            for d = 1:numel(defaults)
                if strcmpi(strtrim(metricSpecs(k).label), strtrim(defaults(d).label))
                    if isempty(strtrim(char(metricSpecs(k).seriesField)))
                        metricSpecs(k).seriesField = defaults(d).seriesField;
                    end
                    if isempty(strtrim(char(metricSpecs(k).timeField)))
                        metricSpecs(k).timeField = defaults(d).timeField;
                    end
                    break;
                end
            end
        end
    end
    hasSeries = arrayfun(@(s) ~isempty(strtrim(char(s.seriesField))) && ...
                              ~isempty(strtrim(char(s.timeField))), metricSpecs);
    if ~any(hasSeries)
        warndlg('No metric has both Series field and Time field.', ...
            'plotRecoveryDynamics');
        return;
    end
    metricSpecs = metricSpecs(hasSeries);

    %% ---- save folder + auto-close ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save recovery-dynamics figures (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir,0), saveDir = '';
    else
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
    end

    defAutoClose = '5'; if isempty(saveDir), defAutoClose = '0'; end
    nFigsExp = numel(metricSpecs) * numel(groups);
    ac = inputdlg({sprintf(['Auto-close each figure after N seconds (0 = keep open).\n' ...
        '%d figures will be produced.'], nFigsExp)}, ...
        'Auto-close', [1 60], {defAutoClose});
    autoCloseSec = 0;
    if ~isempty(ac)
        v = str2double(ac{1});
        if ~isnan(v) && v >= 0, autoCloseSec = v; end
    end

    %% ---- load series (shared cache) ----
    seriesCacheFile = fullfile(here, 'gemsplots_series_cache.mat');
    [seriesByFile, hitCount, loadCount] = ...
        loadAllSeriesCached(state, metricSpecs, seriesCacheFile);
    fprintf('[plotRecoveryDynamics] Series ready: %d cache hit(s), %d fresh load(s).\n', ...
        hitCount, loadCount);

    %% ---- compute dynamic features per (file, metric) ----
    features = computeDynFeatures(state, metricSpecs, seriesByFile);
    fprintf('[plotRecoveryDynamics] Features computed.\n');

    %% ---- render ----
    try
        pubfig_setup('Theme','light', ...
            'BaseFontSize', 18, ...   % smaller than 24 so the 2x2 tiled layout breathes
            'LineWidth',    2.0, ...
            'MarkerSize',   12);
    catch ME
        warning('pubfig_setup failed: %s', ME.message);
    end

    animalsAll = uniqueStable(state.animal);
    totalSteps = numel(metricSpecs) * numel(groups);

    wbR = waitbar(0, 'Rendering recovery-dynamics figures...', ...
        'Name','plotRecoveryDynamics: rendering');
    try, set(wbR,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wbR,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>

    step = 0; nRendered = 0;
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
                figH = renderDynFigure(state, animalsAll, groups{gi}, ...
                    spec, gi, k, features);
                if isempty(figH) || ~ishandle(figH), continue; end
                nRendered = nRendered + 1;
                if ~isempty(saveDir), saveFigureAllFormats(figH, saveDir); end
                if autoCloseSec > 0
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

    fprintf('Rendered %d recovery-dynamics figure(s).\n', nRendered);
end


% ==========================================================================
% =========================== feature extraction ===========================
% ==========================================================================
function features = computeDynFeatures(state, metricSpecs, seriesByFile)
% For every RECOVERY file × every metric, compute the four dynamic
% features from section 5 of the Dashboard_Guide.
%
% Returns a struct with fields .peakExc, .slope30, .final30Mean, .AUCdev,
% each nFiles x nMetrics. Baseline files and any file with no valid
% series stay NaN throughout.

    nFiles   = numel(state.files);
    nMetrics = numel(metricSpecs);
    features.peakExc     = nan(nFiles, nMetrics);
    features.slope30     = nan(nFiles, nMetrics);
    features.final30Mean = nan(nFiles, nMetrics);
    features.AUCdev      = nan(nFiles, nMetrics);

    isRec = strcmpi(state.phase, 'recovery');
    for i = 1:nFiles
        if ~isRec(i), continue; end
        for k = 1:nMetrics
            s = seriesByFile{i,k};
            if isempty(s) || isempty(s.y) || isempty(s.t), continue; end
            y = double(s.y(:));
            t = double(s.t(:));
            valid = ~isnan(y) & ~isnan(t);
            if sum(valid) < 5, continue; end
            y = y(valid);  t = t(valid);
            t = t - t(1);                                     % anchor to 0

            features.peakExc(i,k) = max(y) - min(y);

            % First 30 s slope (OLS)
            sel = t <= 30;
            if sum(sel) >= 3
                p = polyfit(t(sel), y(sel), 1);
                features.slope30(i,k) = p(1);
            end

            % Final 30 s mean
            tEnd = t(end);
            sel = t >= (tEnd - 30);
            if any(sel)
                features.final30Mean(i,k) = mean(y(sel));
            end

            % AUC of |y - final30Mean|
            if ~isnan(features.final30Mean(i,k))
                dev = abs(y - features.final30Mean(i,k));
                features.AUCdev(i,k) = trapz(t, dev);
            end
        end
    end
end


% ==========================================================================
% =========================== render =======================================
% ==========================================================================
function figH = renderDynFigure(state, animalsAll, conds, spec, groupIdx, kMetric, features)
% One figure for one (metric, user-group), 2x2 tiled with the four
% recovery-dynamics features. Each tile: conditions on x, one box per
% condition (no baseline/recovery split -- these are recovery-only
% summaries). ColorBySubject so the same animal keeps its shade rank
% across the four tiles.

    figH = gobjects(0);
    nC = numel(conds);
    if nC == 0, return; end
    nAn = numel(animalsAll);

    featNames  = {'peakExc',          'slope30',           'final30Mean',          'AUCdev'};
    featLabels = {'Peak excursion',   'First-30 s slope',  'Final-30 s mean',      'AUC of deviation'};

    unitsOut = strtrim(char(getf(spec,'unitsOut')));
    if isempty(unitsOut)
        unitFor = repmat({''}, 1, 4);
    else
        unitFor = { ...
            unitsOut, ...                                  % peakExc -> same units
            sprintf('%s / s', unitsOut), ...               % slope30
            unitsOut, ...                                  % final30Mean
            sprintf('%s\\cdot s', unitsOut) };             % AUCdev (kept simple)
        % use plain text instead of LaTeX:
        unitFor{4} = sprintf('%s*s', unitsOut);
    end

    dispLab = displayLabel(spec);
    figName = sprintf('dyn - group%d - %s', groupIdx, dispLab);
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position', [60 60 1400 900]);

    tl = tiledlayout(figH, 2, 2, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s — group %d  (recovery dynamics features)', dispLab, groupIdx), ...
        'Interpreter','none');

    % Build per-feature, per-condition vectors over all queued animals
    % (one entry per animal; NaN where no recovery file present)
    droppedAll = cell(0,1);
    for f = 1:4
        ax = nexttile(tl);
        fname = featNames{f};
        flabel = featLabels{f};
        funit  = unitFor{f};
        yLab   = ternary(isempty(funit), flabel, sprintf('%s (%s)', flabel, funit));

        % Collect data per condition
        data = cell(nC, 1);
        for ci = 1:nC
            cond = conds{ci};
            v = nan(nAn, 1);
            for ai = 1:nAn
                aniName = animalsAll{ai};
                idx = find( strcmpi(state.condition, cond)  & ...
                            strcmpi(state.phase,     'recovery') & ...
                            strcmpi(state.animal,    aniName), 1);
                if isempty(idx), continue; end
                val = features.(fname)(idx, kMetric);
                v(ai) = val;
                if isnan(val) && f == 1   % flag once per (cond, animal), only on the first panel
                    [~, fn, fe] = fileparts(state.files{idx});
                    droppedAll{end+1,1} = sprintf( ...
                        '%s / animal=%s  recovery has no usable %s metric (%s)', ...
                        cond, aniName, displayLabel(spec), [fn fe]); %#ok<AGROW>
                end
            end
            data{ci} = v;
        end

        if all(cellfun(@(c) all(isnan(c)), data))
            text(ax, 0.5, 0.5, 'no data', ...
                'HorizontalAlignment','center', 'Units','normalized', ...
                'Interpreter','none');
            axis(ax, 'off');
            title(ax, flabel, 'Interpreter','none');
            continue;
        end

        boxScatterPlot(data, ...
            'Parent',         ax, ...
            'GroupLabels',    conds, ...
            'YLabel',         yLab, ...
            'XLabel',         '', ...
            'Title',          flabel, ...
            'ColorBySubject', true, ...
            'MarkerSize',     80, ...
            'Legend',         'off');
    end

    if ~isempty(droppedAll)
        fprintf('[%s — group %d] %d dropped (cond, animal) for recovery dynamics:\n', ...
            dispLab, groupIdx, numel(droppedAll));
        for n = 1:min(20, numel(droppedAll))
            fprintf('    %s\n', droppedAll{n});
        end
        if numel(droppedAll) > 20
            fprintf('    ... (+%d more)\n', numel(droppedAll) - 20);
        end
    end
end


% ==========================================================================
% =========================== defaults / helpers ==========================
% ==========================================================================
function specs = defaultWindowedMetricSpecs()
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

function v = getf(s, name)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = ''; end
end

function u = uniqueStable(c)
    if isempty(c), u = {}; return; end
    if ~iscell(c), c = cellstr(c); end
    [~, ia] = unique(c, 'stable');
    u = c(sort(ia));
end

function out = ternary(c, a, b), if c, out = a; else, out = b; end, end


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
        warning('plotRecoveryDynamics:saveFig','.fig save failed: %s', ME.message);
    end

    origFigColor = get(fig, 'Color');
    origInvert   = get(fig, 'InvertHardcopy');
    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    axList = findall(fig, 'Type', 'axes');
    origAxColors = cell(numel(axList),1);
    for k = 1:numel(axList)
        origAxColors{k} = get(axList(k), 'Color');
        set(axList(k), 'Color', 'w');
    end
    cleanupRestore = onCleanup(@() restoreFigColors( ...
        fig, origFigColor, origInvert, axList, origAxColors)); %#ok<NASGU>

    try
        exportgraphics(fig, [base '.png'], ...
            'Resolution', 200, 'BackgroundColor', 'white');
    catch ME
        warning('plotRecoveryDynamics:savePng','.png save failed: %s', ME.message);
    end

    svgOk = false;
    try
        exportgraphics(fig, [base '.svg'], ...
            'ContentType', 'vector', 'BackgroundColor', 'white');
        svgOk = true;
    catch
    end
    if ~svgOk
        wState = warning('off', 'all');
        try, print(fig, [base '.svg'], '-dsvg', '-vector'); catch ME
            warning(wState);
            warning('plotRecoveryDynamics:saveSvg','.svg save failed: %s', ME.message);
        end
        warning(wState);
    end
end

function restoreFigColors(fig, figColor, invert, axList, axColors)
    try
        if isgraphics(fig)
            set(fig, 'Color', figColor, 'InvertHardcopy', invert);
        end
    catch, end
    for k = 1:numel(axList)
        try
            if isgraphics(axList(k))
                set(axList(k), 'Color', axColors{k});
            end
        catch, end
    end
end
