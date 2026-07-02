function plotWindowedViolins()
%PLOTWINDOWEDVIOLINS  Box + KDE-violin plots of windowed averages,
% normalised to baseline. Same data + cache + UI as plotWindowedMetrics,
% but the layout is swapped: WINDOWS are on the x-axis and conditions
% appear as coloured sub-violins inside each window — so you can read
% how the spread of normalised values evolves from the first 1 min of
% recovery through the full duration.
%
% One figure per (metric, user-group). For each figure:
%   x-axis     : windows (1 min, 2 min, 5 min, 10 min, full, ...)
%   subgroup   : conditions in that user-group, side-by-side per window
%   each cell  : kernel-density violin + inner box of the normalised
%                values across animals
%
% Reuses gemsplots_queue.mat for files/groups/metric specs, and
% gemsplots_windowed_cache.mat for the per-(file, spec, window) scalars
% — so if plotWindowedMetrics has already populated the cache, this is
% basically instant.

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

    %% ---- validate groups against queue conditions ----
    queueConds = unique(state.condition);
    allGroupConds = unique(horzcat(groups{:}));
    unknown = setdiff(allGroupConds, queueConds);
    if ~isempty(unknown)
        fprintf('[plotWindowedViolins] !! UNKNOWN condition(s) in groups (not in queue, will plot blank):\n');
        for k = 1:numel(unknown)
            fprintf('     %s\n', unknown{k});
        end
        fprintf('   Detected queue conditions: %s\n', strjoin(queueConds, ', '));
        fprintf('   Likely typo — edit groups via the UI and re-run.\n');
    end

    %% ---- metrics (same as plotWindowedMetrics) ----
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
    hasSeries = false(numel(metricSpecs),1);
    for k = 1:numel(metricSpecs)
        hasSeries(k) = ~isempty(strtrim(char(metricSpecs(k).seriesField))) && ...
                       ~isempty(strtrim(char(metricSpecs(k).timeField)));
    end
    if ~any(hasSeries)
        warndlg('No metric has both Series field and Time field.', ...
            'plotWindowedViolins');
        return;
    end
    metricSpecs = metricSpecs(hasSeries);

    %% ---- windows (same as plotWindowedMetrics) ----
    if ~isfield(state,'windowsStr') || isempty(state.windowsStr)
        state.windowsStr = '60, 120, 300, 600, full';
    end
    [winSecs, winLabels, ok, winStr] = defineWindowsUI(state.windowsStr);
    if ~ok, fprintf('Cancelled.\n'); return; end
    state.windowsStr = winStr;
    save(cacheFile, '-struct', 'state');

    %% ---- output prefs ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save violin figures (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir,0), saveDir = '';
    else
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
    end

    defAutoClose = '5'; if isempty(saveDir), defAutoClose = '0'; end
    ac = inputdlg({sprintf(['Auto-close each figure after N seconds (0 = keep open).\n' ...
        '%d figures will be produced.'], numel(metricSpecs) * numel(groups))}, ...
        'Auto-close', [1 60], {defAutoClose});
    autoCloseSec = 0;
    if ~isempty(ac)
        v = str2double(ac{1});
        if ~isnan(v) && v >= 0, autoCloseSec = v; end
    end

    %% ---- compute (shared series cache with plotWindowedMetrics + plotTimeTraces) ----
    seriesCacheFile = fullfile(here, 'gemsplots_series_cache.mat');
    [baselineByFile, windowByFile, hitCount, loadCount] = ...
        loadAllWindowedCached(state, metricSpecs, winSecs, seriesCacheFile);
    fprintf('[plotWindowedViolins] Series values ready: %d cache hit(s), %d fresh load(s).\n', ...
        hitCount, loadCount);

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

    wbR = waitbar(0, 'Rendering violin plots...', ...
        'Name','plotWindowedViolins: rendering');
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
                figH = renderViolinFigure(state, animalsAll, ...
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

    fprintf('Rendered %d violin figure(s).\n', step);
end


% ==========================================================================
% =========================== render ======================================
% ==========================================================================
function figH = renderViolinFigure(state, animalsAll, conds, spec, ...
        groupIdx, baselineThisMetric, windowThisMetric, winSecs, winLabels)
% One figure for (metric, user-group). X = windows. Subgroups = conditions.
% Each cell holds normalised values (rec_window_mean - baseline_mean)/baseline_mean
% across animals.

    figH = gobjects(0);
    nC = numel(conds);
    if nC == 0, return; end
    nW = numel(winSecs);
    nAn = numel(animalsAll);

    % data{ci,wi}(p) = normalised window value for the p-th recovery TRIAL of the
    % condition (NOT per animal): pair each recovery to ITS OWN baseline by trial
    % stem (fallback same animal+condition), so two trials of the same condition
    % from the same animal stay as SEPARATE observations.
    stems = cellfun(@trialStem, state.files, 'UniformOutput', false);

    data = cell(nC, nW);
    reasons = cell(0,1);
    condUsable = false(nC, 1);
    for ci = 1:nC
        cond   = conds{ci};
        recIdx = find( strcmpi(state.condition, cond) & strcmpi(state.phase, 'recovery') );
        nT     = numel(recIdx);
        for wi = 1:nW, data{ci,wi} = nan(nT,1); end
        for p = 1:nT
            iRec    = recIdx(p);
            aniName = state.animal{iRec};

            % this trial's own baseline: same trial stem, else same animal+condition
            iBase = find( strcmpi(state.condition, cond) & strcmpi(state.phase, 'baseline') ...
                        & strcmp(stems, stems{iRec}), 1 );
            if isempty(iBase)
                iBase = find( strcmpi(state.condition, cond) & strcmpi(state.phase, 'baseline') ...
                            & strcmpi(state.animal, aniName), 1 );
            end
            if isempty(iBase)
                reasons{end+1,1} = sprintf( ...
                    '%s / animal=%s : no BASELINE paired to recovery=%s', ...
                    cond, aniName, briefName(state.files{iRec})); %#ok<AGROW>
                continue;
            end

            bm = baselineThisMetric(iBase);
            if isnan(bm) || bm == 0
                reasons{end+1,1} = sprintf( ...
                    '%s / animal=%s : baseline mean is %s for %s', ...
                    cond, aniName, ternary(isnan(bm),'NaN','zero'), ...
                    briefName(state.files{iBase})); %#ok<AGROW>
                continue;
            end

            anyWindow = false;
            for wi = 1:nW
                wm = windowThisMetric(iRec, 1, wi);
                if isnan(wm), continue; end
                data{ci, wi}(p) = (wm - bm) / bm;
                anyWindow = true;
            end
            if anyWindow
                condUsable(ci) = true;
            else
                reasons{end+1,1} = sprintf( ...
                    '%s / animal=%s : every recovery-window mean is NaN for %s', ...
                    cond, aniName, briefName(state.files{iRec})); %#ok<AGROW>
            end
        end
    end

    blankCols = find(~condUsable);
    if ~isempty(blankCols)
        fprintf('[%s — group %d] !! BLANK column(s) in plot:\n', ...
            displayLabel(spec), groupIdx);
        for n = 1:numel(blankCols)
            fprintf('    condition "%s" has no usable (baseline, recovery) pair\n', ...
                conds{blankCols(n)});
        end
    end

    if ~isempty(reasons)
        fprintf('[%s — group %d] %d animal-condition(s) excluded:\n', ...
            displayLabel(spec), groupIdx, numel(reasons));
        for n = 1:min(20, numel(reasons))
            fprintf('    %s\n', reasons{n});
        end
        if numel(reasons) > 20
            fprintf('    ... (+%d more)\n', numel(reasons) - 20);
        end
    end

    dispLab = displayLabel(spec);
    yLab    = sprintf('%s — normalised (rec − base) / base', dispLab);
    figName = sprintf('windowed_violin - group%d - %s', groupIdx, dispLab);
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position',[100 100 1200 720]);

    if all(cellfun(@(c) all(isnan(c(:))), data(:)))
        text(0.5, 0.5, sprintf('no data for group %d / %s', groupIdx, dispLab), ...
            'HorizontalAlignment','center','Units','normalized', ...
            'Interpreter','none');
        axis off;
        return;
    end

    boxViolinPlot(data, ...
        'GroupLabels',    conds, ...         % x-axis = conditions
        'SubgroupLabels', winLabels, ...     % colored sub-violins = windows
        'YLabel',         yLab, ...
        'XLabel',         sprintf('Group %d', groupIdx), ...
        'Title',          sprintf('%s — group %d (windowed spread, baseline-normalised)', ...
                                  dispLab, groupIdx), ...
        'ShowScatter',    true, ...
        'MarkerSize',     40, ...
        'BoxAlpha',       0.20, ...
        'ViolinAlpha',    0.30);
end


% ==========================================================================
% =========================== windows UI ==================================
% ==========================================================================
function [winSecs, winLabels, ok, winStr] = defineWindowsUI(defaultStr)
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
    keep      = false(numel(items),1);
    for k = 1:numel(items)
        if strcmpi(items{k}, 'full')
            winSecs(k)   = Inf;
            winLabels{k} = 'full';
            keep(k) = true;
        else
            v = str2double(items{k});
            if isnan(v) || v <= 0, continue; end
            winSecs(k) = v;
            if v < 60,     winLabels{k} = sprintf('%g s',   v);
            elseif v < 3600, winLabels{k} = sprintf('%g min', v/60);
            else,            winLabels{k} = sprintf('%g hr',  v/3600);
            end
            keep(k) = true;
        end
    end
    winSecs = winSecs(keep);  winLabels = winLabels(keep);
    ok = ~isempty(winSecs);
    if ~ok, warndlg('No valid windows parsed.','Bad input'); end
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

function s = briefName(fp)
    [~, b, e] = fileparts(char(fp));
    s = [b e];
end

function s = trialStem(fp)
% Trial id shared by a baseline/recovery pair = the filename with the phase
% token (and everything after it) stripped, so repeats of the same condition
% from the same animal (e.g. ..._CME1_... vs ..._CME2_...) get DISTINCT stems.
    [~, b] = fileparts(char(fp));
    s = regexprep(b, '_(stim_rec|stim_bl|recovery|baseline|rec|bl)_.*$', '', 'ignorecase');
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
        warning('plotWindowedViolins:saveFig','.fig save failed: %s', ME.message);
    end

    % white background for raster/vector exports
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
        warning('plotWindowedViolins:savePng','.png save failed: %s', ME.message);
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
            warning('plotWindowedViolins:saveSvg','.svg save failed: %s', ME.message);
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
