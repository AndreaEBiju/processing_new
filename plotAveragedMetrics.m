function plotAveragedMetrics()
%PLOTAVERAGEDMETRICS  Centralised box+scatter plots of any user-defined
% averaged metrics across stim conditions. One figure per (metric, user-
% defined stim group). Each figure shows the conditions within that group
% on the x-axis, two sub-boxes (baseline, recovery) per condition, with
% ColorBySubject so each animal keeps its shade rank across all figures.
%
% Workflow
%   1. buildFileQueue       file queue UI, persisted in gemsplots_queue.mat
%   2. uiDefineGroups       which conditions belong to which user group
%   3. defineMetricsUI      what to plot — list of metric specs (label,
%                           file suffix, field name, aggregator, channel).
%                           Pre-filled with a sensible default set on
%                           first run; editable + persisted thereafter.
%                           Includes an "Inspect .mat file..." helper so
%                           you can browse variable names in any output
%                           file without leaving the UI.
%   4. loadMetric per row   generic loader (see loadMetric.m)
%   5. render               one figure per (metric, group)
%
% No metric is hardcoded — adding a new one is just a row in the metrics
% table, not a code change.

    %% ---- cache + queue UI ----
    here      = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(here, 'gemsplots_queue.mat');

    state = buildFileQueue(cacheFile);
    if isempty(state.files)
        fprintf('No files selected. Exiting.\n');
        return;
    end

    %% ---- groups UI ----
    if ~isfield(state,'groups'), state.groups = {}; end
    conditions = uniqueStable(state.condition);
    [groups, ok] = uiDefineGroups(conditions, state.groups);
    if ~ok
        fprintf('Cancelled.\n');
        return;
    end
    state.groups = groups;
    fprintf('[plotAveragedMetrics] Saving groups to cache (%d group(s)).\n', ...
        numel(groups));
    save(cacheFile, '-struct', 'state');
    fprintf('[plotAveragedMetrics] Cache saved. About to call defineMetricsUI...\n');

    %% ---- metrics UI ----
    metricSpecs = defineMetricsUI(cacheFile, defaultMetricSpecs());
    fprintf('[plotAveragedMetrics] defineMetricsUI returned %d spec(s).\n', ...
        numel(metricSpecs));
    if isempty(metricSpecs)
        fprintf('No metrics defined. Exiting.\n');
        return;
    end

    %% ---- pre-load every metric value (with on-disk cache) ----
    % First run is slow (each loadMetric pulls a result file through
    % Google Drive); subsequent runs are instant because we cache the
    % extracted scalars to a small local .mat next to the queue cache.
    % Cache entries are invalidated automatically when the source file's
    % modification time changes. Delete gemsplots_metrics_cache.mat to
    % force a full reload.
    metricsCacheFile = fullfile(here, 'gemsplots_metrics_cache.mat');
    nMetrics = numel(metricSpecs);
    [values, hitCount, loadCount] = loadAllMetricsCached( ...
        state, metricSpecs, metricsCacheFile);
    fprintf( ...
        '[plotAveragedMetrics] Metric values ready: %d cache hit(s), %d fresh load(s).\n', ...
        hitCount, loadCount);

    %% ---- render ----
    try
        % Presentation-grade defaults — fonts and lines are larger so the
        % saved PNG/SVG read clearly when dropped into PowerPoint or
        % shrunk into a paper figure. Override per-call below if needed.
        pubfig_setup('Theme','light', ...
            'BaseFontSize', 24, ...
            'LineWidth',    2.4, ...
            'MarkerSize',   14);
    catch ME
        warning('pubfig_setup failed: %s', ME.message);
    end

    animalsAll = uniqueStable(state.animal);
    totalSteps = nMetrics * numel(groups);

    %% ---- ask for save folder (cancel = don't save) ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save every figure (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir, 0)
        saveDir = '';
        fprintf('[plotAveragedMetrics] Save cancelled — figures will be displayed but not written to disk.\n');
    else
        if ~exist(saveDir, 'dir'), mkdir(saveDir); end
        fprintf('[plotAveragedMetrics] Saving figures to: %s\n', saveDir);
    end

    %% ---- ask for auto-close delay ----
    % After each figure is rendered (and saved if a folder was picked),
    % pause this many seconds then close it. 0 = keep every figure open.
    defAutoClose = '5';
    if isempty(saveDir), defAutoClose = '0'; end   % nothing saved -> keep open by default
    ac = inputdlg( ...
        {sprintf(['Auto-close each figure after N seconds (0 = keep all open).\n' ...
                  '%d figures will be produced.'], nMetrics * numel(groups))}, ...
        'Auto-close', [1 60], {defAutoClose});
    if isempty(ac)
        autoCloseSec = 0;
    else
        autoCloseSec = str2double(ac{1});
        if isnan(autoCloseSec) || autoCloseSec < 0
            autoCloseSec = 0;
        end
    end
    if autoCloseSec > 0
        fprintf('[plotAveragedMetrics] Each figure will auto-close after %g s.\n', autoCloseSec);
    end

    wbR = waitbar(0, 'Rendering plots...', ...
        'Name','plotAveragedMetrics: rendering');
    try, set(wbR,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wbR,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>
    step = 0;
    try
        for k = 1:nMetrics
            spec = metricSpecs(k);
            for gi = 1:numel(groups)
                step = step + 1;
                if ishandle(wbR)
                    msg = sprintf('[%d / %d]  %s — group %d', ...
                        step, totalSteps, displayLabel(spec), gi);
                    if ~isempty(saveDir), msg = [msg ' (saving)']; end %#ok<AGROW>
                    waitbar(step/totalSteps, wbR, msg);
                end
                figH = renderMetricGroupFigure(state, animalsAll, groups{gi}, ...
                    spec, gi, values(:,k));
                if ~isempty(saveDir) && ~isempty(figH) && ishandle(figH)
                    saveFigureAllFormats(figH, saveDir);
                end
                if autoCloseSec > 0 && ~isempty(figH) && ishandle(figH)
                    drawnow;
                    pause(autoCloseSec);
                    if ishandle(figH), close(figH); end
                end
            end
        end
    catch ME
        if ishandle(wbR), close(wbR); end
        rethrow(ME);
    end
    if ishandle(wbR), close(wbR); end
    if ~isempty(saveDir)
        fprintf('[plotAveragedMetrics] Saved %d figure(s) to %s\n', step, saveDir);
    end

    fprintf('Rendered %d figures (%d metrics x %d groups).\n', ...
        step, nMetrics, numel(groups));
end


% ==========================================================================
% =========================== groups UI ===================================
% ==========================================================================
function [groups, ok] = uiDefineGroups(conditions, priorGroups)
% uifigure-based group definition. One row per group with an editable
% comma-separated condition list, pre-filled from priorGroups if any.
% Avoids the legacy inputdlg, which on macOS R2025a sometimes never
% surfaces when chained after a uifigure close.

    fprintf('[uiDefineGroups] Opening. Detected conditions: %s\n', ...
        strjoin(conditions, ', '));

    fig = uifigure('Name','Define groups', 'Position',[220 220 760 460]);
    g = uigridlayout(fig,[4 1]);
    g.RowHeight   = {'fit','fit','1x','fit'};
    g.ColumnWidth = {'1x'};

    uilabel(g, 'Text', ...
        sprintf('Detected conditions: %s', strjoin(conditions, ', ')), ...
        'WordWrap','on', 'FontWeight','bold');

    btnRow = uipanel(g,'BorderType','none');
    bl = uigridlayout(btnRow,[1 4]);
    bl.ColumnWidth = {110,140,'1x',110};
    bl.Padding = [0 4 0 4];
    btnAdd = uibutton(bl,'Text','Add group');
    btnRem = uibutton(bl,'Text','Remove selected');
    uilabel(bl,'Text','');
    btnGo  = uibutton(bl,'Text','Continue', ...
        'BackgroundColor',[0.7 0.9 0.7]);

    initData = priorGroupsToTable(priorGroups);
    if isempty(initData)
        initData = {'Group 1', ''};
    end
    tbl = uitable(g, ...
        'ColumnName',    {'Label','Conditions (comma-separated)'}, ...
        'ColumnEditable',[true true], ...
        'ColumnWidth',   {130,'1x'}, ...
        'Data',          initData, ...
        'CellSelectionCallback', @(s,e) setappdata(s,'sel',e.Indices));

    status = uilabel(g, 'Text', sprintf('%d group(s).', size(initData,1)));

    btnAdd.ButtonPushedFcn = @(~,~) onAdd();
    btnRem.ButtonPushedFcn = @(~,~) onRem();
    btnGo.ButtonPushedFcn  = @(~,~) onGo();
    fig.CloseRequestFcn    = @(~,~) onCancel();

    cancelled = false;
    uiwait(fig);

    if cancelled
        groups = {}; ok = false;
    else
        groups = tableToGroups(tbl.Data);
        ok = ~isempty(groups);
        if ~ok
            uialert(fig, 'All groups were empty.', 'Bad input');
        end
    end
    if isvalid(fig), delete(fig); end

    % ---- nested ----
    function onAdd()
        D = tbl.Data;
        D(end+1,:) = {sprintf('Group %d', size(D,1)+1), ''}; %#ok<AGROW>
        tbl.Data = D;
        status.Text = sprintf('%d group(s).', size(D,1));
    end
    function onRem()
        sel = getappdata(tbl,'sel');
        if isempty(sel), return; end
        D = tbl.Data;
        D(unique(sel(:,1)),:) = [];
        tbl.Data = D;
        status.Text = sprintf('%d group(s).', size(D,1));
    end
    function onGo()
        fprintf('[uiDefineGroups] Continue clicked. %d row(s).\n', ...
            size(tbl.Data,1));
        uiresume(fig);
    end
    function onCancel()
        fprintf('[uiDefineGroups] Cancel/close.\n');
        cancelled = true;
        uiresume(fig);
    end
end


function D = priorGroupsToTable(priorGroups)
    if isempty(priorGroups), D = cell(0,2); return; end
    n = numel(priorGroups);
    D = cell(n,2);
    for k = 1:n
        D{k,1} = sprintf('Group %d', k);
        if iscell(priorGroups{k})
            D{k,2} = strjoin(priorGroups{k}, ', ');
        else
            D{k,2} = char(priorGroups{k});
        end
    end
end

function groups = tableToGroups(D)
    groups = {};
    if isempty(D), return; end
    for i = 1:size(D,1)
        items = strsplit(strtrim(D{i,2}), ',');
        items = strtrim(items);
        items = items(~cellfun('isempty', items));
        if ~isempty(items)
            groups{end+1,1} = items; %#ok<AGROW>
        end
    end
end


% ==========================================================================
% =========================== rendering ===================================
% ==========================================================================
function figH = renderMetricGroupFigure(state, animalsAll, conds, spec, groupIdx, valuesForMetric)
% Render one figure for a single (metric, user-group). Conditions on x,
% baseline/recovery as the two subgroups. animalsAll fixes row ordering
% so the same animal lands at the same subjIdx in every figure.
% valuesForMetric is a length-nFiles vector in spec.unitsIn — it gets
% converted to spec.unitsOut here before being plotted. Returns the
% figure handle so callers can save it.

    figH = gobjects(0);
    nC = numel(conds);
    if nC == 0, return; end
    nAn = numel(animalsAll);

    % Convert raw stored values into the requested plot units once.
    unitsIn  = getf(spec,'unitsIn');
    unitsOut = getf(spec,'unitsOut');
    if ~isempty(unitsIn) && ~isempty(unitsOut) && ~strcmpi(unitsIn, unitsOut)
        valuesForMetric = convertUnits(valuesForMetric, unitsIn, unitsOut);
    end

    data = cell(nC, 2);                  % rows = conditions, cols = phases
    missing = cell(0,1);                 % flag queued sources whose value was NaN
    for ci = 1:nC
        cond = conds{ci};
        for ph = 1:2
            if ph == 1, phaseName = 'baseline'; else, phaseName = 'recovery'; end
            v = nan(nAn,1);
            for ai = 1:nAn
                aniName = animalsAll{ai};
                idx = find( strcmpi(state.condition, cond)     & ...
                            strcmpi(state.phase,     phaseName) & ...
                            strcmpi(state.animal,    aniName), 1);
                if ~isempty(idx)
                    v(ai) = valuesForMetric(idx);
                    if isnan(v(ai))
                        % Source was in queue but loader returned NaN —
                        % flag it. (Animals with no queued source for
                        % this combo are silently skipped — those are
                        % expected absences, not load failures.)
                        [~, fb, fe] = fileparts(state.files{idx});
                        missing{end+1,1} = sprintf( ...
                            '%s / %s / animal=%s  (%s)', ...
                            cond, phaseName, aniName, [fb fe]); %#ok<AGROW>
                    end
                end
            end
            data{ci, ph} = v;
        end
    end

    if ~isempty(missing)
        fprintf('[%s — group %d] %d missing data point(s):\n', ...
            displayLabel(spec), groupIdx, numel(missing));
        nShow = min(20, numel(missing));
        for n = 1:nShow
            fprintf('    %s\n', missing{n});
        end
        if numel(missing) > nShow
            fprintf('    ... (+%d more)\n', numel(missing) - nShow);
        end
    end

    dispLab = displayLabel(spec);
    figName = sprintf('group%d - %s', groupIdx, dispLab);
    % Force normal window state — pubfig_setup defaults to 'maximized'
    % which would stack 88 figures on top of each other. Size chosen
    % so the larger presentation-grade fonts breathe.
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position',[100 100 1100 700]);

    if all(cellfun(@(c) all(isnan(c(:))), data(:)))
        text(0.5, 0.5, sprintf('no data for group %d / %s', groupIdx, dispLab), ...
            'HorizontalAlignment','center','Units','normalized', ...
            'Interpreter','none');
        axis off;
        return;
    end

    boxScatterPlot(data, ...
        'GroupLabels',    conds, ...
        'SubgroupLabels', {'baseline','recovery'}, ...
        'YLabel',         dispLab, ...
        'XLabel',         sprintf('Group %d', groupIdx), ...
        'Title',          sprintf('%s — group %d', dispLab, groupIdx), ...
        'ColorBySubject', true, ...
        'ConnectPaired',  true, ...
        'PairLineStyle',  ':', ...
        'MarkerSize',     110);
end


% ==========================================================================
% =========================== defaults ====================================
% ==========================================================================
function specs = defaultMetricSpecs()
% Sensible defaults seeded into the metric table on first run.
    specs = [ ...
        ms('HR',             'bpm', 'bpm', '_HRBR.mat',         'avgHeartRate',       'auto', NaN); ...
        ms('Breathing rate', 'bpm', 'bpm', '_HRBR.mat',         'avgBreathRate',      'auto', NaN); ...
        ms('HRV',            's',   'ms',  '_HRVMeasures.mat',  'hrv',                'auto', NaN); ...
        ms('pNN5',           '%',   '%',   '_HRVMeasures.mat',  'pnn5',               'auto', NaN); ...
        ms('RMSSD',          's',   'ms',  '_HRVMeasures.mat',  'rmssd',              'auto', NaN); ...
        ms('Sample entropy', '',    '',    '_HRVMeasures.mat',  'sampEn',             'auto', NaN); ...
        ms('SD1',            's',   'ms',  '_HRVMeasures.mat',  'sd1',                'auto', NaN); ...
        ms('SD2',            's',   'ms',  '_HRVMeasures.mat',  'sd2',                'auto', NaN); ...
        ms('SW rate ANT1',   'cpm', 'cpm', '_slowWaves.mat',    'slowWaveRateSeries', 'mean', 1);   ...
        ms('SW rate ANT2',   'cpm', 'cpm', '_slowWaves.mat',    'slowWaveRateSeries', 'mean', 2);   ...
        ms('SW rate ANT3',   'cpm', 'cpm', '_slowWaves.mat',    'slowWaveRateSeries', 'mean', 3) ];
end

function s = ms(label, unitsIn, unitsOut, suffix, field, aggregator, channel)
    s = struct('label',label,'unitsIn',unitsIn,'unitsOut',unitsOut, ...
               'suffix',suffix,'field',field, ...
               'aggregator',aggregator,'channel',channel);
end

function s = displayLabel(spec)
% '<label> (<unitsOut>)' if unitsOut non-empty, else bare label.
    lab = strtrim(char(getf(spec,'label')));
    un  = strtrim(char(getf(spec,'unitsOut')));
    if isempty(un), s = lab; else, s = sprintf('%s (%s)', lab, un); end
end

function v = getf(s, name)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = ''; end
end


% ==========================================================================
% =========================== misc ========================================
% ==========================================================================
function saveFigureAllFormats(fig, outDir)
% Save one figure in .fig, .png, and .svg into outDir, using a sanitised
% version of the figure's Name as the basename. Each format is wrapped
% in its own try/catch so a single failure doesn't abort the loop.
    name = char(get(fig, 'Name'));
    if isempty(name), name = sprintf('figure_%d', fig.Number); end
    safe = regexprep(name, '[^\w\-.()% ]', '_');
    safe = regexprep(safe, '\s+', '_');
    safe = regexprep(safe, '_+', '_');
    safe = strtrim(safe);
    if isempty(safe), safe = sprintf('figure_%d', fig.Number); end
    base = fullfile(outDir, safe);

    try
        savefig(fig, [base '.fig']);
    catch ME
        warning('plotAveragedMetrics:saveFig', ...
            '.fig save failed for %s: %s', safe, ME.message);
    end
    try
        exportgraphics(fig, [base '.png'], 'Resolution', 200);
    catch ME
        warning('plotAveragedMetrics:savePng', ...
            '.png save failed for %s: %s', safe, ME.message);
    end
    try
        print(fig, [base '.svg'], '-dsvg', '-vector');
    catch ME
        warning('plotAveragedMetrics:saveSvg', ...
            '.svg save failed for %s: %s', safe, ME.message);
    end
end

function u = uniqueStable(c)
% unique() preserving the order of first occurrence; tolerates empty input.
    if isempty(c), u = {}; return; end
    if ~iscell(c), c = cellstr(c); end
    [~, ia] = unique(c, 'stable');
    u = c(sort(ia));
end

function s = shortName(fp)
% Last-folder/basename label for progress dialogs (no scary long paths).
    [d, b, e] = fileparts(char(fp));
    [~, parent] = fileparts(d);
    s = [parent filesep b e];
    if numel(s) > 70, s = ['...' s(end-66:end)]; end
end
