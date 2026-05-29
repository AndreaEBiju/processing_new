function plotTimeTraces()
%PLOTTIMETRACES  Per-(metric, condition) recovery-time-series plot with
% mean ± SEM band and individual animal traces.
%
% For each metric and each condition in the queue, find every recovery
% file matching that condition. Interpolate every animal's recovery
% time-series onto a common time grid (intersection range, median dt),
% then plot:
%   - shaded ± SEM band (theme color)
%   - light per-animal traces inside the band
%   - bold mean trace on top, with sample size in the legend
%
% One figure per (metric, condition). Saves to a user-picked folder if
% chosen, auto-closes after a user-defined delay.
%
% Caches the FULL time series via loadAllSeriesCached.m so subsequent
% runs are instant (separate cache from the windowed-scalar one).

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

    queueConds = unique(state.condition);
    allGroupConds = unique(horzcat(groups{:}));
    unknown = setdiff(allGroupConds, queueConds);
    if ~isempty(unknown)
        fprintf('[plotTimeTraces] !! UNKNOWN condition(s) in groups (not in queue):\n');
        for k = 1:numel(unknown), fprintf('     %s\n', unknown{k}); end
        fprintf('   Detected: %s\n', strjoin(queueConds, ', '));
    end

    %% ---- metrics (filter to series-capable, same pattern as plotWindowedMetrics) ----
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
        warndlg('No metric has both Series field and Time field.', 'plotTimeTraces');
        return;
    end
    metricSpecs = metricSpecs(hasSeries);
    fprintf('[plotTimeTraces] %d metric(s) x %d condition(s) -> %d figure(s) max.\n', ...
        numel(metricSpecs), numel(queueConds), numel(metricSpecs)*numel(queueConds));

    %% ---- output prefs ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save time-trace figures (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir,0), saveDir = '';
    else
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
    end

    nFigsExp = numel(metricSpecs) * numel(queueConds);
    defAutoClose = '5'; if isempty(saveDir), defAutoClose = '0'; end
    ac = inputdlg({sprintf(['Auto-close each figure after N seconds (0 = keep all open).\n' ...
        'Up to %d figures will be produced.'], nFigsExp)}, ...
        'Auto-close', [1 60], {defAutoClose});
    autoCloseSec = 0;
    if ~isempty(ac)
        v = str2double(ac{1});
        if ~isnan(v) && v >= 0, autoCloseSec = v; end
    end

    %% ---- load full series (cached + parallel) ----
    seriesCacheFile = fullfile(here, 'gemsplots_series_cache.mat');
    [seriesByFile, hitCount, loadCount] = ...
        loadAllSeriesCached(state, metricSpecs, seriesCacheFile);
    fprintf('[plotTimeTraces] Series ready: %d cache hit(s), %d fresh load(s).\n', ...
        hitCount, loadCount);

    %% ---- render ----
    try
        pubfig_setup('Theme','light', ...
            'BaseFontSize', 22, ...
            'LineWidth',    2.0, ...
            'MarkerSize',   10);
    catch ME
        warning('pubfig_setup failed: %s', ME.message);
    end

    totalSteps = numel(metricSpecs) * numel(queueConds);
    wbR = waitbar(0, 'Rendering time traces...', ...
        'Name','plotTimeTraces: rendering');
    try, set(wbR,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wbR,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>

    step = 0; nRendered = 0;
    try
        for k = 1:numel(metricSpecs)
            spec = metricSpecs(k);
            for ci = 1:numel(queueConds)
                cond = queueConds{ci};
                step = step + 1;
                if ishandle(wbR)
                    msg = sprintf('[%d / %d]  %s — %s', ...
                        step, totalSteps, displayLabel(spec), cond);
                    if ~isempty(saveDir), msg = [msg ' (saving)']; end %#ok<AGROW>
                    waitbar(step/totalSteps, wbR, msg);
                end
                figH = renderTraceFigure(state, spec, cond, seriesByFile(:,k));
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

    fprintf('Rendered %d time-trace figure(s) (skipped %d for missing data).\n', ...
        nRendered, totalSteps - nRendered);
end


% ==========================================================================
% =========================== render ======================================
% ==========================================================================
function figH = renderTraceFigure(state, spec, cond, seriesThisMetric)
    figH = gobjects(0);

    recIdx = find( strcmpi(state.condition, cond) & ...
                   strcmpi(state.phase,     'recovery') );
    if isempty(recIdx)
        fprintf('[plotTimeTraces] %s / %s : no recovery files in queue\n', ...
            displayLabel(spec), cond);
        return;
    end

    % Gather each animal's series
    Ys = {}; Ts = {}; animLabels = {};
    for ri = 1:numel(recIdx)
        i = recIdx(ri);
        s = seriesThisMetric{i};
        if isempty(s) || isempty(s.y) || isempty(s.t), continue; end
        Ys{end+1}         = double(s.y); %#ok<AGROW>
        Ts{end+1}         = double(s.t); %#ok<AGROW>
        animLabels{end+1} = state.animal{i}; %#ok<AGROW>
    end
    nTr = numel(Ys);
    if nTr == 0
        fprintf('[plotTimeTraces] %s / %s : no valid series loaded\n', ...
            displayLabel(spec), cond);
        return;
    end

    % Normalise each animal's time vector so it starts at zero (recovery
    % files often start at the cut point, not at 0).
    for k = 1:nTr
        Ts{k} = Ts{k} - Ts{k}(1);
    end

    % Common grid: from 0 to the SHORTEST animal's max (so every animal
    % contributes everywhere on the grid). Median dt across animals.
    tEnds = cellfun(@(t) t(end), Ts);
    tEnd = min(tEnds);
    if tEnd <= 0
        fprintf('[plotTimeTraces] %s / %s : zero common time range\n', ...
            displayLabel(spec), cond);
        return;
    end
    dts = arrayfun(@(k) median(diff(Ts{k})), 1:nTr);
    dt = median(dts(dts > 0));
    if ~isfinite(dt) || dt <= 0, dt = tEnd / 200; end   % fallback
    tGrid = (0:dt:tEnd).';

    Ymat = nan(numel(tGrid), nTr);
    for k = 1:nTr
        try
            Ymat(:,k) = interp1(Ts{k}, Ys{k}, tGrid, 'linear', NaN);
        catch
        end
    end

    nValid = sum(~isnan(Ymat), 2);
    meanY  = mean(Ymat, 2, 'omitnan');
    if nTr > 1
        sdY  = std(Ymat, 0, 2, 'omitnan');
        semY = sdY ./ sqrt(max(nValid, 1));
    else
        semY = zeros(size(meanY));
    end

    % ----- figure -----
    dispLab = displayLabel(spec);
    figName = sprintf('trace - %s - %s', cond, dispLab);
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position', [100 100 1100 600]);
    ax = axes('Parent', figH); hold(ax,'on');

    % Theme colour
    palette = getappdata(0, 'PUB_PALETTE');
    if isempty(palette) || size(palette,2) ~= 3
        col = [0 0.447 0.741];
    else
        col = palette(1,:);
    end

    % SEM band (only when n > 1)
    if nTr > 1
        keep = ~isnan(meanY) & ~isnan(semY);
        tBand = tGrid(keep);
        lo    = meanY(keep) - semY(keep);
        hi    = meanY(keep) + semY(keep);
        if numel(tBand) >= 2
            fill(ax, [tBand; flipud(tBand)], [lo; flipud(hi)], col, ...
                'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
                'HandleVisibility','off');
        end
    end

    % Individual traces (light)
    for k = 1:nTr
        plot(ax, tGrid, Ymat(:,k), '-', ...
            'Color', [col 0.40], 'LineWidth', 0.8, ...
            'HandleVisibility','off');
    end

    % Mean line (bold)
    plot(ax, tGrid, meanY, '-', ...
        'Color', col, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('mean (n=%d animals)', nTr));

    % Cosmetics
    xlabel(ax, 'Recovery time (s)', 'Interpreter','none');
    ylabel(ax, dispLab,             'Interpreter','none');
    title (ax, sprintf('%s — %s recovery', dispLab, cond), 'Interpreter','none');
    set(ax, 'TickLabelInterpreter','none');
    legend(ax, 'show', 'Location','best', 'Interpreter','none');
    grid(ax, 'on');
    xlim(ax, [0 tEnd]);

    % Footnote: which animals
    if ~isempty(animLabels)
        annotation(figH, 'textbox', [0.01 0.0 0.98 0.04], ...
            'String', sprintf('animals: %s', strjoin(unique(animLabels,'stable'), ', ')), ...
            'EdgeColor','none', 'FontSize', 10, ...
            'HorizontalAlignment','left', 'Interpreter','none');
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
        warning('plotTimeTraces:saveFig','.fig save failed: %s', ME.message);
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
        warning('plotTimeTraces:savePng','.png save failed: %s', ME.message);
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
            warning('plotTimeTraces:saveSvg','.svg save failed: %s', ME.message);
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
