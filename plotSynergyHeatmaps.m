function plotSynergyHeatmaps()
%PLOTSYNERGYHEATMAPS  Synergy heatmaps for each metric, using our
% baseline-normalised effect (recovery − baseline)/baseline (NOT the raw
% Δ used in the Dashboard_Guide doc).
%
% For each metric and each combined (M, E) condition (3 × 3 = 9 cells):
%   norm(cond)    = mean over animals of (rec_mean − base_mean)/base_mean
%   synergy(E,M)  = norm(M,E)  −  [ norm(M-alone) + norm(E-alone) ]
%
% Working on baseline-normalised values makes control_norm = 0 (since
% (base − base)/base = 0), so the standard "+ control" term in the
% additive expectation drops out and the formula reduces to the one
% above (same reasoning the Dashboard_Guide gives for raw Δ).
%
% Cell colour: diverging palette centred at 0
%   red  = super-additive (combined > sum of parts)
%   blue = sub-additive  (combined < sum of parts)
%   white ≈ additive    (combined ≈ E-alone + M-alone)
% Numeric value annotated in each cell. NaN cells are greyed out and
% labelled 'NaN'.
%
% Sharing the same series cache as plotWindowedMetrics /
% plotWindowedViolins / plotTimeTraces, so if the data has already been
% loaded by any of them this runs in seconds.

    %% ---- queue + metrics ----
    here      = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(here, 'gemsplots_queue.mat');

    state = buildFileQueue(cacheFile);
    if isempty(state.files), fprintf('No files. Exiting.\n'); return; end

    metricSpecs = defineMetricsUI(cacheFile, defaultWindowedMetricSpecs());
    if isempty(metricSpecs), fprintf('No metrics. Exiting.\n'); return; end

    % Auto-fill series/time fields from label defaults (same as other wrappers)
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
        warndlg('No metric has both Series field and Time field.', 'plotSynergyHeatmaps');
        return;
    end
    metricSpecs = metricSpecs(hasSeries);

    %% ---- parse (M, E) levels from condition strings ----
    [Mlev, Elev] = parseStimLevels(state.condition);
    if any(isnan(Mlev) | isnan(Elev))
        bad = state.condition(isnan(Mlev) | isnan(Elev));
        bad = uniqueStable(bad);
        fprintf('[plotSynergyHeatmaps] !! Could not parse (M, E) from these conditions (skipped):\n');
        for k = 1:numel(bad), fprintf('     %s\n', bad{k}); end
    end

    % E-alone conditions have Mlev=0; M-alone have Elev=0; combined have both > 0
    eLevels = sort(unique(Elev(Elev > 0 & Mlev == 0 & ~isnan(Elev))));
    mLevels = sort(unique(Mlev(Mlev > 0 & Elev == 0 & ~isnan(Mlev))));
    if isempty(eLevels) || isempty(mLevels)
        warndlg(['Need both E-alone and M-alone conditions in the queue to ' ...
                 'compute synergy. Found ' num2str(numel(eLevels)) ' E-alone, ' ...
                 num2str(numel(mLevels)) ' M-alone.'], 'plotSynergyHeatmaps');
        return;
    end
    fprintf('[plotSynergyHeatmaps] E levels: %s   M levels: %s\n', ...
        mat2str(eLevels), mat2str(mLevels));

    %% ---- save folder + auto-close ----
    saveDir = uigetdir(pwd, ...
        'Pick a folder to save synergy heatmaps (.fig/.png/.svg). Cancel = display only.');
    if isequal(saveDir,0), saveDir = '';
    else
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
    end

    defAutoClose = '5'; if isempty(saveDir), defAutoClose = '0'; end
    ac = inputdlg({sprintf(['Auto-close each figure after N seconds (0 = keep open).\n' ...
        '%d figures will be produced.'], numel(metricSpecs))}, ...
        'Auto-close', [1 60], {defAutoClose});
    autoCloseSec = 0;
    if ~isempty(ac)
        v = str2double(ac{1});
        if ~isnan(v) && v >= 0, autoCloseSec = v; end
    end

    %% ---- load series (shared cache + parallel) ----
    seriesCacheFile = fullfile(here, 'gemsplots_series_cache.mat');
    [seriesByFile, hitCount, loadCount] = ...
        loadAllSeriesCached(state, metricSpecs, seriesCacheFile);
    fprintf('[plotSynergyHeatmaps] Series ready: %d cache hit(s), %d fresh load(s).\n', ...
        hitCount, loadCount);

    %% ---- compute baseline-normalised values per (file) ----
    % For each (animal, condition) pair (i.e. one baseline + one recovery
    % file), compute (rec_mean - base_mean) / base_mean and keep per condition.
    [normByCond, nByCond] = computeNormalisedByCondition( ...
        state, metricSpecs, seriesByFile);

    %% ---- render ----
    try
        pubfig_setup('Theme','light', ...
            'BaseFontSize', 22, ...
            'LineWidth',    2.0, ...
            'MarkerSize',   10);
    catch ME
        warning('pubfig_setup failed: %s', ME.message);
    end

    wbR = waitbar(0, 'Rendering synergy heatmaps...', ...
        'Name','plotSynergyHeatmaps: rendering');
    try, set(wbR,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wbR,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>

    nRendered = 0;
    try
        for k = 1:numel(metricSpecs)
            spec = metricSpecs(k);
            if ishandle(wbR)
                msg = sprintf('[%d / %d]  %s', k, numel(metricSpecs), displayLabel(spec));
                if ~isempty(saveDir), msg = [msg ' (saving)']; end %#ok<AGROW>
                waitbar(k/numel(metricSpecs), wbR, msg);
            end
            figH = renderSynergyHeatmap(spec, normByCond(k,:), nByCond(k,:), ...
                state.condition, mLevels, eLevels);
            if isempty(figH) || ~ishandle(figH), continue; end
            nRendered = nRendered + 1;
            if ~isempty(saveDir), saveFigureAllFormats(figH, saveDir); end
            if autoCloseSec > 0
                drawnow; pause(autoCloseSec);
                if ishandle(figH), close(figH); end
            end
        end
    catch ME
        if ishandle(wbR), close(wbR); end
        rethrow(ME);
    end
    if ishandle(wbR), close(wbR); end

    fprintf('Rendered %d synergy heatmap(s).\n', nRendered);
end


% ==========================================================================
% =========================== parse + aggregate ============================
% ==========================================================================
function [Mlev, Elev] = parseStimLevels(conditions)
% Parse condition strings into numeric M and E stimulation levels.
%   'M10'       -> M=10, E=0
%   'E100'      -> E=100, M=0
%   'M10E100'   -> M=10, E=100
%   'E100M10'   -> E=100, M=10
%   anything else -> NaN
    n = numel(conditions);
    Mlev = nan(n,1);
    Elev = nan(n,1);
    for k = 1:n
        cond = upper(strtrim(char(conditions{k})));
        m = regexp(cond, '^M(\d+)E(\d+)$', 'tokens', 'once');
        if ~isempty(m)
            Mlev(k) = str2double(m{1});
            Elev(k) = str2double(m{2});
            continue;
        end
        m = regexp(cond, '^E(\d+)M(\d+)$', 'tokens', 'once');
        if ~isempty(m)
            Elev(k) = str2double(m{1});
            Mlev(k) = str2double(m{2});
            continue;
        end
        m = regexp(cond, '^M(\d+)$', 'tokens', 'once');
        if ~isempty(m)
            Mlev(k) = str2double(m{1});
            Elev(k) = 0;
            continue;
        end
        m = regexp(cond, '^E(\d+)$', 'tokens', 'once');
        if ~isempty(m)
            Elev(k) = str2double(m{1});
            Mlev(k) = 0;
        end
    end
end


function [normByCond, nByCond] = computeNormalisedByCondition(state, metricSpecs, seriesByFile)
% For each metric and each (animal, condition) pair, compute
% (rec_mean - base_mean) / base_mean using the FULL-series means, then
% return a per-condition mean across animals.
%
% Returns
%   normByCond : nMetrics x nConds  (mean normalised value per cond)
%   nByCond    : nMetrics x nConds  (animal count contributing)
% Conditions are state-condition's unique list (uniqueStable order).
%
% NOTE: if a single animal contributes multiple trials of the same
% condition (multiple baseline AND recovery files), each trial's
% normalised value is averaged in.

    nFiles    = numel(state.files);
    nMetrics  = numel(metricSpecs);
    conds     = state.condition;
    phases    = state.phase;
    animals   = state.animal;

    % file-level mean of the full series
    fileMean = nan(nFiles, nMetrics);
    for i = 1:nFiles
        for k = 1:nMetrics
            s = seriesByFile{i,k};
            if isempty(s) || isempty(s.y), continue; end
            fileMean(i,k) = mean(s.y, 'omitnan');
        end
    end

    uconds = uniqueStable(conds);
    normByCond = nan(nMetrics, numel(uconds));
    nByCond    = zeros(nMetrics, numel(uconds));

    for c = 1:numel(uconds)
        cond = uconds{c};
        % every (animal, trial) pair in this condition
        baseRows = find(strcmpi(conds, cond) & strcmpi(phases, 'baseline'));
        recRows  = find(strcmpi(conds, cond) & strcmpi(phases, 'recovery'));
        if isempty(baseRows) || isempty(recRows), continue; end

        % match each recovery row to an animal-matched baseline row
        for rr = 1:numel(recRows)
            iRec = recRows(rr);
            iBase = baseRows(find(strcmpi(animals(baseRows), animals{iRec}), 1));
            if isempty(iBase), continue; end
            for k = 1:nMetrics
                bm = fileMean(iBase,k);
                rm = fileMean(iRec,k);
                if isnan(bm) || isnan(rm) || bm == 0, continue; end
                v = (rm - bm) / bm;
                if isnan(normByCond(k,c))
                    normByCond(k,c) = v;
                    nByCond(k,c)    = 1;
                else
                    n = nByCond(k,c);
                    normByCond(k,c) = (normByCond(k,c)*n + v) / (n+1);
                    nByCond(k,c)    = n + 1;
                end
            end
        end
    end
end


% ==========================================================================
% =========================== render =======================================
% ==========================================================================
function figH = renderSynergyHeatmap(spec, normRow, nRow, allConds, mLevels, eLevels)
% normRow / nRow are 1 x nConds slices (one metric) of the matrices
% returned by computeNormalisedByCondition.

    figH = gobjects(0);
    uconds = uniqueStable(allConds);
    nM = numel(mLevels);
    nE = numel(eLevels);

    % index lookup: condition string -> position in uconds
    idxFor = @(cstr) find(strcmpi(uconds, cstr), 1);

    synergy   = nan(nM, nE);
    observed  = nan(nM, nE);
    expected  = nan(nM, nE);
    nCounts   = zeros(nM, nE);   % min animal count across (obs, mAlone, eAlone)

    for i = 1:nM
        for j = 1:nE
            combName  = sprintf('M%dE%d', mLevels(i), eLevels(j));
            mAlone    = sprintf('M%d',    mLevels(i));
            eAlone    = sprintf('E%d',    eLevels(j));
            ic = idxFor(combName);  im = idxFor(mAlone);  ie = idxFor(eAlone);
            if isempty(ic) || isempty(im) || isempty(ie), continue; end

            obs = normRow(ic);  mA = normRow(im);  eA = normRow(ie);
            if isnan(obs) || isnan(mA) || isnan(eA), continue; end

            observed(i,j) = obs;
            expected(i,j) = mA + eA;
            synergy(i,j)  = obs - (mA + eA);
            nCounts(i,j)  = min([nRow(ic), nRow(im), nRow(ie)]);
        end
    end

    dispLab = displayLabel(spec);
    figName = sprintf('synergy - %s', dispLab);
    figH = figure('Name', figName, 'WindowState','normal', ...
        'Position', [100 100 1000 750]);
    ax = axes('Parent', figH);

    % symmetric range for diverging colormap
    maxAbs = max(abs(synergy(:)), [], 'omitnan');
    if isempty(maxAbs) || ~isfinite(maxAbs) || maxAbs == 0, maxAbs = 1; end

    % Use a NaN-friendly imagesc: replace NaN with sentinel for plotting,
    % then mask the colour. Simplest: imagesc with AlphaData.
    plotMat = synergy;
    plotMat(isnan(plotMat)) = 0;
    h = imagesc(ax, plotMat, [-maxAbs, maxAbs]);
    set(h, 'AlphaData', ~isnan(synergy));

    colormap(ax, rdbuColormap(64));
    cb = colorbar(ax);
    cb.Label.String       = 'synergy = norm(combined) − [norm(M-alone) + norm(E-alone)]';
    cb.Label.Interpreter  = 'none';
    cb.TickLabelInterpreter = 'none';

    set(ax, ...
        'XTick',                 1:nE, ...
        'XTickLabel',            arrayfun(@(x) sprintf('E%d', x), eLevels, 'UniformOutput', false), ...
        'YTick',                 1:nM, ...
        'YTickLabel',            arrayfun(@(x) sprintf('M%d', x), mLevels, 'UniformOutput', false), ...
        'TickLabelInterpreter',  'none', ...
        'YDir',                  'reverse', ...
        'XAxisLocation',         'top', ...
        'Color',                 [0.92 0.92 0.92]);   % background for NaN cells
    axis(ax, 'equal'); axis(ax, 'tight');

    xlabel(ax, 'Electrical level (E)', 'Interpreter','none');
    ylabel(ax, 'Mechanical level (M)', 'Interpreter','none');
    title(ax, sprintf('%s — synergy heatmap  (normalised: (rec − base) / base)', dispLab), ...
        'Interpreter','none');

    % cell annotations: synergy value + (n = ...)
    for i = 1:nM
        for j = 1:nE
            if isnan(synergy(i,j))
                text(ax, j, i, 'NaN', ...
                    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                    'Color', [0.45 0.45 0.45], 'FontSize', 16);
                continue;
            end
            % text colour: dark for pale cells, white for saturated
            relMag = abs(synergy(i,j)) / maxAbs;
            if relMag > 0.55, txtCol = [1 1 1]; else, txtCol = [0 0 0]; end
            text(ax, j, i, sprintf('%+.3f\n(n=%d)', synergy(i,j), nCounts(i,j)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'Color', txtCol, 'FontSize', 16, 'Interpreter','none');
        end
    end
end


% ==========================================================================
% =========================== helpers ======================================
% ==========================================================================
function cm = rdbuColormap(n)
% Diverging red-white-blue colormap. Mid colour is white.
    if nargin < 1, n = 64; end
    half = floor(n/2);
    % blue end -> white middle
    rB = linspace(0.05, 1.00, half).';
    gB = linspace(0.27, 1.00, half).';
    bB = linspace(0.55, 1.00, half).';
    % white middle -> red end
    rR = linspace(1.00, 0.65, n - half).';
    gR = linspace(1.00, 0.10, n - half).';
    bR = linspace(1.00, 0.13, n - half).';
    cm = [rB gB bB; rR gR bR];
end


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
        warning('plotSynergyHeatmaps:saveFig','.fig save failed: %s', ME.message);
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
        warning('plotSynergyHeatmaps:savePng','.png save failed: %s', ME.message);
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
            warning('plotSynergyHeatmaps:saveSvg','.svg save failed: %s', ME.message);
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
