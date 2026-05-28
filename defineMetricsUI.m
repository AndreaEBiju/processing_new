function specs = defineMetricsUI(cacheFile, defaults)
%DEFINEMETRICSUI  Interactive UI to define the set of metrics to plot.
%
%   specs = defineMetricsUI(cacheFile, defaults)
%     If cacheFile (.mat) has a `metrics` field, restore that table.
%     Otherwise pre-populate the UI with `defaults` (a struct array using
%     the same fields described below). On Continue, the resulting struct
%     array is saved back into cacheFile under the field `metrics`,
%     preserving any other fields already in the cache. Returns the spec
%     struct array.
%
%   specs = defineMetricsUI(cacheFile)
%     As above but with an empty default set (you'll start with no rows).
%
%   specs = defineMetricsUI('', defaults)
%     No persistence; show the UI seeded with `defaults`, return whatever
%     the user clicks through with.
%
% Each metric spec has fields
%   .label       string  — metric display name (without units)
%   .unitsIn     string  — unit the value is STORED in (e.g. 's' if the
%                          .mat saved HRV in seconds). Empty = "as-is".
%   .unitsOut    string  — unit to PLOT in (e.g. 'ms'). If unitsIn !=
%                          unitsOut, values are converted via
%                          convertUnits() at plot time. The plot y-axis
%                          reads '<label> (<unitsOut>)'.
%                          Leave both blank for unit-less metrics
%                          (e.g. sample entropy).
%   .suffix      string  — appended to the source file's stem to locate
%                          the output .mat (e.g. '_HRBR.mat'); '.mat' is
%                          auto-appended if missing
%   .field       string  — variable name inside the .mat
%   .aggregator  string  — one of {auto,mean,median,max,min,first,last,sum,scalar}
%                          applied after channel selection. 'auto' (default)
%                          returns the value as-is if already scalar, else
%                          takes the omitnan mean — so precomputed averages
%                          like avgHeartRate don't need a redundant op.
%   .channel     numeric — column index for multi-channel matrices;
%                          NaN/empty = no channel selection
%
% UI controls
%   Add row             - inserts an empty row at the bottom
%   Remove selected     - drops rows whose cells are currently selected
%   Reset to defaults   - replaces the table with `defaults`
%   Inspect .mat file…  - file picker, then lists the variables in that
%                         file in a dialog so you can see what fields
%                         exist before naming them in the table
%   Continue            - return the table contents

    if nargin < 1, cacheFile = ''; end
    if nargin < 2, defaults  = emptySpecs(); end

    fprintf('[defineMetricsUI] Entry. Cache=%s\n', cacheFile);

    % ----- restore -----
    specsIn = defaults;
    cacheState = struct();
    if ~isempty(cacheFile) && exist(cacheFile,'file')
        try
            cacheState = load(cacheFile);
            if isfield(cacheState,'metrics') && ~isempty(cacheState.metrics)
                specsIn = normaliseSpecs(cacheState.metrics);
            end
        catch ME
            warning('defineMetricsUI:cache', ...
                'Cache restore failed (%s); using defaults.', ME.message);
        end
    end

    % ----- run UI -----
    fprintf('[defineMetricsUI] Opening UI with %d existing spec(s).\n', ...
        numel(specsIn));
    specs = uiDefineMetrics(specsIn, defaults);
    fprintf('[defineMetricsUI] UI returned %d spec(s).\n', numel(specs));

    % ----- save -----
    if ~isempty(cacheFile) && ~isempty(specs)
        try
            cacheState.metrics = specs;       %#ok<NASGU>
            save(cacheFile, '-struct', 'cacheState');
            fprintf('[defineMetricsUI] Saved to %s\n', cacheFile);
        catch ME
            warning('defineMetricsUI:cache', ...
                'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
% =========================== UI ==========================================
% ==========================================================================
function specs = uiDefineMetrics(specsIn, defaults)
    fig = uifigure('Name','Define metrics', ...
        'Position',[200 200 1100 480]);
    grid = uigridlayout(fig,[3 1]);
    grid.RowHeight   = {'fit','1x','fit'};
    grid.ColumnWidth = {'1x'};

    % --- buttons ---
    btnRow = uipanel(grid,'BorderType','none');
    btnRow.Layout.Row = 1;
    btnLayout = uigridlayout(btnRow,[1 6]);
    btnLayout.ColumnWidth = {120,140,140,160,'1x',120};
    btnLayout.Padding = [0 4 0 4];
    btnAdd     = uibutton(btnLayout,'Text','Add row');
    btnRem     = uibutton(btnLayout,'Text','Remove selected');
    btnReset   = uibutton(btnLayout,'Text','Reset to defaults');
    btnInspect = uibutton(btnLayout,'Text','Inspect .mat file...');
    uilabel(btnLayout,'Text','');
    btnGo      = uibutton(btnLayout,'Text','Continue', ...
        'BackgroundColor',[0.7 0.9 0.7]);

    % --- table ---
    aggOpts = {'auto','mean','median','max','min','first','last','sum','scalar'};
    tbl = uitable(grid, ...
        'ColumnName',     {'Label','Stored units','Plot units', ...
                           'File suffix','Field','Aggregator','Channel'}, ...
        'ColumnEditable', [true true true true true true true], ...
        'ColumnFormat',   {'char','char','char','char','char',aggOpts,'numeric'}, ...
        'ColumnWidth',    {150,90,90,150,170,110,75}, ...
        'Data',           specsToTable(specsIn), ...
        'CellSelectionCallback', @(s,e) setappdata(s,'sel',e.Indices));
    tbl.Layout.Row = 2;

    % --- status ---
    status = uilabel(grid, ...
        'Text', sprintf('%d metric(s).', size(tbl.Data,1)));
    status.Layout.Row = 3;

    btnAdd.ButtonPushedFcn     = @(~,~) onAdd();
    btnRem.ButtonPushedFcn     = @(~,~) onRem();
    btnReset.ButtonPushedFcn   = @(~,~) onReset();
    btnInspect.ButtonPushedFcn = @(~,~) onInspect();
    btnGo.ButtonPushedFcn      = @(~,~) onGo();
    fig.CloseRequestFcn        = @(~,~) onCancel();

    cancelled = false;
    uiwait(fig);

    if cancelled
        specs = emptySpecs();
    else
        specs = tableToSpecs(tbl.Data);
    end
    if isvalid(fig), delete(fig); end

    % ---- nested ----
    function onAdd()
        D = tbl.Data;
        D(end+1,:) = {'', '', '', '', '', 'auto', NaN}; %#ok<AGROW>
        tbl.Data = D;
        status.Text = sprintf('%d metric(s).', size(D,1));
    end
    function onRem()
        sel = getappdata(tbl,'sel');
        if isempty(sel), return; end
        D = tbl.Data;
        D(unique(sel(:,1)),:) = [];
        tbl.Data = D;
        status.Text = sprintf('%d metric(s).', size(D,1));
    end
    function onReset()
        tbl.Data = specsToTable(defaults);
        status.Text = sprintf('%d metric(s).', size(tbl.Data,1));
    end
    function onInspect()
        [name, folder] = uigetfile({'*.mat','MAT files'}, 'Inspect a .mat file');
        if isequal(name,0), return; end
        fp = fullfile(folder,name);
        try
            w = whos('-file', fp);
        catch ME
            uialert(fig, sprintf('Could not read %s: %s', fp, ME.message), ...
                'Inspect failed');
            return;
        end
        if isempty(w)
            uialert(fig, sprintf('%s has no variables.', name), 'Inspect');
            return;
        end
        showVarsDialog(name, w);
    end
    function onGo()
        fprintf('[defineMetricsUI] Continue clicked. %d row(s).\n', ...
            size(tbl.Data,1));
        uiresume(fig);
    end
    function onCancel()
        fprintf('[defineMetricsUI] Cancel/close.\n');
        cancelled = true;
        uiresume(fig);
    end
end


% ==========================================================================
% =========================== conversions =================================
% ==========================================================================
function D = specsToTable(specs)
% Emit a 7-column cell table:
%   Label | StoredUnits | PlotUnits | Suffix | Field | Aggregator | Channel
    n = numel(specs);
    if n == 0, D = cell(0,7); return; end
    D = cell(n,7);
    for i = 1:n
        s = specs(i);
        D{i,1} = strOrEmpty(getf(s,'label'));
        D{i,2} = strOrEmpty(getf(s,'unitsIn'));
        D{i,3} = strOrEmpty(getf(s,'unitsOut'));
        D{i,4} = strOrEmpty(getf(s,'suffix'));
        D{i,5} = strOrEmpty(getf(s,'field'));
        D{i,6} = strOrEmpty(getf(s,'aggregator'));
        if isempty(D{i,6}), D{i,6} = 'auto'; end
        D{i,7} = numOrNaN(getf(s,'channel'));
    end
end

function specs = tableToSpecs(D)
    if isempty(D), specs = emptySpecs(); return; end
    rows = emptySpecs();
    for i = 1:size(D,1)
        lab = strOrEmpty(D{i,1});
        suf = strOrEmpty(D{i,4});
        fld = strOrEmpty(D{i,5});
        if isempty(lab) || isempty(suf) || isempty(fld)
            continue;       % skip incomplete rows
        end
        rows(end+1).label    = lab; %#ok<AGROW>
        rows(end).unitsIn    = strOrEmpty(D{i,2});
        rows(end).unitsOut   = strOrEmpty(D{i,3});
        rows(end).suffix     = suf;
        rows(end).field      = fld;
        agg = strOrEmpty(D{i,6});
        if isempty(agg), agg = 'auto'; end
        rows(end).aggregator = agg;
        rows(end).channel    = numOrNaN(D{i,7});
    end
    specs = rows(:);
end

function specs = normaliseSpecs(raw)
% Accept current + legacy shapes; emit a clean column struct array.
% Migration order (each fills any still-empty fields, never overwrites):
%   1. Convert legacy 5- or 6-col cell tables into struct array.
%   2. If 'units' field exists (older split-out form), copy to unitsOut.
%   3. If unitsOut is still empty, try splitting 'HRV (ms)' style labels.
%   4. If unitsIn is still empty, default it to unitsOut (no conversion).
    if isstruct(raw)
        specs = raw(:);
    elseif iscell(raw) && size(raw,2) == 7
        specs = tableToSpecs(raw);
    elseif iscell(raw) && size(raw,2) == 6
        specs = legacyTableToSpecs6(raw);
    elseif iscell(raw) && size(raw,2) == 5
        specs = legacyTableToSpecs5(raw);
    else
        specs = emptySpecs();
    end
    if isempty(specs), return; end

    % Add missing fields to the struct array (no-op if already there)
    if ~isfield(specs,'unitsIn'),    [specs.unitsIn]    = deal(''); end
    if ~isfield(specs,'unitsOut'),   [specs.unitsOut]   = deal(''); end
    if ~isfield(specs,'aggregator'), [specs.aggregator] = deal('auto'); end
    if ~isfield(specs,'channel'),    [specs.channel]    = deal(NaN); end

    for k = 1:numel(specs)
        if isempty(specs(k).unitsOut)
            if isfield(specs(k),'units') && ~isempty(specs(k).units)
                specs(k).unitsOut = specs(k).units;
            else
                [bare, u] = splitLabelUnits(specs(k).label);
                specs(k).label    = bare;
                specs(k).unitsOut = u;
            end
        end
        if isempty(specs(k).unitsIn)
            specs(k).unitsIn = specs(k).unitsOut;  % default: no conversion
        end
        if isempty(specs(k).aggregator)
            specs(k).aggregator = 'auto';
        end
    end
end

function specs = legacyTableToSpecs5(D)
% Very old 5-col table: Label | Suffix | Field | Aggregator | Channel
    rows = emptySpecs();
    for i = 1:size(D,1)
        lab = strOrEmpty(D{i,1});
        suf = strOrEmpty(D{i,2});
        fld = strOrEmpty(D{i,3});
        if isempty(lab) || isempty(suf) || isempty(fld), continue; end
        [bare, u] = splitLabelUnits(lab);
        rows(end+1).label    = bare; %#ok<AGROW>
        rows(end).unitsIn    = u;
        rows(end).unitsOut   = u;
        rows(end).suffix     = suf;
        rows(end).field      = fld;
        agg = strOrEmpty(D{i,4});
        if isempty(agg), agg = 'auto'; end
        rows(end).aggregator = agg;
        rows(end).channel    = numOrNaN(D{i,5});
    end
    specs = rows(:);
end

function specs = legacyTableToSpecs6(D)
% Earlier 6-col table: Label | Units | Suffix | Field | Aggregator | Channel
    rows = emptySpecs();
    for i = 1:size(D,1)
        lab = strOrEmpty(D{i,1});
        un  = strOrEmpty(D{i,2});
        suf = strOrEmpty(D{i,3});
        fld = strOrEmpty(D{i,4});
        if isempty(lab) || isempty(suf) || isempty(fld), continue; end
        rows(end+1).label    = lab; %#ok<AGROW>
        rows(end).unitsIn    = un;
        rows(end).unitsOut   = un;
        rows(end).suffix     = suf;
        rows(end).field      = fld;
        agg = strOrEmpty(D{i,5});
        if isempty(agg), agg = 'auto'; end
        rows(end).aggregator = agg;
        rows(end).channel    = numOrNaN(D{i,6});
    end
    specs = rows(:);
end

function [bare, units] = splitLabelUnits(label)
% 'HRV (ms)' -> ('HRV', 'ms');  'Sample entropy' -> ('Sample entropy', '')
    label = char(label);
    m = regexp(label, '^(.*?)\s*\(([^)]+)\)\s*$', 'tokens', 'once');
    if ~isempty(m) && numel(m) == 2
        bare  = strtrim(m{1});
        units = strtrim(m{2});
    else
        bare  = strtrim(label);
        units = '';
    end
end

function s = emptySpecs()
    s = struct('label',{},'unitsIn',{},'unitsOut',{},'suffix',{},'field',{}, ...
               'aggregator',{},'channel',{});
end


% ==========================================================================
% =========================== misc ========================================
% ==========================================================================
function showVarsDialog(filename, w)
% Modeless uifigure showing every variable in a .mat file as a sortable,
% scrollable uitable. Sorting (toggle headers) makes long lists easier
% to scan; clicking a row copies that variable name to the clipboard
% so it can be pasted straight into the Field column of the main table.
    f = uifigure('Name', sprintf('Variables in %s', filename), ...
        'Position', [200 200 720 540]);
    g = uigridlayout(f, [2 1]);
    g.RowHeight   = {'1x','fit'};
    g.ColumnWidth = {'1x'};

    rows = cell(numel(w), 4);
    for k = 1:numel(w)
        rows{k,1} = w(k).name;
        rows{k,2} = w(k).class;
        rows{k,3} = mat2str(w(k).size);
        rows{k,4} = w(k).bytes;
    end
    tbl = uitable(g, ...
        'ColumnName',     {'Name','Class','Size','Bytes'}, ...
        'ColumnWidth',    {260, 100, 160, 100}, ...
        'ColumnSortable', [true true true true], ...
        'Data',           rows, ...
        'CellSelectionCallback', @(~,e) copyName(rows, e));
    tbl.Layout.Row = 1;  tbl.Layout.Column = 1;

    btnRow = uipanel(g,'BorderType','none');
    btnRow.Layout.Row = 2;
    bl = uigridlayout(btnRow,[1 2]);
    bl.ColumnWidth = {'1x','fit'};
    bl.Padding = [0 4 0 4];
    info = uilabel(bl, ...
        'Text', sprintf('%d variable(s). Click a row to copy the name.', numel(w)));
    info.Layout.Column = 1;
    btnClose = uibutton(bl,'Text','Close','ButtonPushedFcn', @(~,~) delete(f));
    btnClose.Layout.Column = 2;

    function copyName(rows, e)
        if isempty(e.Indices), return; end
        r = e.Indices(1,1);
        if r >= 1 && r <= size(rows,1)
            try
                clipboard('copy', rows{r,1});
                info.Text = sprintf('Copied "%s" to clipboard.', rows{r,1});
            catch
                % clipboard unavailable on some platforms — silent
            end
        end
    end
end


function v = getf(s, name)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = []; end
end
function s = strOrEmpty(v)
    if isstring(v), s = char(v);
    elseif ischar(v), s = v;
    elseif isnumeric(v) && isscalar(v) && ~isnan(v), s = num2str(v);
    else, s = '';
    end
    s = strtrim(s);
end
function n = numOrNaN(v)
    if isempty(v), n = NaN;
    elseif isnumeric(v) && isscalar(v), n = double(v);
    elseif ischar(v) || isstring(v)
        n = str2double(v);
        if isempty(n), n = NaN; end
    else, n = NaN;
    end
end
