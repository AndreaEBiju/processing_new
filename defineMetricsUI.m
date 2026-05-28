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
%   .label       string  — y-axis label / display name
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
        'ColumnName',     {'Label','File suffix','Field','Aggregator','Channel'}, ...
        'ColumnEditable', [true true true true true], ...
        'ColumnFormat',   {'char','char','char',aggOpts,'numeric'}, ...
        'ColumnWidth',    {220,180,200,120,90}, ...
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
        D(end+1,:) = {'',  '',  '',  'auto',  NaN}; %#ok<AGROW>
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
    n = numel(specs);
    if n == 0, D = cell(0,5); return; end
    D = cell(n,5);
    for i = 1:n
        s = specs(i);
        D{i,1} = strOrEmpty(getf(s,'label'));
        D{i,2} = strOrEmpty(getf(s,'suffix'));
        D{i,3} = strOrEmpty(getf(s,'field'));
        D{i,4} = strOrEmpty(getf(s,'aggregator'));
        if isempty(D{i,4}), D{i,4} = 'auto'; end
        D{i,5} = numOrNaN(getf(s,'channel'));
    end
end

function specs = tableToSpecs(D)
    if isempty(D), specs = emptySpecs(); return; end
    rows = struct('label',{},'suffix',{},'field',{}, ...
                  'aggregator',{},'channel',{});
    for i = 1:size(D,1)
        lab = strOrEmpty(D{i,1});
        suf = strOrEmpty(D{i,2});
        fld = strOrEmpty(D{i,3});
        if isempty(lab) || isempty(suf) || isempty(fld)
            continue;       % skip incomplete rows
        end
        rows(end+1).label    = lab; %#ok<AGROW>
        rows(end).suffix     = suf;
        rows(end).field      = fld;
        agg = strOrEmpty(D{i,4});
        if isempty(agg), agg = 'auto'; end
        rows(end).aggregator = agg;
        rows(end).channel    = numOrNaN(D{i,5});
    end
    specs = rows(:);
end

function specs = normaliseSpecs(raw)
% Accept various legacy shapes; emit a clean column struct array.
    if isstruct(raw)
        specs = raw(:);
    elseif iscell(raw) && size(raw,2) >= 5
        specs = tableToSpecs(raw);
    else
        specs = emptySpecs();
    end
    % fill any missing fields
    for k = 1:numel(specs)
        if ~isfield(specs,'aggregator') || isempty(specs(k).aggregator)
            specs(k).aggregator = 'auto';
        end
        if ~isfield(specs,'channel')
            specs(k).channel = NaN;
        end
    end
end

function s = emptySpecs()
    s = struct('label',{},'suffix',{},'field',{},'aggregator',{},'channel',{});
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
