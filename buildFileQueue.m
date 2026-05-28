function state = buildFileQueue(cacheFile)
%BUILDFILEQUEUE  Interactive UI to build/edit a persistent file queue.
%
%   state = buildFileQueue(cacheFile)
%     Restores any saved queue from cacheFile (.mat), shows a UI to add /
%     remove / edit files, then saves the updated queue back to cacheFile
%     on Continue. cacheFile is created if it doesn't yet exist.
%
%   state = buildFileQueue()
%     No persistence — starts with an empty queue, returns the entered
%     queue without saving.
%
%   Returns struct with fields
%     .files     cell column, full source-file paths
%     .animal    cell column, one-letter animal IDs
%     .condition cell column, stim-condition strings (e.g. 'M10')
%     .phase     cell column, 'baseline' or 'recovery'
%   Any other fields previously stored in the cache (e.g. 'groups') are
%   preserved and returned untouched, so callers can layer extra metadata.
%
%   If the user cancels (closes the window), state.files is {}.
%
%   Auto-parse heuristics for each added file
%     animal    = first letter of the 2nd underscore-token of the basename
%     condition = first token that matches /^[ME]\d+([ME]\d+)?$/ (e.g.
%                 M10, M50, M10E100); falls back to token 1
%     phase     = 'recovery' if the basename contains 'recovery', else
%                 'baseline' (so a missing label is fine to leave as-is
%                 if the file really is baseline)
%   All four columns are editable inline in the table — fix anything the
%   parser got wrong before clicking Continue.
%
%   UI buttons
%     Add files...     standard file picker; multi-select inside one folder
%     Add folder...    pick a root folder + glob pattern (default
%                      '*_v*_blankmotion.mat'), then recursively walks the
%                      tree and adds every match. Duplicates against the
%                      existing queue are skipped automatically.
%     Remove selected  drop the rows whose cells are currently selected
%     Clear queue      empty the table
%     Continue         return the table contents
%
%   Version-filter rule
%     Both Add files... and Add folder... only accept basenames that
%     contain a semver-style version tag of the form _v<major>(.<minor>)+_
%     (e.g. _v0.1.0_, _v1.2_, _v0.10.3_). Files like
%     'foo_blankmotion.mat' or 'foo_HRBR.mat' without a version tag are
%     silently dropped; the status bar shows how many were skipped. To
%     change the rule, edit `isVersionedName` at the bottom of this file.

    if nargin < 1, cacheFile = ''; end

    %% restore cache (or start empty)
    state = struct('files',{{}}, 'animal',{{}}, 'condition',{{}}, ...
                   'phase',{{}});
    if ~isempty(cacheFile) && exist(cacheFile,'file')
        try
            S = load(cacheFile);
            f = fieldnames(S);
            for k = 1:numel(f)
                state.(f{k}) = S.(f{k});      % including extras like .groups
            end
        catch ME
            warning('buildFileQueue:cache', ...
                'Cache restore failed (%s); starting fresh.', ME.message);
        end
    end

    %% run UI
    state = uiBuildQueue(state);

    %% save back on success
    if ~isempty(state.files) && ~isempty(cacheFile)
        try
            save(cacheFile, '-struct', 'state');
        catch ME
            warning('buildFileQueue:cache', ...
                'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
% =========================== UI ==========================================
% ==========================================================================
function state = uiBuildQueue(state)
    fig = uifigure('Name','Build file queue', ...
        'Position',[200 200 1100 560]);
    grid = uigridlayout(fig,[3 1]);
    grid.RowHeight   = {'fit','1x','fit'};
    grid.ColumnWidth = {'1x'};

    % --- buttons row ---
    btnRow = uipanel(grid,'BorderType','none');
    btnRow.Layout.Row = 1; btnRow.Layout.Column = 1;
    btnLayout = uigridlayout(btnRow,[1 6]);
    btnLayout.ColumnWidth = {120,150,140,120,'1x',120};
    btnLayout.Padding = [0 4 0 4];
    btnAdd     = uibutton(btnLayout,'Text','Add files...');
    btnAddDir  = uibutton(btnLayout,'Text','Add folder...');
    btnRem     = uibutton(btnLayout,'Text','Remove selected');
    btnClr     = uibutton(btnLayout,'Text','Clear queue');
    uilabel(btnLayout,'Text','');
    btnGo      = uibutton(btnLayout,'Text','Continue', ...
        'BackgroundColor',[0.7 0.9 0.7]);

    % --- table ---
    initData = stateToTable(state);
    tbl = uitable(grid, ...
        'ColumnName',     {'File','Animal','Condition','Phase'}, ...
        'ColumnEditable', [false true true true], ...
        'ColumnWidth',    {'auto',80,120,100}, ...
        'Data',           initData, ...
        'CellSelectionCallback', @(s,e) setappdata(s,'sel',e.Indices));
    tbl.Layout.Row = 2;  tbl.Layout.Column = 1;

    % --- status ---
    status = uilabel(grid, ...
        'Text', sprintf('%d file(s) queued.', size(initData,1)));
    status.Layout.Row = 3;  status.Layout.Column = 1;

    btnAdd.ButtonPushedFcn    = @(~,~) onAdd();
    btnAddDir.ButtonPushedFcn = @(~,~) onAddFolder();
    btnRem.ButtonPushedFcn    = @(~,~) onRem();
    btnClr.ButtonPushedFcn    = @(~,~) onClear();
    btnGo.ButtonPushedFcn     = @(~,~) onGo();
    fig.CloseRequestFcn       = @(~,~) onCancel();

    cancelled = false;
    uiwait(fig);

    if cancelled
        state.files = {}; state.animal = {};
        state.condition = {}; state.phase = {};
    else
        D = tbl.Data;
        if isempty(D)
            state.files = {}; state.animal = {};
            state.condition = {}; state.phase = {};
        else
            state.files     = D(:,1);
            state.animal    = D(:,2);
            state.condition = D(:,3);
            state.phase     = D(:,4);
        end
    end
    if isvalid(fig), delete(fig); end

    % ---- nested ----
    function onAdd()
        [names, folder] = uigetfile({'*.mat','MAT files'}, ...
            'Add files to queue','MultiSelect','on');
        if isequal(names,0), return; end
        if ischar(names), names = {names}; end
        D = tbl.Data;
        nAdded = 0; nSkippedVer = 0; nSkippedDup = 0;
        for k = 1:numel(names)
            fpath = fullfile(folder, names{k});
            if ~isVersionedName(names{k})
                nSkippedVer = nSkippedVer + 1;
                continue;
            end
            if ~isempty(D) && any(strcmpi(D(:,1), fpath))
                nSkippedDup = nSkippedDup + 1;
                continue;
            end
            [ani, cond, ph] = parseFilename(names{k});
            D(end+1,:) = {fpath, ani, cond, ph}; %#ok<AGROW>
            nAdded = nAdded + 1;
        end
        tbl.Data = D;
        status.Text = sprintf( ...
            '%d file(s) queued. Added %d, skipped %d unversioned, %d duplicate(s).', ...
            size(D,1), nAdded, nSkippedVer, nSkippedDup);
        drawnow;
    end
    function onAddFolder()
        rootDir = uigetdir(pwd, 'Pick a folder to search recursively');
        if isequal(rootDir, 0), return; end
        p = inputdlg({['Filename pattern (glob).  Default keeps only ' ...
                      'versioned files of the form *_v0.x.x_blankmotion.mat:']}, ...
            'Recursive search', [1 80], {'*_v*_blankmotion.mat'});
        if isempty(p), return; end
        pattern = strtrim(p{1});
        if isempty(pattern), pattern = '*_v*_blankmotion.mat'; end

        % --- progress dialog so the long recursive dir() is visible ---
        pd = uiprogressdlg(fig, ...
            'Title',       'Scanning folder', ...
            'Message',     sprintf('Scanning %s for "%s" ...', rootDir, pattern), ...
            'Indeterminate','on', ...
            'Cancelable',  'off');
        drawnow;

        try
            hits = dir(fullfile(rootDir, '**', pattern));
            hits = hits(~[hits.isdir]);

            if isempty(hits)
                safeClose(pd);
                uialert(fig, sprintf('No files matching "%s" under %s', ...
                    pattern, rootDir), 'Recursive search', 'Icon','warning');
                return;
            end

            % switch to determinate progress for the per-file pass
            pd.Indeterminate = 'off';
            pd.Title         = 'Adding files';

            D = tbl.Data;
            nAdded = 0; nSkippedVer = 0; nSkippedDup = 0;
            for k = 1:numel(hits)
                pd.Value   = k / numel(hits);
                pd.Message = sprintf('[%d / %d]  %s', k, numel(hits), hits(k).name);
                if mod(k, 25) == 1, drawnow limitrate; end

                fpath = fullfile(hits(k).folder, hits(k).name);
                if ~isVersionedName(hits(k).name)
                    nSkippedVer = nSkippedVer + 1;
                    continue;
                end
                if ~isempty(D) && any(strcmpi(D(:,1), fpath))
                    nSkippedDup = nSkippedDup + 1;
                    continue;
                end
                [ani, cond, ph] = parseFilename(hits(k).name);
                D(end+1,:) = {fpath, ani, cond, ph}; %#ok<AGROW>
                nAdded = nAdded + 1;
            end

            tbl.Data = D;
            status.Text = sprintf( ...
                '%d file(s) queued. Added %d, skipped %d unversioned, %d duplicate(s).', ...
                size(D,1), nAdded, nSkippedVer, nSkippedDup);
        catch ME
            safeClose(pd);
            uialert(fig, sprintf('Folder scan failed: %s', ME.message), ...
                'Error', 'Icon','error');
            return;
        end

        safeClose(pd);
        drawnow;
    end

    function onRem()
        sel = getappdata(tbl,'sel');
        if isempty(sel), return; end
        D = tbl.Data;
        D(unique(sel(:,1)),:) = [];
        tbl.Data = D;
        status.Text = sprintf('%d file(s) queued.', size(D,1));
    end
    function onClear()
        tbl.Data = cell(0,4);
        status.Text = '0 file(s) queued.';
    end
    function onGo()
        fprintf('[buildFileQueue] Continue clicked. %d row(s) in table.\n', ...
            size(tbl.Data,1));
        uiresume(fig);
    end
    function onCancel()
        fprintf('[buildFileQueue] Cancel/close requested.\n');
        cancelled = true;
        uiresume(fig);
    end
end


function D = stateToTable(state)
    n = numel(state.files);
    if n == 0, D = cell(0,4); return; end
    D = cell(n,4);
    for i = 1:n
        D{i,1} = state.files{i};
        D{i,2} = getOrDefault(state.animal,    i, '');
        D{i,3} = getOrDefault(state.condition, i, '');
        D{i,4} = getOrDefault(state.phase,     i, '');
    end
end

function v = getOrDefault(c, i, default)
    if iscell(c) && i <= numel(c) && ~isempty(c{i})
        v = c{i};
    else
        v = default;
    end
end


function safeClose(pd)
% Close a uiprogressdlg if it's still alive (no-op otherwise).
    try
        if ~isempty(pd) && isvalid(pd)
            close(pd);
        end
    catch
        % swallow — the dialog might already be torn down
    end
end


function tf = isVersionedName(name)
% True if the file's basename contains a semver-style version tag of the
% form _v<digits>(.<digits>)+_  (e.g. _v0.1.0_, _v1.2_, _v0.10.3_).
% Edit the regex here to relax or tighten the rule.
    [~, base, ~] = fileparts(name);
    tf = ~isempty(regexp(base, '_v\d+(\.\d+)+_', 'once'));
end


function [ani, cond, ph] = parseFilename(filename)
% Heuristic auto-parse used when adding a file.
    [~, base, ~] = fileparts(filename);
    tokens = strsplit(base, '_');

    % animal: first letter of 2nd underscore-token
    if numel(tokens) >= 2 && ~isempty(tokens{2})
        ani = lower(tokens{2}(1));
    elseif ~isempty(tokens) && ~isempty(tokens{1})
        ani = lower(tokens{1}(1));
    else
        ani = '?';
    end

    % condition: first token that matches the M/E stim pattern
    cond = '';
    condRe = '^[ME]\d+([ME]\d+)?$';
    for k = 1:numel(tokens)
        t = upper(tokens{k});
        if ~isempty(regexp(t, condRe, 'once'))
            cond = t;
            break;
        end
    end
    if isempty(cond) && ~isempty(tokens)
        cond = tokens{1};
    end

    % phase
    if contains(base, 'recovery', 'IgnoreCase', true)
        ph = 'recovery';
    elseif contains(base, 'baseline','IgnoreCase', true) || ...
           ~isempty(regexp(base, '(^|_)bl(_|$)', 'once'))
        ph = 'baseline';
    else
        ph = 'baseline';
    end
end
