function [groups, ok] = defineGroupsUI(conditions, priorGroups)
%DEFINEGROUPSUI  Interactive UI to define condition groupings.
%
%   [groups, ok] = defineGroupsUI(conditions, priorGroups)
%     conditions  : cellstr of detected conditions (shown at the top)
%     priorGroups : cell of cellstr — previous groupings (each cell is a
%                   group, holding the condition names in it). Pre-fills
%                   the table; pass {} on first run.
%
%   Returns
%     groups : cell of cellstr (each cell is one group's conditions)
%     ok     : true if the user clicked Continue with at least one
%              non-empty group, false on cancel / empty input.

    fprintf('[defineGroupsUI] Opening. Detected conditions: %s\n', ...
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
    if isempty(initData), initData = {'Group 1', ''}; end
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
        fprintf('[defineGroupsUI] Continue clicked. %d row(s).\n', ...
            size(tbl.Data,1));
        uiresume(fig);
    end
    function onCancel()
        fprintf('[defineGroupsUI] Cancel/close.\n');
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
