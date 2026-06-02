function Q = bulk_build_queue(cacheFile)
% BULK_BUILD_QUEUE  Interactive queue of (neural, R-peak, animal, condition, phase).
%
%   Q = bulk_build_queue(cacheFile)
%
% Neural files are *_blankmotion.mat (yOut, cols 1-2 = RVN/LVN). The matching
% R-peak file is auto-guessed as the sibling <neuralbase>_HRBR.mat (heartlocs);
% edit it inline if the guess is wrong. animal/condition/phase are auto-parsed
% from the neural basename (editable). Persists to cacheFile on Continue.
%
% Returns Q with cell-column fields: neural, rpeak, animal, condition, phase
% (Q.neural = {} if cancelled).

    if nargin < 1; cacheFile = ''; end
    Q = struct('neural',{{}},'rpeak',{{}},'animal',{{}},'condition',{{}},'phase',{{}});
    if ~isempty(cacheFile) && exist(cacheFile,'file')
        try; S = load(cacheFile); if isfield(S,'Q'); Q = S.Q; end; catch; end
    end
    Q = ui(Q);
    if ~isempty(Q.neural) && ~isempty(cacheFile)
        try; save(cacheFile,'Q'); catch; end
    end
end

% ======================================================================
function Q = ui(Q)
    fig = uifigure('Name','Bulk queue — neural + R-peak files','Position',[150 150 1180 560]);
    g = uigridlayout(fig,[3 1]); g.RowHeight = {'fit','1x','fit'};

    bp = uipanel(g,'BorderType','none'); bp.Layout.Row=1;
    bl = uigridlayout(bp,[1 6]); bl.ColumnWidth={120,140,140,120,'1x',110}; bl.Padding=[0 4 0 4];
    bAdd = uibutton(bl,'Text','Add files...');
    bDir = uibutton(bl,'Text','Add folder...');
    bRem = uibutton(bl,'Text','Remove selected');
    bClr = uibutton(bl,'Text','Clear');
    uilabel(bl,'Text','');
    bGo  = uibutton(bl,'Text','Continue','BackgroundColor',[0.7 0.9 0.7]);

    tbl = uitable(g, 'ColumnName',{'Neural file','R-peak file','Animal','Condition','Phase'}, ...
        'ColumnEditable',[false true true true true], ...
        'ColumnWidth',{'auto',260,70,110,90}, 'Data', to_table(Q), ...
        'CellSelectionCallback', @(s,e) setappdata(s,'sel',e.Indices));
    tbl.Layout.Row = 2;
    st = uilabel(g,'Text',sprintf('%d file(s) queued.', size(tbl.Data,1))); st.Layout.Row=3;

    bAdd.ButtonPushedFcn = @(~,~) onAdd();
    bDir.ButtonPushedFcn = @(~,~) onDir();
    bRem.ButtonPushedFcn = @(~,~) onRem();
    bClr.ButtonPushedFcn = @(~,~) onClr();
    bGo.ButtonPushedFcn  = @(~,~) uiresume(fig);
    fig.CloseRequestFcn  = @(~,~) onCancel();

    cancelled = false; uiwait(fig);
    if cancelled || isempty(tbl.Data)
        Q = struct('neural',{{}},'rpeak',{{}},'animal',{{}},'condition',{{}},'phase',{{}});
    else
        D = tbl.Data;
        Q.neural=D(:,1); Q.rpeak=D(:,2); Q.animal=D(:,3); Q.condition=D(:,4); Q.phase=D(:,5);
    end
    if isvalid(fig); delete(fig); end

    function onAdd()
        [names,folder] = uigetfile({'*blankmotion.mat;*.mat','Neural .mat'}, ...
            'Add neural files','MultiSelect','on');
        if isequal(names,0); return; end
        if ischar(names); names={names}; end
        D = tbl.Data;
        for i=1:numel(names)
            np = fullfile(folder,names{i});
            D(end+1,:) = make_row(np); %#ok<AGROW>
        end
        tbl.Data = D; st.Text = sprintf('%d file(s) queued.', size(D,1));
    end
    function onDir()
        root = uigetdir(pwd,'Pick a folder to scan recursively'); if isequal(root,0); return; end
        hits = dir(fullfile(root,'**','*blankmotion.mat'));
        hits = hits(~[hits.isdir] & ~contains({hits.name},'_HRBR'));
        D = tbl.Data;
        for i=1:numel(hits)
            np = fullfile(hits(i).folder,hits(i).name);
            if ~isempty(D) && any(strcmpi(D(:,1),np)); continue; end
            D(end+1,:) = make_row(np); %#ok<AGROW>
        end
        tbl.Data = D; st.Text = sprintf('%d file(s) queued.', size(D,1));
    end
    function onRem()
        sel = getappdata(tbl,'sel'); if isempty(sel); return; end
        D = tbl.Data; D(unique(sel(:,1)),:) = []; tbl.Data=D;
        st.Text = sprintf('%d file(s) queued.', size(D,1));
    end
    function onClr(); tbl.Data = cell(0,5); st.Text='0 file(s) queued.'; end
    function onCancel(); cancelled=true; uiresume(fig); end
end

% ======================================================================
function row = make_row(neuralPath)
    [folder,base,~] = fileparts(neuralPath);
    rpk = fullfile(folder,[base '_HRBR.mat']);
    if ~exist(rpk,'file')
        alt = dir(fullfile(folder,[base '*HRBR*.mat']));
        if ~isempty(alt); rpk = fullfile(alt(1).folder,alt(1).name); else; rpk = ''; end
    end
    [ani,cond,ph] = parse_name(base);
    row = {neuralPath, rpk, ani, cond, ph};
end

function [ani,cond,ph] = parse_name(base)
    tok = strsplit(base,'_');
    if numel(tok)>=2 && ~isempty(tok{2}); ani=lower(tok{2}(1));
    elseif ~isempty(tok); ani=lower(tok{1}(1)); else; ani='?'; end
    cond = '';
    for k=1:numel(tok)
        if ~isempty(regexp(upper(tok{k}),'^[ME]\d+([ME]\d+)?$','once')); cond=upper(tok{k}); break; end
    end
    if isempty(cond) && ~isempty(tok); cond=tok{1}; end
    if contains(base,'recovery','IgnoreCase',true) || contains(base,'stim_rec','IgnoreCase',true)
        ph='recovery';
    elseif contains(base,'baseline','IgnoreCase',true) || ~isempty(regexp(base,'(^|_)bl(_|$)','once'))
        ph='baseline';
    else; ph='baseline'; end
end

function D = to_table(Q)
    n = numel(Q.neural); D = cell(n,5);
    for i=1:n
        D{i,1}=Q.neural{i};
        D{i,2}=cell_or(Q.rpeak,i); D{i,3}=cell_or(Q.animal,i);
        D{i,4}=cell_or(Q.condition,i); D{i,5}=cell_or(Q.phase,i);
    end
end
function v = cell_or(c,i); if iscell(c)&&i<=numel(c)&&~isempty(c{i}); v=c{i}; else; v=''; end; end
