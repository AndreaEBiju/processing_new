function files = bulk_queue_ui(conv, cacheFile)
% BULK_QUEUE_UI  Table-based file queue (like your buildFileQueue): "Add
% folder..." recursively scans (via bulk_scan_files + the conventions in conv)
% and drops every matched file into one editable table. Review, fix the
% phase/condition/animal/HRBR inline, remove unwanted rows, then Continue.
% Persists the queue to cacheFile.
%
% Returns a struct array (neural, hrbr, stem, animal, condition, phase);
% empty if cancelled.

    if nargin < 2; cacheFile = ''; end
    init = bulk_cache_get(cacheFile, 'queue');
    if isempty(init) || ~iscell(init); init = cell(0,5); end
    files = struct('neural',{},'hrbr',{},'stem',{},'animal',{},'condition',{},'phase',{});

    fig = uifigure('Name','Bulk file queue','Position',[150 150 1200 560]);
    g = uigridlayout(fig,[3 1]); g.RowHeight = {'fit','1x','fit'};
    bp = uipanel(g,'BorderType','none'); bp.Layout.Row = 1;
    bl = uigridlayout(bp,[1 7]); bl.ColumnWidth = {120,140,120,140,110,'1x',110}; bl.Padding = [0 4 0 4];
    bDir = uibutton(bl,'Text','Add folder...');
    bImp = uibutton(bl,'Text','Import queue...');
    bFil = uibutton(bl,'Text','Add files...');
    bRem = uibutton(bl,'Text','Remove selected');
    bClr = uibutton(bl,'Text','Clear');
    uilabel(bl,'Text','');
    bGo  = uibutton(bl,'Text','Continue','BackgroundColor',[0.7 0.9 0.7]);

    tbl = uitable(g, 'ColumnName',{'Neural file','HRBR file','Animal','Condition','Phase'}, ...
        'ColumnEditable',[false true true true true], 'ColumnWidth',{'auto',260,70,110,90}, ...
        'Data',init, 'CellSelectionCallback',@(s,e) setappdata(s,'sel',e.Indices));
    tbl.Layout.Row = 2;
    stt = uilabel(g,'Text',sprintf('%d file(s).', size(init,1))); stt.Layout.Row = 3;

    bDir.ButtonPushedFcn = @(~,~) onDir();
    bImp.ButtonPushedFcn = @(~,~) onImport();
    bFil.ButtonPushedFcn = @(~,~) onFiles();
    bRem.ButtonPushedFcn = @(~,~) onRem();
    bClr.ButtonPushedFcn = @(~,~) onClr();
    bGo.ButtonPushedFcn  = @(~,~) uiresume(fig);
    fig.CloseRequestFcn  = @(~,~) onCancel();

    cancelled = false; uiwait(fig);
    if ~cancelled && ~isempty(tbl.Data)
        D = tbl.Data;
        for i = 1:size(D,1)
            [~, bn] = fileparts(D{i,1});
            files(end+1) = struct('neural',D{i,1},'hrbr',D{i,2},'stem',bn, ...
                'animal',D{i,3},'condition',D{i,4},'phase',D{i,5}); %#ok<AGROW>
        end
        bulk_cache_set(cacheFile, 'queue', tbl.Data);
    end
    if isvalid(fig); delete(fig); end

    % ---- nested ----
    function appendScan(sf)
        D = tbl.Data; nAdd = 0;
        for i = 1:numel(sf)
            if ~isempty(D) && any(strcmpi(D(:,1), sf(i).neural)); continue; end
            D(end+1,:) = {sf(i).neural, sf(i).hrbr, sf(i).animal, sf(i).condition, sf(i).phase}; %#ok<AGROW>
            nAdd = nAdd + 1;
        end
        tbl.Data = D; stt.Text = sprintf('%d file(s) (added %d).', size(D,1), nAdd);
    end
    function onDir()
        root = uigetdir(pwd,'Pick a folder to scan recursively');
        if isequal(root,0); return; end
        appendScan(bulk_scan_files({root}, conv));
    end
    function onImport()
        % import an existing buildFileQueue cache (gemsplots_queue.mat): take its
        % files/animal/condition/phase, derive each HRBR sibling from conv.
        [fn,fp] = uigetfile({'*.mat','MAT-files'}, ...
            'Pick a buildFileQueue cache (e.g. gemsplots_queue.mat)');
        if isequal(fn,0); return; end
        try; S = load(fullfile(fp,fn)); catch ME; uialert(fig, ME.message, 'Import failed'); return; end
        if ~isfield(S,'files') || isempty(S.files)
            uialert(fig, 'No "files" list found in that .mat.', 'Import'); return;
        end
        nf = numel(S.files);
        sf = struct('neural',{},'hrbr',{},'stem',{},'animal',{},'condition',{},'phase',{});
        for i = 1:nf
            np = char(S.files{i});
            ph = field_at(S,'phase',i,'baseline');
            an = field_at(S,'animal',i,'');
            cd = field_at(S,'condition',i,'');
            [~, bn] = fileparts(np);
            sf(end+1) = struct('neural',np,'hrbr',resolve_hrbr(np,ph,conv),'stem',bn, ...
                'animal',an,'condition',cd,'phase',ph); %#ok<AGROW>
        end
        appendScan(sf);
    end
    function onFiles()
        [names,folder] = uigetfile({'*blankmotion.mat;*.mat','Neural .mat'}, ...
            'Add neural files','MultiSelect','on');
        if isequal(names,0); return; end
        if ischar(names); names = {names}; end
        sf = struct('neural',{},'hrbr',{},'stem',{},'animal',{},'condition',{},'phase',{});
        for i = 1:numel(names); sf(end+1) = guess_one(fullfile(folder,names{i}), conv); end %#ok<AGROW>
        appendScan(sf);
    end
    function onRem()
        sel = getappdata(tbl,'sel'); if isempty(sel); return; end
        D = tbl.Data; D(unique(sel(:,1)),:) = []; tbl.Data = D;
        stt.Text = sprintf('%d file(s).', size(D,1));
    end
    function onClr(); tbl.Data = cell(0,5); stt.Text = '0 file(s).'; end
    function onCancel(); cancelled = true; uiresume(fig); end
end

% ----------------------------------------------------------------------
function row = guess_one(np, conv)
    [~, base] = fileparts(np);
    phase = ''; stem = base;
    if endsWith(base, conv.sufRecNeural)
        phase = 'recovery'; stem = base(1:end-numel(conv.sufRecNeural));
    elseif endsWith(base, conv.sufBaseNeural) && ~contains(base,'stim_rec','IgnoreCase',true)
        phase = 'baseline'; stem = base(1:end-numel(conv.sufBaseNeural));
    end
    if ~isempty(phase); hrbr = resolve_hrbr(np, phase, conv); else; hrbr = ''; end
    [ani, cond] = parse_stem_local(stem);
    row = struct('neural',np,'hrbr',hrbr,'stem',stem,'animal',ani,'condition',cond,'phase',phase);
end

function v = field_at(S, name, i, default)
    v = default;
    if ~isfield(S,name); return; end
    c = S.(name);
    if iscell(c) && i<=numel(c) && ~isempty(c{i}); v = c{i};
    elseif ischar(c); v = c;
    elseif ~iscell(c) && i<=numel(c); v = c(i); end
end

function [ani, cond] = parse_stem_local(stem)
    tok = strsplit(stem,'_');
    if numel(tok)>=2 && ~isempty(tok{2}); ani = lower(tok{2}(1));
    elseif ~isempty(tok) && ~isempty(tok{1}); ani = lower(tok{1}(1)); else; ani = '?'; end
    cond = '';
    for k = 1:numel(tok)
        if ~isempty(regexp(upper(tok{k}),'^[ME]\d+([ME]\d+)?$','once')); cond = upper(tok{k}); break; end
    end
    if isempty(cond) && ~isempty(tok); cond = tok{1}; end
end
