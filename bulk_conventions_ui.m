function conv = bulk_conventions_ui(cacheFile)
% BULK_CONVENTIONS_UI  Ask once for the file-naming conventions applied to every
% file in the bulk run. Persisted to cacheFile so it pre-fills next time.
%
% Returns conv with fields:
%   sufBaseNeural, sufRecNeural  neural-file suffixes (before '.mat')
%   sufBaseHRBR,   sufRecHRBR    matching heartbeat-file suffixes
%   rvnCol, lvnCol               column index (1-5) of RVN and LVN in yOut
%   heartVar                     variable name holding R-peak sample indices
% conv = [] if cancelled.

    if nargin < 1; cacheFile = ''; end
    d = default_conv();
    cached = bulk_cache_get(cacheFile, 'conv');
    if isstruct(cached) && ~isempty(fieldnames(cached)); d = cached; end

    fig = uifigure('Name','Bulk naming conventions','Position',[250 250 620 420]);
    g = uigridlayout(fig,[9 2]);
    g.RowHeight = repmat({'fit'},1,9); g.ColumnWidth = {230,'1x'};

    lbl(g,'Baseline neural suffix:');   eBN = ef(g, d.sufBaseNeural);
    lbl(g,'Recovery neural suffix:');   eRN = ef(g, d.sufRecNeural);
    lbl(g,'Baseline HRBR suffix:');     eBH = ef(g, d.sufBaseHRBR);
    lbl(g,'Recovery HRBR suffix:');     eRH = ef(g, d.sufRecHRBR);
    lbl(g,'RVN channel (col 1-5):');    eRV = efn(g, d.rvnCol);
    lbl(g,'LVN channel (col 1-5):');    eLV = efn(g, d.lvnCol);
    lbl(g,'Heartbeat variable name:');  eHV = ef(g, d.heartVar);

    note = uilabel(g,'Text','(suffixes exclude the .mat extension)','FontAngle','italic');
    note.Layout.Column = [1 2];
    bGo = uibutton(g,'Text','Continue','BackgroundColor',[0.7 0.9 0.7]);
    bGo.Layout.Column = 2;

    conv = []; bGo.ButtonPushedFcn = @(~,~) onGo(); fig.CloseRequestFcn = @(~,~) uiresume(fig);
    uiwait(fig);
    if isvalid(fig); delete(fig); end
    if ~isempty(conv); bulk_cache_set(cacheFile, 'conv', conv); end

    function onGo()
        conv = struct('sufBaseNeural',strtrim(eBN.Value), 'sufRecNeural',strtrim(eRN.Value), ...
                      'sufBaseHRBR',strtrim(eBH.Value),   'sufRecHRBR',strtrim(eRH.Value), ...
                      'rvnCol',round(eRV.Value), 'lvnCol',round(eLV.Value), ...
                      'heartVar',strtrim(eHV.Value));
        uiresume(fig);
    end
end

function d = default_conv()
    d = struct('sufBaseNeural','_v0.2.2_blankmotion', ...
               'sufRecNeural','_v0.2.2_recovery_blankmotion', ...
               'sufBaseHRBR','_v0.2.2_blankmotion_HRBR', ...
               'sufRecHRBR','_v0.2.2_recovery_HRBR', ...
               'rvnCol',1, 'lvnCol',2, 'heartVar','heartlocs');
end

function lbl(g,t); u = uilabel(g,'Text',t); u.Layout.Column = 1; end
function e = ef(g,v);  e = uieditfield(g,'text','Value',v);  e.Layout.Column = 2; end
function e = efn(g,v); e = uieditfield(g,'numeric','Value',v,'Limits',[1 5],'RoundFractionalValues',true); e.Layout.Column = 2; end
