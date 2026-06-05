function fig = bulk_plot_synergy(normRows, metric, chLabel)
% BULK_PLOT_SYNERGY  3x3 M x E synergy heatmap for one metric, one channel.
%
%   fig = bulk_plot_synergy(normRows, metric, chLabel)
%
% normRows come from bulk_compile and already carry the per-(animal,condition)
% baseline-normalized scalar effect in .scalar.(metric). For each cell:
%   norm(cond)     = mean over animals of scalar(metric)
%   synergy(Mi,Ej) = norm(MiEj) - [ norm(Mi) + norm(Ej) ]
% red = super-additive, blue = sub-additive, white ~ additive.
%
%   metric : 'rate' | 'excess' | 'vpp' | 'fwhm' | 'cv2' | 'fano'

    conds = unique({normRows.condition});
    [Mlev, Elev] = parse_ME(conds);
    eLevels = sort(unique(Elev(Elev>0 & Mlev==0 & ~isnan(Elev))));
    mLevels = sort(unique(Mlev(Mlev>0 & Elev==0 & ~isnan(Mlev))));
    if isempty(eLevels) || isempty(mLevels)
        warning('bulk_plot_synergy:levels','Need both M-alone and E-alone conditions.'); fig = []; return;
    end
    nM = numel(mLevels); nE = numel(eLevels);

    synergy = nan(nM, nE);
    for i = 1:nM
        for j = 1:nE
            nME = norm_effect(normRows, sprintf('M%dE%d', mLevels(i), eLevels(j)), chLabel, metric);
            nMa = norm_effect(normRows, sprintf('M%d', mLevels(i)),  chLabel, metric);
            nEa = norm_effect(normRows, sprintf('E%d', eLevels(j)),  chLabel, metric);
            if all(isfinite([nME nMa nEa])); synergy(i,j) = nME - (nMa + nEa); end
        end
    end

    fig = figure('Color','w','Name',sprintf('synergy %s — %s', metric, chLabel), ...
        'Position',[200 200 720 620]);
    ax = axes(fig);
    mx = max(abs(synergy(:)),[],'omitnan'); if isempty(mx)||~isfinite(mx)||mx==0; mx = 1; end
    Pm = synergy; Pm(isnan(Pm)) = 0;
    h = imagesc(ax, Pm, [-mx mx]); set(h,'AlphaData',~isnan(synergy));
    colormap(ax, rdbu(64)); cb = colorbar(ax);
    cb.Label.String = 'synergy = norm(M+E) - [norm(M)+norm(E)]'; cb.Label.Interpreter='none';
    set(ax,'XTick',1:nE,'XTickLabel',arrayfun(@(x)sprintf('E%d',x),eLevels,'UniformOutput',false), ...
           'YTick',1:nM,'YTickLabel',arrayfun(@(x)sprintf('M%d',x),mLevels,'UniformOutput',false), ...
           'TickLabelInterpreter','none','YDir','reverse','XAxisLocation','top','Color',[0.92 0.92 0.92]);
    axis(ax,'equal'); axis(ax,'tight');
    xlabel(ax,'Electrical level (E)'); ylabel(ax,'Mechanical level (M)');
    title(ax, sprintf('%s — synergy  |  %s  (norm = (rec-base)/base)', upper(metric), chLabel), 'Interpreter','none');
    for i = 1:nM
        for j = 1:nE
            if isnan(synergy(i,j))
                text(ax,j,i,'NaN','HorizontalAlignment','center','Color',[0.45 0.45 0.45]);
            else
                tc = [0 0 0]; if abs(synergy(i,j))/mx > 0.55; tc = [1 1 1]; end
                text(ax,j,i,sprintf('%+.3f',synergy(i,j)),'HorizontalAlignment','center','Color',tc,'FontSize',13);
            end
        end
    end
    whiten_figure(fig);   % uniform style: white bg, black text, font 20
end

% ======================================================================
function v = norm_effect(normRows, cond, chLabel, metric)
    mm = lower(metric); acc = [];
    for r = 1:numel(normRows)
        R = normRows(r);
        if strcmpi(R.condition,cond) && strcmpi(R.label,chLabel) && isfield(R.scalar,mm)
            acc(end+1) = R.scalar.(mm); %#ok<AGROW>
        end
    end
    acc = acc(isfinite(acc));
    if isempty(acc); v = NaN; else; v = mean(acc); end
end

function [Mlev, Elev] = parse_ME(conds)
    n = numel(conds); Mlev = nan(n,1); Elev = nan(n,1);
    for k = 1:n
        c = upper(strtrim(conds{k}));
        t = regexp(c,'^M(\d+)E(\d+)$','tokens','once');
        if ~isempty(t); Mlev(k)=str2double(t{1}); Elev(k)=str2double(t{2}); continue; end
        t = regexp(c,'^E(\d+)M(\d+)$','tokens','once');
        if ~isempty(t); Elev(k)=str2double(t{1}); Mlev(k)=str2double(t{2}); continue; end
        t = regexp(c,'^M(\d+)$','tokens','once');
        if ~isempty(t); Mlev(k)=str2double(t{1}); Elev(k)=0; continue; end
        t = regexp(c,'^E(\d+)$','tokens','once');
        if ~isempty(t); Elev(k)=str2double(t{1}); Mlev(k)=0; end
    end
end

function cm = rdbu(n)
    if nargin<1; n=64; end; h=floor(n/2);
    cm = [ [linspace(0.05,1,h).' linspace(0.27,1,h).' linspace(0.55,1,h).']; ...
           [linspace(1,0.65,n-h).' linspace(1,0.10,n-h).' linspace(1,0.13,n-h).'] ];
end
