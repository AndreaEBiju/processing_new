function fig = me_heatmap_render(Z, mL, eL, titleStr, cbLabel)
% ME_HEATMAP_RENDER  Draw an M x E mean-%-change heatmap with the marginals
% placed OUTSIDE the axis labels: E-alone as the top row (row 1), M-alone as the
% left column (col 1); interior is the M x E combinations.
%
%   Z   : (nM+1) x (nE+1) matrix. Z(1,1)=NaN corner; Z(1,2:end)=E-alone;
%         Z(2:end,1)=M-alone; Z(2:end,2:end)=combined. NaN cells are blanked.
%   mL  : mechanical levels (length nM); eL : electrical levels (length nE)

    nM = numel(mL); nE = numel(eL);
    fig = figure('Color','w','Name',titleStr,'Position',[200 200 740 660]);
    ax  = axes(fig);
    mx = max(abs(Z(:)),[],'omitnan'); if isempty(mx)||~isfinite(mx)||mx==0; mx = 1; end
    Pm = Z; Pm(isnan(Pm)) = 0;
    h = imagesc(ax, Pm, [-mx mx]); set(h,'AlphaData',~isnan(Z));
    colormap(ax, rdbu(64)); cb = colorbar(ax);
    cb.Label.String = cbLabel; cb.Label.Interpreter = 'none';

    xl = [{'M alone'}, arrayfun(@(x)sprintf('E%d',x), eL, 'UniformOutput',false)];
    yl = [{'E alone'}, arrayfun(@(x)sprintf('M%d',x), mL, 'UniformOutput',false)];
    set(ax, 'XTick',1:nE+1, 'XTickLabel',xl, 'YTick',1:nM+1, 'YTickLabel',yl, ...
        'TickLabelInterpreter','none', 'XAxisLocation','top', 'Color',[0.92 0.92 0.92]);
    axis(ax,'equal'); axis(ax,'tight');
    title(ax, titleStr, 'Interpreter','none');

    for r = 1:nM+1
        for c = 1:nE+1
            if isnan(Z(r,c)); continue; end
            tc = [0 0 0]; if abs(Z(r,c))/mx > 0.55; tc = [1 1 1]; end
            text(ax, c, r, sprintf('%+.0f%%', Z(r,c)), 'HorizontalAlignment','center', ...
                'Color',tc, 'FontSize',13);
        end
    end
    hold(ax,'on');
    plot(ax, [1.5 1.5], [0.5 nM+1.5], 'k-', 'LineWidth',1.2);
    plot(ax, [0.5 nE+1.5], [1.5 1.5], 'k-', 'LineWidth',1.2);
end

function cm = rdbu(n)
    if nargin<1; n=64; end; h=floor(n/2);
    cm = [ [linspace(0.05,1,h).' linspace(0.27,1,h).' linspace(0.55,1,h).']; ...
           [linspace(1,0.65,n-h).' linspace(1,0.10,n-h).' linspace(1,0.13,n-h).'] ];
end
