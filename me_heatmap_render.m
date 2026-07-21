function fig = me_heatmap_render(Z, mL, eL, titleStr, cbLabel, sigMask, valueFmt)
% ME_HEATMAP_RENDER  M x E mean-%-change heatmap drawn as separated tiles, with
% the marginals OUTSIDE the interior 3x3 grid (matching the mockup):
%   * interior 3x3 (M rows x E cols)
%   * E-alone row sits ABOVE the E column labels (top), separated by a gap
%   * M-alone column sits to the LEFT of the M row labels, separated by a gap
%
%   Z : (nM+1) x (nE+1).  Z(1,1)=corner (NaN); Z(1,2:end)=E-alone;
%       Z(2:end,1)=M-alone; Z(2:end,2:end)=combined. NaN cells are blanked.
%   sigMask (optional): logical, same shape as Z. Significant cells get a bold
%       black border and a trailing '*' on the value.
%   valueFmt (optional): sprintf format for the per-cell value label, default
%       '%+.0f%%' (existing % change from baseline convention -- unchanged for
%       every pre-existing caller). Pass e.g. '%.2f' for non-percent
%       quantities (a model's predicted absolute rate_hat, say) so cells don't
%       get a misleading trailing '%'.

    mL = mL(:)'; eL = eL(:)';
    nM = numel(mL); nE = numel(eL);
    if nargin < 6 || isempty(sigMask); sigMask = false(size(Z)); end
    if nargin < 7 || isempty(valueFmt); valueFmt = '%+.0f%%'; end
    Ealone = Z(1, 2:end);       % 1 x nE
    Malone = Z(2:end, 1);       % nM x 1
    Comb   = Z(2:end, 2:end);   % nM x nE
    SgE = sigMask(1, 2:end);  SgM = sigMask(2:end, 1);  SgC = sigMask(2:end, 2:end);

    fig = figure('Color','w','Name',titleStr,'Position',[200 200 780 660]);
    ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    set(ax,'YDir','reverse');

    allv = [Ealone(:); Malone(:); Comb(:)];
    mx = max(abs(allv),[],'omitnan'); if isempty(mx)||~isfinite(mx)||mx==0; mx = 1; end
    cmap = rdbu(256);

    % interior: E columns at x = 1..nE, M rows at y = 1..nM
    for i = 1:nM
        for j = 1:nE; draw_cell(ax, j, i, Comb(i,j), mx, cmap, SgC(i,j), valueFmt); end
    end
    % M-alone column at x = -1 (gap at x = 0 = the M row-label column)
    for i = 1:nM; draw_cell(ax, -1, i, Malone(i), mx, cmap, SgM(i), valueFmt); end
    % E-alone row at y = -1 (gap at y = 0 = the E column-label row)
    for j = 1:nE; draw_cell(ax, j, -1, Ealone(j), mx, cmap, SgE(j), valueFmt); end

    % axis labels live in the gap row/column (y=0 for E, x=0 for M)
    for j = 1:nE; text(ax, j, 0, sprintf('E%d',eL(j)), 'HorizontalAlignment','center','FontSize',20,'Color',[0 0 0]); end
    for i = 1:nM; text(ax, 0, i, sprintf('M%d',mL(i)), 'HorizontalAlignment','center','FontSize',20,'Color',[0 0 0]); end
    text(ax, -1, 0, 'M alone', 'HorizontalAlignment','center','FontSize',20,'FontAngle','italic','Color',[0 0 0]);
    text(ax,  0,-1, 'E alone', 'HorizontalAlignment','center','FontSize',20,'FontAngle','italic','Color',[0 0 0]);

    xlim(ax, [-1.7, nE+0.7]); ylim(ax, [-1.7, nM+0.7]);
    title(ax, titleStr, 'Interpreter','none', 'FontSize',20, 'Color',[0 0 0]);

    colormap(ax, cmap); caxis(ax, [-mx mx]);
    cb = colorbar(ax); cb.Label.String = cbLabel; cb.Label.Interpreter = 'none';
    whiten_figure(fig);   % uniform style: white bg, black text, font 20 (incl. colorbar)
end

% ----------------------------------------------------------------------
function draw_cell(ax, xc, yc, v, mx, cmap, isSig, valueFmt)
    if nargin < 7 || isempty(isSig); isSig = false; end
    if nargin < 8 || isempty(valueFmt); valueFmt = '%+.0f%%'; end
    w = 0.46;
    if isnan(v)
        col = [0.92 0.92 0.92];
    else
        t = (v + mx)/(2*mx); t = min(max(t,0),1);
        col = cmap(max(1, round(t*(size(cmap,1)-1))+1), :);
    end
    if isSig && ~isnan(v); ec = [0 0 0]; lw = 3; else; ec = [1 1 1]; lw = 1; end
    patch(ax, xc+[-w -w w w], yc+[-w w w -w], col, 'EdgeColor',ec, 'LineWidth',lw);
    if ~isnan(v)
        tc = [0 0 0]; if abs(v)/mx > 0.55; tc = [1 1 1]; end
        str = sprintf(valueFmt, v); if isSig; str = [str '*']; end
        text(ax, xc, yc, str, 'HorizontalAlignment','center', 'Color',tc, 'FontSize',20);
    end
end

function cm = rdbu(n)
    if nargin<1; n=64; end; h=floor(n/2);
    cm = [ [linspace(0.05,1,h).' linspace(0.27,1,h).' linspace(0.55,1,h).']; ...
           [linspace(1,0.65,n-h).' linspace(1,0.10,n-h).' linspace(1,0.13,n-h).'] ];
end
