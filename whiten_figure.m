function restoreFn = whiten_figure(fig)
% WHITEN_FIGURE  Apply the one uniform plot style to FIG (permanent):
%   * white figure + axes background (InvertHardcopy off)
%   * black axis rulers, axis title and x/y/z labels
%   * black tiledlayout super-title / labels
%   * white legend box with black text and edge
%   * black colorbar ruler + label
%   * font size 20 on ALL text (ticks, labels, titles, legend, colorbar,
%     and every manual text() object, e.g. heatmap cell/edge labels)
%
% Data colors (lines, patches, scatter, box fills, heatmap cell fills) are left
% alone, and manual text() keeps its color ONLY when it is deliberately white
% (white-on-saturated heatmap numbers). Everything else becomes black.
%
% Returns a no-op restore handle for backward compatibility (style is permanent).

    FS = 20;
    if nargin < 1 || isempty(fig) || ~isgraphics(fig); restoreFn = @()[]; return; end

    set(fig,'Color','w','InvertHardcopy','off');

    % ---- axes ----
    axs = findall(fig,'Type','axes');
    for i = 1:numel(axs)
        a = axs(i);
        try
            set(a,'Color','w','XColor','k','YColor','k','ZColor','k', ...
                  'FontSize',FS,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1);
            a.Title.Color = 'k'; a.XLabel.Color = 'k'; a.YLabel.Color = 'k';
            a.Title.FontSize = FS; a.XLabel.FontSize = FS; a.YLabel.FontSize = FS;
        catch, end %#ok<CTCH>
    end

    % ---- tiledlayout super-title / labels ----
    tls = findall(fig,'Type','tiledlayout');
    for i = 1:numel(tls)
        for pr = {'Title','Subtitle','XLabel','YLabel'}
            try
                h = tls(i).(pr{1});
                if ~isempty(h) && isgraphics(h); h.Color = 'k'; h.FontSize = FS; end
            catch, end %#ok<CTCH>
        end
    end

    % ---- legends ----
    lg = findall(fig,'Type','legend');
    for i = 1:numel(lg)
        try, set(lg(i),'Color','w','TextColor','k','EdgeColor',[0.15 0.15 0.15],'FontSize',FS); catch, end %#ok<CTCH>
    end

    % ---- every text object: font 20, black unless deliberately white ----
    txts = findall(fig,'Type','text');
    for i = 1:numel(txts)
        t = txts(i);
        try
            set(t,'FontSize',FS);
            c = get(t,'Color');
            if ~all(c >= 0.95); set(t,'Color','k'); end   % keep near-white text white
        catch, end %#ok<CTCH>
    end

    % ---- colorbars ----
    cbs = findall(fig,'Type','colorbar');
    for i = 1:numel(cbs)
        try, set(cbs(i),'Color','k','FontSize',FS); cbs(i).Label.Color = 'k'; cbs(i).Label.FontSize = FS; catch, end %#ok<CTCH>
    end

    restoreFn = @()[];   % style is permanent; no-op
end
