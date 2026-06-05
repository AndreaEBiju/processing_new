function restoreFn = whiten_figure(fig)
% WHITEN_FIGURE  Enforce the shared export style on FIG, then return a handle
% restoreFn() that puts the on-screen appearance back afterwards.
%
% Enforced (uniform across every saved plot):
%   * white figure + axes background (InvertHardcopy off)
%   * black axis rulers, axis title and x/y/z labels, colorbar ruler + label
%   * font size 20 for ALL text: tick labels, axis labels, titles, colorbar,
%     and every manually-placed text() object (e.g. heatmap cell/edge labels)
%
% Left untouched: data colors (lines, patches, scatter, box fills, heatmap
% cell fills) and the COLOR of manual text() objects -- so white-on-saturated
% heatmap cell numbers stay white and the dark-cell-on-light style is kept.

    FS = 20;

    figColor = get(fig,'Color'); invert = get(fig,'InvertHardcopy');
    set(fig,'Color','w','InvertHardcopy','off');

    % ---- axes: background, ruler/label colors, font size ----
    axs = findall(fig,'Type','axes');
    axState = cell(numel(axs),1);
    for i = 1:numel(axs)
        a = axs(i);
        axState{i} = {get(a,'Color'), get(a,'XColor'), get(a,'YColor'), get(a,'ZColor'), ...
                      a.Title.Color, a.XLabel.Color, a.YLabel.Color, ...
                      get(a,'FontSize'), get(a,'LabelFontSizeMultiplier'), get(a,'TitleFontSizeMultiplier'), ...
                      a.Title.FontSize, a.XLabel.FontSize, a.YLabel.FontSize};
        set(a,'Color','w','XColor','k','YColor','k','ZColor','k', ...
              'FontSize',FS,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1);
        a.Title.Color = 'k'; a.XLabel.Color = 'k'; a.YLabel.Color = 'k';
        a.Title.FontSize = FS; a.XLabel.FontSize = FS; a.YLabel.FontSize = FS;
    end

    % ---- every text object (incl. manual heatmap labels): font size only ----
    txts = findall(fig,'Type','text');
    txtState = cell(numel(txts),1);
    for i = 1:numel(txts)
        t = txts(i);
        txtState{i} = get(t,'FontSize');
        try, set(t,'FontSize',FS); catch, end %#ok<CTCH>
    end

    % ---- colorbars: ruler + label color, font size ----
    cbs = findall(fig,'Type','colorbar');
    cbState = cell(numel(cbs),1);
    for i = 1:numel(cbs)
        c = cbs(i);
        cbState{i} = {get(c,'Color'), c.Label.Color, get(c,'FontSize'), c.Label.FontSize};
        set(c,'Color','k','FontSize',FS); c.Label.Color = 'k'; c.Label.FontSize = FS;
    end

    restoreFn = @doRestore;

    function doRestore()
        if isgraphics(fig)
            try, set(fig,'Color',figColor,'InvertHardcopy',invert); catch, end %#ok<CTCH>
        end
        for ii = 1:numel(axs)
            a = axs(ii); if ~isgraphics(a), continue; end; s = axState{ii};
            try
                set(a,'Color',s{1},'XColor',s{2},'YColor',s{3},'ZColor',s{4}, ...
                      'FontSize',s{8},'LabelFontSizeMultiplier',s{9},'TitleFontSizeMultiplier',s{10});
                a.Title.Color = s{5}; a.XLabel.Color = s{6}; a.YLabel.Color = s{7};
                a.Title.FontSize = s{11}; a.XLabel.FontSize = s{12}; a.YLabel.FontSize = s{13};
            catch, end %#ok<CTCH>
        end
        for ii = 1:numel(txts)
            t = txts(ii); if ~isgraphics(t), continue; end
            try, set(t,'FontSize',txtState{ii}); catch, end %#ok<CTCH>
        end
        for ii = 1:numel(cbs)
            c = cbs(ii); if ~isgraphics(c), continue; end; s = cbState{ii};
            try, set(c,'Color',s{1},'FontSize',s{3}); c.Label.Color = s{2}; c.Label.FontSize = s{4}; catch, end %#ok<CTCH>
        end
    end
end
