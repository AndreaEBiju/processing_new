function pub = pubfig_setup(varargin)
%PUBFIG_SETUP  Publication-quality plotting defaults + save helpers (self-contained).

% ---------- Parse inputs ----------
p = inputParser;
addParameter(p,'FontName','Times New Roman',@(s)ischar(s)||isstring(s));
addParameter(p,'BaseFontSize',9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LineWidth',1.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MarkerSize',7,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'PlotsDir','plots',@(s)ischar(s)||isstring(s));
addParameter(p,'Renderer','painters',@(s)ischar(s)||isstring(s));
addParameter(p,'ColorOrder',[],@(x)isempty(x)||(isnumeric(x)&&size(x,2)==3));
addParameter(p,'Theme','light',@(s) any(strcmpi(string(s),["light","dark"])));
addParameter(p,'EnableLaTeX',true,@islogical);
parse(p,varargin{:});
opt = p.Results;

% Resolve theme -> background/foreground triplet + categorical palette
themeName = lower(char(opt.Theme));
[bgCol, fgCol, gridCol, minorGridCol, paletteCat] = theme_palette(themeName);
if isempty(opt.ColorOrder)
    opt.ColorOrder = paletteCat;
end

% ---------- Ensure plots directory exists ----------
plotsDir = char(opt.PlotsDir);
if ~exist(plotsDir,'dir'), mkdir(plotsDir); end
setappdata(0,'PUB_PLOTS_DIR',plotsDir);

% ---------- Figure sizes (cm) ----------
phi = (1+sqrt(5))/2;
W1 = 8.5;  H1 = W1/phi;   % single column
W2 = 18.0; H2 = 6.0;      % double column (shorter)
sizes = struct('single',[W1 H1], 'double',[W2 H2]);
setappdata(0,'PUBSIZE',sizes);   % <-- store as one struct (no dots in key)

% ---------- Root defaults ----------
set(0,'defaultFigureWindowStyle','normal');
try, set(0,'defaultFigureWindowState','maximized'); end %#ok
set(groot,'defaultFigureColor',bgCol);
set(groot,'defaultFigureRenderer',char(opt.Renderer));
set(groot,'defaultFigureInvertHardcopy','off');

% Fonts
set(groot,'defaultAxesFontName',char(opt.FontName));
set(groot,'defaultTextFontName',char(opt.FontName));
set(groot,'defaultLegendFontName',char(opt.FontName));
set(groot,'defaultAxesFontSize',opt.BaseFontSize);
set(groot,'defaultTextFontSize',opt.BaseFontSize);
set(groot,'defaultAxesLabelFontSizeMultiplier',1.0);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.0);
set(groot,'defaultAxesTitleFontWeight','normal');

% Interpreters
if opt.EnableLaTeX
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultColorbarTickLabelInterpreter','latex');
else
    set(groot,'defaultTextInterpreter','none');
    set(groot,'defaultAxesTickLabelInterpreter','none');
    set(groot,'defaultLegendInterpreter','none');
    set(groot,'defaultColorbarTickLabelInterpreter','none');
end

% Lines/markers
set(groot,'defaultLineLineWidth',opt.LineWidth);
set(groot,'defaultLineMarkerSize',opt.MarkerSize);
set(groot,'defaultStairLineWidth',opt.LineWidth);
set(groot,'defaultStemLineWidth',opt.LineWidth);
set(groot,'defaultErrorbarLineWidth',opt.LineWidth);
set(groot,'defaultConstantLineLineWidth',opt.LineWidth);

% Axes & colors (theme-aware)
set(groot,'defaultAxesLineWidth',1.0);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLength',[0.015 0.015]);
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesColor',bgCol);
set(groot,'defaultAxesXColor',fgCol);
set(groot,'defaultAxesYColor',fgCol);
set(groot,'defaultAxesZColor',fgCol);
set(groot,'defaultAxesGridColor',gridCol);
set(groot,'defaultAxesMinorGridColor',minorGridCol);
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot,'defaultAxesGridAlpha',0.25);
set(groot,'defaultAxesMinorGridAlpha',0.25);
set(groot,'defaultAxesXMinorGrid','off');
set(groot,'defaultAxesYMinorGrid','off');
set(groot,'defaultAxesColorOrder',opt.ColorOrder);
set(groot,'defaultAxesLooseInset',[0.02 0.02 0.02 0.02]);
set(groot,'defaultAxesSortMethod','childorder');
set(groot,'defaultTextColor',fgCol);

% Legends & tiling
set(groot,'defaultLegendBox','off');
set(groot,'defaultLegendLocation','best');
set(groot,'defaultLegendAutoUpdate','off');
set(groot,'defaultLegendColor',bgCol);
set(groot,'defaultLegendTextColor',fgCol);
set(groot,'defaultLegendEdgeColor',fgCol);
set(groot,'defaultTiledlayoutTileSpacing','compact');
set(groot,'defaultTiledlayoutPadding','compact');

% Stash theme/palette so downstream code (e.g. boxScatterPlot) can read
setappdata(0,'PUB_THEME',  themeName);
setappdata(0,'PUB_PALETTE',opt.ColorOrder);
setappdata(0,'PUB_FGCOL',  fgCol);
setappdata(0,'PUB_BGCOL',  bgCol);

% Colormap
set(groot,'defaultFigureColormap',parula);

% ---------- Return helpers ----------
pub.size         = @sizefig;           % apply size to CURRENT figure
pub.save         = @save_all_open;     % save ALL open figures
pub.savecurrent  = @save_current_only; % save only CURRENT figure
pub.autoname     = @save_autoname;     % auto-name CURRENT figure
pub.dir          = plotsDir;
pub.sizes        = sizes;

% ---------- Nested: sizing ----------
    function sizefig(which)
        if nargin<1 || isempty(which), which = 'single'; end
        s = getappdata(0,'PUBSIZE');
        if ~isstruct(s) || ~isfield(s,which), s = sizes; end
        sz = s.(which);
        set(gcf,'Units','centimeters','Position',[2 2 sz], ...
                 'PaperUnits','centimeters','PaperPositionMode','auto');
    end

% ---------- Nested: SAVE HELPERS ----------
    function save_all_open(basename,which,dpi)
        if nargin<1 || isempty(basename), basename = 'figure'; end
        if nargin<2 || isempty(which),    which    = 'single'; end
        if nargin<3 || isempty(dpi),      dpi      = 300; end
        assert(any(strcmp(which,{'single','double'})),'which must be ''single'' or ''double''');
        assert(isnumeric(dpi)&&isscalar(dpi)&&dpi>0,'dpi must be positive scalar');

        figs = findall(0,'Type','figure');
        if isempty(figs), warning('No open figures to save.'); return; end
        [~,idx] = sort([figs.Number]); figs = figs(idx);
        used = containers.Map('KeyType','char','ValueType','logical');

        for k = 1:numel(figs)
            fh = figs(k);
            figure(fh);
            sizefig(which);

            nm = get(fh,'Name');
            if ~isempty(nm), base = sanitize(nm);
            else,            base = sprintf('%s_%02d', sanitize(basename), k);
            end
            if isKey(used,base), base = sprintf('%s_%02d', base, k); end
            used(base) = true;

            base = uniquify(base, getappdata(0,'PUB_PLOTS_DIR'));
            save_one_fig(fh, base, dpi);
        end
    end

    function save_current_only(name,which,dpi)
        if nargin<1 || isempty(name), name = 'figure'; end
        if nargin<2 || isempty(which), which = 'single'; end
        if nargin<3 || isempty(dpi),   dpi   = 300; end
        assert(any(strcmp(which,{'single','double'})),'which must be ''single'' or ''double''');
        assert(isnumeric(dpi)&&isscalar(dpi)&&dpi>0,'dpi must be positive scalar');

        sizefig(which);
        name = uniquify(sanitize(name), getappdata(0,'PUB_PLOTS_DIR'));
        save_one_fig(gcf, name, dpi);
    end

    function save_autoname(which,dpi)
        if nargin<1 || isempty(which), which = 'single'; end
        if nargin<2 || isempty(dpi),   dpi   = 300; end
        assert(any(strcmp(which,{'single','double'})),'which must be ''single'' or ''double''');

        sizefig(which);
        ax = gca;
        t  = get_text(ax, 'Title');
        xl = get_text(ax, 'XLabel');
        yl = get_text(ax, 'YLabel');

        parts = {};
        if ~isempty(t),  parts{end+1} = t; end %#ok<AGROW>
        if ~isempty(yl) || ~isempty(xl)
            if ~isempty(yl) && ~isempty(xl)
                parts{end+1} = sprintf('%s_vs_%s', yl, xl); %#ok<AGROW>
            elseif ~isempty(yl)
                parts{end+1} = yl; %#ok<AGROW>
            else
                parts{end+1} = xl; %#ok<AGROW>
            end
        end
        if isempty(parts), parts = {'figure'}; end

        base = sanitize(strjoin(parts,'__'));
        base = truncate_to(base, 100);
        base = uniquify(base, getappdata(0,'PUB_PLOTS_DIR'));
        save_one_fig(gcf, base, dpi);
    end

% ---------- Utilities ----------
    function save_one_fig(fh, base, dpi)
        out = @(ext) fullfile(getappdata(0,'PUB_PLOTS_DIR'),[base '.' ext]);
        set(fh,'Renderer','painters');
        try, savefig(fh, out('fig'));                               catch, warning('Could not save .fig: %s', base); end
        try, print(fh,'-dpdf','-painters', out('pdf'));             catch, warning('PDF save failed: %s', base); end
        try, print(fh,'-dpng',sprintf('-r%d',dpi), out('png'));     catch, warning('PNG save failed: %s', base); end
        try, print(fh,'-depsc','-painters', out('eps'));            catch, warning('EPS save failed: %s', base); end
        try, print(fh,'-dsvg','-painters', out('svg'));             catch, warning('SVG save failed: %s', base); end
    end

    function s = sanitize(str)
        s = char(str);
        s = regexprep(s,'[\\{}$^]','');
        s = regexprep(s,'[^\w\-.()%]','_');
        s = regexprep(s,'_+','_');
        s = strtrim(s);
        if isempty(s), s = 'figure'; end
    end

    function s = truncate_to(s, N)
        if numel(s) > N, s = s(1:N); end
    end

    function base = uniquify(base, dirpath)
        exts = {'fig','pdf','png','eps','svg'};
        k = 0; candidate = base;
        while true
            exists_any = false;
            for ii = 1:numel(exts)
                if exist(fullfile(dirpath,[candidate '.' exts{ii}]),'file')
                    exists_any = true; break;
                end
            end
            if ~exists_any, base = candidate; return; end
            k = k + 1;
            candidate = sprintf('%s_%02d', base, k);
        end
    end

    function txt = get_text(ax, which)
        obj = get(ax, which);
        if isprop(obj,'String'), txt = obj.String; else, txt = ''; end
        if isstring(txt), txt = strjoin(cellstr(txt),' ');
        elseif iscellstr(txt), txt = strjoin(txt,' ');
        elseif ~ischar(txt), txt = '';
        end
        txt = strtrim(txt);
    end
end

% ---------- Theme-aware palettes ----------
function [bg, fg, gridCol, minorGridCol, palette] = theme_palette(name)
% Resolve theme name -> background / foreground / grid colours and a
% categorical palette tuned for that background.
    switch lower(name)
        case 'dark'
            bg           = [0.12 0.12 0.12];
            fg           = [0.92 0.92 0.92];
            gridCol      = [0.80 0.80 0.80];
            minorGridCol = [0.60 0.60 0.60];
            palette = [0.30 0.70 1.00;     % bright blue
                       0.95 0.55 0.30;     % bright orange
                       0.30 0.85 0.55;     % bright teal-green
                       0.85 0.55 0.95;     % bright purple
                       1.00 0.85 0.30;     % bright gold
                       0.55 0.85 1.00;     % light cyan
                       1.00 0.45 0.50;     % bright coral-red
                       0.85 0.85 0.85];    % light grey
        otherwise   % 'light'
            bg           = [1.00 1.00 1.00];
            fg           = [0.15 0.15 0.15];
            gridCol      = [0.15 0.15 0.15];
            minorGridCol = [0.30 0.30 0.30];
            palette = [0.000 0.447 0.741;  % blue
                       0.850 0.325 0.098;  % vermillion
                       0.000 0.620 0.450;  % teal-green
                       0.494 0.184 0.556;  % purple
                       0.929 0.694 0.125;  % deep gold
                       0.301 0.745 0.933;  % sky blue
                       0.635 0.078 0.184;  % maroon
                       0.250 0.250 0.250]; % dark grey
    end
end
