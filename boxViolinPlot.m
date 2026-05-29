function [ax, h] = boxViolinPlot(data, varargin)
%BOXVIOLINPLOT  Kernel-density violin + box overlay, optionally grouped.
%
%   ax = boxViolinPlot(data) renders one box + symmetric KDE violin per
%   non-empty cell. Same input shape rules as boxScatterPlot:
%     - cell vector 1xG / Gx1  -> G ungrouped violins
%     - cell matrix G x S      -> G groups of S sub-violins
%
%   [ax, h] = boxViolinPlot(...) returns a struct of handles:
%     h.violin (G x S patch handles)
%     h.box    (G x S boxchart handles)
%     h.scatter(G x S scatter handles; gobjects if ShowScatter=false)
%
%   Name-value options
%     'GroupLabels'    cellstr, x tick labels                      []
%     'SubgroupLabels' cellstr, legend entries (S > 1)             []
%     'YLabel'         y-axis label                                ''
%     'XLabel'         x-axis label                                ''
%     'Title'          figure title                                ''
%     'Colors'         S x 3 RGB                                   Okabe-Ito-ish
%     'BoxAlpha'       face alpha of the inner box                 0.20
%     'ViolinAlpha'    face alpha of the violin patch              0.30
%     'BoxWidth'       width of one violin (== one subgroup slot)  see below
%     'BoxInnerFrac'   inner box width as a fraction of BoxWidth   0.30
%     'GroupGap'       x-spacing between primary groups            1.0
%     'ShowScatter'    overlay individual data points              false
%     'JitterWidth'    scatter jitter                              0.4*BoxWidth
%     'MarkerSize'     scatter marker size                         28
%     'MarkerAlpha'    scatter marker alpha                        0.7
%     'Notch'          'on' | 'off' (passed to boxchart)           'off'
%     'KdeBandwidth'   ksdensity Bandwidth ([] = auto)             []
%     'KdePoints'      KDE resolution                              100
%     'Parent'         target axes handle                          gca
%     'Legend'         'auto' | 'on' | 'off'                       'auto'
%
%   Requires Statistics and Machine Learning Toolbox for ksdensity.
%   Cells with fewer than 2 finite values fall back to a horizontal
%   tick at that y-value (no violin can be estimated from 1 point).

% ---------- Normalise input ----------
if isnumeric(data) && isvector(data), data = {data(:)}; end
assert(iscell(data), 'data must be a cell array of numeric vectors.');
if isvector(data), data = data(:); end
[G, S] = size(data);

% ---------- Parse options ----------
p = inputParser;
p.FunctionName = 'boxViolinPlot';
addParameter(p,'GroupLabels',{},    @(x)iscellstr(x)||isstring(x));
addParameter(p,'SubgroupLabels',{}, @(x)iscellstr(x)||isstring(x));
addParameter(p,'YLabel','',         @(s)ischar(s)||isstring(s));
addParameter(p,'XLabel','',         @(s)ischar(s)||isstring(s));
addParameter(p,'Title','',          @(s)ischar(s)||isstring(s));
addParameter(p,'Colors',[],         @(x)isnumeric(x)&&size(x,2)==3);
addParameter(p,'BoxAlpha',0.20,     @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'ViolinAlpha',0.30,  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'BoxWidth',[],       @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'BoxInnerFrac',0.30, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'GroupGap',1.0,      @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'ShowScatter',false, @islogical);
addParameter(p,'JitterWidth',[],    @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(p,'MarkerSize',28,     @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MarkerAlpha',0.7,   @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'Notch','off',       @(s)any(strcmpi(s,{'on','off'})));
addParameter(p,'KdeBandwidth',[],   @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'KdePoints',100,     @(x)isnumeric(x)&&isscalar(x)&&x>=8);
addParameter(p,'Parent',[],         @(x)isempty(x)||isgraphics(x,'axes'));
addParameter(p,'Legend','auto',     @(s)any(strcmpi(s,{'auto','on','off'})));
parse(p, varargin{:});
opt = p.Results;

% ---------- Axes ----------
if isempty(opt.Parent), ax = gca; else, ax = opt.Parent; end
holdState = ishold(ax); hold(ax,'on');

% ---------- Theme foreground (for tick / edge fallback) ----------
fgCol = getappdata(0,'PUB_FGCOL');
if isempty(fgCol) || numel(fgCol) ~= 3, fgCol = [0.15 0.15 0.15]; end

% ---------- Colors ----------
if isempty(opt.Colors)
    C = okabeIto();
    if S == 1, colors = C(1,:);
    else
        reps = ceil(S / size(C,1));
        colors = repmat(C, reps, 1);
        colors = colors(1:S,:);
    end
else
    colors = opt.Colors;
    if size(colors,1) < S
        colors = repmat(colors, ceil(S/size(colors,1)), 1);
        colors = colors(1:S,:);
    end
end

% ---------- Geometry ----------
if isempty(opt.BoxWidth)
    if S == 1, boxW = 0.65; else, boxW = 0.8 / S; end
else
    boxW = opt.BoxWidth;
end
if isempty(opt.JitterWidth), jitW = 0.4 * boxW; else, jitW = opt.JitterWidth; end
gap = opt.GroupGap;

groupCenter = (1:G) * gap;
subOffset   = ((1:S) - (S+1)/2) * boxW;

% ---------- Draw ----------
h.violin  = gobjects(G,S);
h.box     = gobjects(G,S);
h.scatter = gobjects(G,S);

rng_state = rng; rng(0,'twister');  % deterministic jitter

for g = 1:G
    for s = 1:S
        y = data{g,s}(:);
        y = y(isfinite(y));
        if isempty(y), continue; end
        xc = groupCenter(g) + subOffset(s);

        % ----- KDE violin (needs >= 2 distinct values) -----
        if numel(y) >= 2 && numel(unique(y)) >= 2
            try
                if isempty(opt.KdeBandwidth)
                    [f, xi] = ksdensity(y, 'NumPoints', opt.KdePoints);
                else
                    [f, xi] = ksdensity(y, 'NumPoints', opt.KdePoints, ...
                        'Bandwidth', opt.KdeBandwidth);
                end
                % scale density to half the violin slot width
                if max(f) > 0
                    fScaled = f(:) / max(f) * (boxW * 0.48);
                else
                    fScaled = zeros(size(f(:)));
                end
                xi = xi(:);
                px = [xc + fScaled; flipud(xc - fScaled)];
                py = [xi;          flipud(xi)];
                h.violin(g,s) = patch(ax, px, py, colors(s,:), ...
                    'FaceAlpha', opt.ViolinAlpha, ...
                    'EdgeColor', colors(s,:)*0.55, ...
                    'LineWidth', 1.0, ...
                    'HandleVisibility','off');
            catch
                % ksdensity not available — skip the violin patch silently
            end
        end

        % ----- inner box (boxchart at this x position) -----
        h.box(g,s) = boxchart(ax, xc*ones(size(y)), y, ...
            'BoxWidth',         boxW * opt.BoxInnerFrac, ...
            'BoxFaceColor',     colors(s,:), ...
            'BoxFaceAlpha',     opt.BoxAlpha, ...
            'WhiskerLineColor', colors(s,:) * 0.45, ...
            'LineWidth',        1.6, ...
            'MarkerStyle',      'none', ...
            'Notch',            opt.Notch);

        % legend entry: one per subgroup (use the box from g==1)
        if g == 1 && S > 1
            set(h.box(g,s),'DisplayName', char(getLabel(opt.SubgroupLabels, s, sprintf('sub %d', s))));
        else
            set(h.box(g,s),'HandleVisibility','off');
        end

        % ----- optional scatter overlay -----
        if opt.ShowScatter
            xs = xc + (rand(numel(y),1) - 0.5) * jitW;
            h.scatter(g,s) = scatter(ax, xs, y, opt.MarkerSize, ...
                colors(s,:) * 0.85, 'filled', ...
                'MarkerEdgeColor', fgCol, ...
                'MarkerFaceAlpha', opt.MarkerAlpha, ...
                'LineWidth', 0.4, ...
                'HandleVisibility','off');
        end
    end
end

rng(rng_state);

% ---------- Cosmetics ----------
xlim(ax, [groupCenter(1) - (max(boxW,0.6)*S/2 + 0.5*gap), ...
          groupCenter(end) + (max(boxW,0.6)*S/2 + 0.5*gap)]);
set(ax, 'XTick', groupCenter);
if ~isempty(opt.GroupLabels)
    set(ax, 'XTickLabel', cellstr(opt.GroupLabels));
else
    set(ax, 'XTickLabel', arrayfun(@(k)sprintf('%d',k), 1:G,'uni',0));
end
try, set(ax, 'TickLabelInterpreter','none'); catch, end %#ok<CTCH>

if ~isempty(opt.XLabel), xlabel(ax, char(opt.XLabel), 'Interpreter','none'); end
if ~isempty(opt.YLabel), ylabel(ax, char(opt.YLabel), 'Interpreter','none'); end
if ~isempty(opt.Title),  title (ax, char(opt.Title),  'Interpreter','none'); end

showLeg = strcmpi(opt.Legend,'on') || (strcmpi(opt.Legend,'auto') && S > 1);
if showLeg, legend(ax,'show','Location','best','Interpreter','none'); end

if ~holdState, hold(ax,'off'); end
end


% ---------- helpers ----------
function s = getLabel(labels, k, fallback)
    if isempty(labels) || numel(labels) < k, s = fallback; return; end
    if isstring(labels), s = char(labels(k)); else, s = labels{k}; end
end

function C = okabeIto()
    C = getappdata(0,'PUB_PALETTE');
    if isempty(C) || size(C,2) ~= 3
        C = [0.000 0.447 0.741;
             0.850 0.325 0.098;
             0.000 0.620 0.450;
             0.494 0.184 0.556;
             0.929 0.694 0.125;
             0.301 0.745 0.933;
             0.635 0.078 0.184;
             0.250 0.250 0.250];
    end
end
