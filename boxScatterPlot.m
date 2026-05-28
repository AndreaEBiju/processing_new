function [ax, h] = boxScatterPlot(data, varargin)
%BOXSCATTERPLOT  Box plots with scatter overlay, optionally grouped.
%
%   ax = boxScatterPlot(data) renders one box per cell, with all data
%   points overlaid as a jittered scatter.
%
%   data can be either:
%     - a cell vector  1 x G  (or G x 1)      -> G ungrouped boxes
%     - a cell matrix  G x S                  -> G groups of S sub-boxes
%       Each cell holds a numeric vector (NaNs are dropped). Rows are
%       the primary groups (e.g. cond1/cond2/cond3) and columns are the
%       sub-conditions to compare within each group (e.g. 1min/2min/5min).
%
%   [ax, h] = boxScatterPlot(...) also returns a struct of handles
%   (h.box, h.scatter, h.pair) for downstream tweaking.
%
%   Name-value options:
%     'GroupLabels'    cellstr length G,  x-axis tick labels       []
%     'SubgroupLabels' cellstr length S,  legend entries (if S>1)  []
%     'YLabel'         y-axis label                                 ''
%     'XLabel'         x-axis label                                 ''
%     'Title'          figure title                                 ''
%     'Colors'         S x 3 RGB, one row per subgroup (or 1 x 3
%                      for ungrouped). Defaults to Okabe-Ito.
%     'BoxAlpha'       face alpha of the boxes                      0.08
%     'BoxWidth'       width of a single box                        0.6 (S=1)
%                                                                   0.7/S (S>1)
%     'GroupGap'       gap (in x-units) between primary groups      1.0
%     'JitterWidth'    horizontal jitter of scatter points          0.6*BoxWidth
%     'MarkerSize'     scatter marker size                          28
%     'MarkerAlpha'    scatter marker face alpha                    0.85
%     'ShowOutliers'   show boxchart's outlier markers              false
%                      (scatter shows every point anyway)
%     'ConnectPaired'  draw faint lines connecting matching indices
%                      across subgroups within each group           false
%     'ColorBySubject' colour scatter dots by subject index so the
%                      same shade across S boxes within a group is
%                      the same row of data. Box face colour still
%                      encodes subgroup.                             false
%     'CycleMarkers'   cycle scatter marker shape per subject. Off by
%                      default; turn on (true) for >10 subjects where
%                      shade alone is not enough to identify subjects.  false
%     'SubjectColormap' 'auto' uses a lightness ramp of each box's
%                       own colour (so dots share the box hue, only
%                       differing in shade per subject -- recommended
%                       for scientific figures). Otherwise pass a
%                       colormap name ('gray','copper','parula',...)
%                       or an N x 3 RGB matrix.                       'auto'
%     'Notch'          'on' | 'off'                                  'off'
%     'Parent'         target axes handle                           gca
%     'Legend'         'auto' | 'on' | 'off'                         'auto'
%
%   Example (ungrouped, 3 conditions):
%     d = {randn(20,1), 0.3+randn(20,1), -0.4+0.8*randn(20,1)};
%     boxScatterPlot(d, 'GroupLabels',{'cond1','cond2','cond3'}, ...
%                       'YLabel','HR (bpm)','Title','Heart rate');
%
%   Example (grouped, 3 conditions x 3 time windows):
%     d = cell(3,3);
%     for g=1:3, for s=1:3, d{g,s} = g + 0.2*s + 0.4*randn(15,1); end; end
%     boxScatterPlot(d, ...
%        'GroupLabels',{'cond1','cond2','cond3'}, ...
%        'SubgroupLabels',{'1 min','2 min','5 min'}, ...
%        'YLabel','\Delta HR (bpm)','ConnectPaired',true);

% ---------- Normalise input ----------
if isnumeric(data) && isvector(data)
    data = {data(:)};
end
assert(iscell(data), 'data must be a cell array of numeric vectors.');
if isvector(data), data = data(:).'; end           % row -> 1 x G  treat as ungrouped
% By spec rows are groups, cols are subgroups. A 1xG row vector is one
% group with G subgroups, but users almost always mean G ungrouped boxes.
% Detect intent: if it's a vector AND no SubgroupLabels supplied, treat as G x 1.
[G, S] = size(data);
vectorInput = (G == 1 || S == 1);
if vectorInput
    data = reshape(data, [], 1);                  % G x 1
    [G, S] = size(data);
end

% ---------- Parse options ----------
p = inputParser;
p.FunctionName = 'boxScatterPlot';
addParameter(p,'GroupLabels',{},    @(x)iscellstr(x)||isstring(x));
addParameter(p,'SubgroupLabels',{}, @(x)iscellstr(x)||isstring(x));
addParameter(p,'YLabel','',         @(s)ischar(s)||isstring(s));
addParameter(p,'XLabel','',         @(s)ischar(s)||isstring(s));
addParameter(p,'Title','',          @(s)ischar(s)||isstring(s));
addParameter(p,'Colors',[],         @(x)isnumeric(x)&&size(x,2)==3);
addParameter(p,'BoxAlpha',0.08,     @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'BoxWidth',[],       @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'GroupGap',1.0,      @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'JitterWidth',[],    @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
addParameter(p,'MarkerSize',28,     @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MarkerAlpha',0.85,  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'ShowOutliers',false,@islogical);
addParameter(p,'ConnectPaired',false,@islogical);
addParameter(p,'ColorBySubject',false,@islogical);
addParameter(p,'CycleMarkers',false, ...
    @(x)islogical(x)||(ischar(x)||isstring(x)));
addParameter(p,'SubjectColormap','auto', ...
    @(x)ischar(x)||isstring(x)||(isnumeric(x)&&size(x,2)==3));
addParameter(p,'Notch','off',       @(s)any(strcmpi(s,{'on','off'})));
addParameter(p,'Parent',[],         @(x)isempty(x)||isgraphics(x,'axes'));
addParameter(p,'Legend','auto',     @(s)any(strcmpi(s,{'auto','on','off'})));
parse(p, varargin{:});
opt = p.Results;

% Resolve CycleMarkers (accept string for compatibility; off by default)
if ischar(opt.CycleMarkers) || isstring(opt.CycleMarkers)
    opt.CycleMarkers = strcmpi(char(opt.CycleMarkers),'on') ...
                    || strcmpi(char(opt.CycleMarkers),'true');
end
markerSet = {'o','s','^','d','v'};        % cycled per subject when CycleMarkers=true

% ---------- Axes ----------
if isempty(opt.Parent)
    ax = gca;
else
    ax = opt.Parent;
end
holdState = ishold(ax); hold(ax,'on');

% ---------- Theme foreground (paired-line / edge colours) ----------
fgCol = getappdata(0,'PUB_FGCOL');
if isempty(fgCol) || numel(fgCol) ~= 3
    fgCol = [0.15 0.15 0.15];                % light-theme fallback
end

% ---------- Colors ----------
if isempty(opt.Colors)
    C = okabeIto();
    if S == 1
        colors = C(1,:);
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
    if S == 1, boxW = 0.6; else, boxW = 0.7 / S; end
else
    boxW = opt.BoxWidth;
end
if isempty(opt.JitterWidth)
    jitW = 0.6 * boxW;
else
    jitW = opt.JitterWidth;
end
gap = opt.GroupGap;

groupCenter = (1:G) * gap;                % x location of each group
subOffset   = ((1:S) - (S+1)/2) * boxW;   % offset of each subgroup within a group

% ---------- Draw ----------
h.box     = gobjects(G,S);
h.scatter = gobjects(G,S);
h.pair    = cell(G,1);

rng_state = rng;                          % preserve caller's RNG
rng(0,'twister');                         % reproducible jitter

for g = 1:G
    % per-subject palette for this group
    if opt.ColorBySubject
        nMaxGroup = max(cellfun(@(v) numel(v), data(g,:)));
        useAutoShade = (ischar(opt.SubjectColormap) || isstring(opt.SubjectColormap)) ...
                        && strcmpi(char(opt.SubjectColormap),'auto');
        if useAutoShade
            cmapSubj = [];                              % built per-s from box colour
        else
            cmapSubj = makeSubjectCmap(opt.SubjectColormap, nMaxGroup);
        end
    else
        useAutoShade = false;
        cmapSubj     = [];
        nMaxGroup    = 0;
    end

    % paired-line layer first (so it sits behind boxes/scatter)
    if opt.ConnectPaired && S >= 2
        % collect a value matrix N x S for this group, NaN where missing
        nMax = max(cellfun(@(v) numel(v), data(g,:)));
        Y = nan(nMax, S);
        for s = 1:S
            v = data{g,s}(:);
            Y(1:numel(v), s) = v;
        end
        Xc = groupCenter(g) + subOffset;          % 1 x S box centers
        % per-subject jittered x (same jitter across all S so lines stay tidy)
        jx = (rand(nMax,1) - 0.5) * jitW;
        Xmat = bsxfun(@plus, Xc, jx);             % nMax x S
        valid = all(~isnan(Y), 2);
        if any(valid)
            h.pair{g} = plot(ax, Xmat(valid,:).', Y(valid,:).', ...
                '-', 'Color', [fgCol(:).' 0.35], 'LineWidth', 0.6, ...
                'HandleVisibility','off');
        end
    end

    for s = 1:S
        y_raw = data{g,s}(:);
        keep  = ~isnan(y_raw);
        y     = y_raw(keep);
        subjIdx = find(keep);            % original subject indices for ColorBySubject
        if isempty(y), continue; end
        xc = groupCenter(g) + subOffset(s);

        % box
        h.box(g,s) = boxchart(ax, xc*ones(size(y)), y, ...
            'BoxWidth',      boxW * 0.95, ...
            'BoxFaceColor',  colors(s,:), ...
            'BoxFaceAlpha',  opt.BoxAlpha, ...
            'WhiskerLineColor', colors(s,:)*0.6, ...
            'LineWidth',     1.0, ...
            'MarkerStyle',   ternary(opt.ShowOutliers,'o','none'), ...
            'MarkerColor',   colors(s,:), ...
            'Notch',         opt.Notch);

        % put one box per subgroup in the legend (the one in group 1)
        if g == 1 && S > 1
            set(h.box(g,s),'DisplayName', char(getLabel(opt.SubgroupLabels, s, sprintf('sub %d', s))));
        else
            set(h.box(g,s),'HandleVisibility','off');
        end

        % scatter overlay
        if opt.ConnectPaired && S >= 2
            % reuse the same jitter as the paired lines so dots sit on the joins
            % rebuild here by re-seeding deterministically per group is unnecessary;
            % we recompute fresh jitter that is independent of the paired layer
            xs = xc + (rand(numel(y),1) - 0.5) * (jitW * 0.4);
        else
            xs = xc + (rand(numel(y),1) - 0.5) * jitW;
        end
        if opt.ColorBySubject
            if useAutoShade
                ramp = shadeRamp(colors(s,:), nMaxGroup);
                faceCol = ramp(subjIdx,:);
            else
                faceCol = cmapSubj(subjIdx,:);
            end
            edgeCol = colors(s,:) * 0.2;          % near-black tint of box hue
        else
            faceCol = colors(s,:);
            edgeCol = fgCol;                 % follows theme foreground
        end

        if opt.CycleMarkers
            % each subject gets a unique (marker shape, shade) tuple
            nM = numel(markerSet);
            lastH = gobjects(0);
            for shp = 1:nM
                sel = (mod(subjIdx-1, nM) + 1) == shp;
                if ~any(sel), continue; end
                if size(faceCol,1) == 1
                    fc = faceCol;
                else
                    fc = faceCol(sel,:);
                end
                lastH = scatter(ax, xs(sel), y(sel), opt.MarkerSize, fc, ...
                    'filled', ...
                    'Marker',           markerSet{shp}, ...
                    'MarkerEdgeColor',  edgeCol, ...
                    'MarkerFaceAlpha',  opt.MarkerAlpha, ...
                    'LineWidth',        0.4, ...
                    'HandleVisibility','off');
            end
            h.scatter(g,s) = lastH;               % only last handle tracked
        else
            h.scatter(g,s) = scatter(ax, xs, y, opt.MarkerSize, faceCol, ...
                'filled', ...
                'MarkerEdgeColor', edgeCol, ...
                'MarkerFaceAlpha', opt.MarkerAlpha, ...
                'LineWidth',       0.4, ...
                'HandleVisibility','off');
        end
    end
end

rng(rng_state);                           % restore caller RNG

% ---------- Cosmetics ----------
xlim(ax, [groupCenter(1) - (max(boxW,0.6)*S/2 + 0.5*gap), ...
          groupCenter(end) + (max(boxW,0.6)*S/2 + 0.5*gap)]);
set(ax,'XTick', groupCenter);
if ~isempty(opt.GroupLabels)
    set(ax,'XTickLabel', cellstr(opt.GroupLabels));
else
    set(ax,'XTickLabel', arrayfun(@(k)sprintf('%d',k), 1:G,'uni',0));
end

% Force literal interpretation on tick labels (group names may contain
% '_' or '%' which LaTeX would mangle if pubfig_setup set that default).
try, set(ax, 'TickLabelInterpreter','none'); catch, end %#ok<CTCH>

if ~isempty(opt.XLabel), xlabel(ax, char(opt.XLabel), 'Interpreter','none'); end
if ~isempty(opt.YLabel), ylabel(ax, char(opt.YLabel), 'Interpreter','none'); end
if ~isempty(opt.Title),  title (ax, char(opt.Title), 'Interpreter','none'); end

% legend
showLeg = strcmpi(opt.Legend,'on') || (strcmpi(opt.Legend,'auto') && S > 1);
if showLeg
    legend(ax,'show','Location','best','Interpreter','none');
end

if ~holdState, hold(ax,'off'); end
end

% ---------- helpers ----------
function s = getLabel(labels, k, fallback)
    if isempty(labels) || numel(labels) < k
        s = fallback;
    else
        if isstring(labels), s = char(labels(k));
        else,                s = labels{k};
        end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function C = shadeRamp(base, n)
% n-step shade ramp of a single hue, traversed in HSV: lightness (V)
% goes up while saturation (S) drops, giving a real dark-saturated to
% pale-tinted arc within the box hue. Designed for up to ~10 subjects
% where shade encoding reads cleanly.
    base = base(:).';
    if n <= 1
        C = base * 0.6;
        return;
    end
    hsvb = rgb2hsv(base);
    H    = hsvb(1);
    V    = linspace(0.25, 0.95, n).';        % dark -> light
    S    = linspace(1.00, 0.40, n).';        % saturated -> tinted
    C    = hsv2rgb([repmat(H, n, 1), S, V]);
end

function C = makeSubjectCmap(spec, n)
% Return an n x 3 RGB palette for per-subject colouring.
%   spec : colormap name string ('turbo','parula','hsv','viridis',...)
%          or an N x 3 RGB matrix.
%   n    : number of colours required.
% Avoids the extreme dark/light ends of named colormaps so dots stay
% visible against the box face.
    if n < 1, C = zeros(0,3); return; end
    if isnumeric(spec)
        base = spec;
    else
        name = char(spec);
        try
            base = feval(name, 256);
        catch
            base = parula(256);
        end
    end
    if n == 1
        C = base(round(size(base,1)/2),:);
        return;
    end
    t = linspace(0.1, 0.9, n);                          % trim the ends
    C = interp1(linspace(0,1,size(base,1)), base, t, 'linear', 'extrap');
end

function C = okabeIto()
% Default palette. Prefers the active pubfig_setup palette (stored in
% appdata) so the theme set globally flows through to box plots. Falls
% back to a light-theme palette if pubfig_setup was never called.
    C = getappdata(0,'PUB_PALETTE');
    if isempty(C) || size(C,2) ~= 3
        C = [0.000 0.447 0.741;     % blue
             0.850 0.325 0.098;     % vermillion / orange-red
             0.000 0.620 0.450;     % teal-green
             0.494 0.184 0.556;     % purple
             0.929 0.694 0.125;     % deep gold
             0.301 0.745 0.933;     % sky blue
             0.635 0.078 0.184;     % maroon
             0.250 0.250 0.250];    % dark grey
    end
end
