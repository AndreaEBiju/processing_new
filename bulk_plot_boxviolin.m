function fig = bulk_plot_boxviolin(normRows, groups, metric, chLabel, yMax)
% BULK_PLOT_BOXVIOLIN  Box + violin of one baseline-normalized metric, one
% channel, laid out by group, using YOUR boxViolinPlot.m for styling.
%
%   normRows : compiled rows from bulk_compile (per-sample normalized,
%              (recovery - mean_baseline)/mean_baseline)
%   groups   : cell of groups; each a cell of condition names (defineGroupsUI form)
%   metric   : 'rate' | 'vpp' | 'fwhm' | 'cv2'
%   chLabel  : channel label ('RVN' | 'LVN')
%   yMax     : optional scalar; caps the upper y-axis limit on every panel
%              (lower limit stays automatic). [] / omitted = auto.
%
% Each violin holds the normalized recovery distribution for that condition,
% pooled across animals.

    if nargin < 5; yMax = []; end
    assert(exist('boxViolinPlot','file')==2, ['bulk_plot_boxviolin:noDep ' ...
        'boxViolinPlot.m not on path. addpath(<your processing_new folder>) first.']);

    ylab = sprintf('%s  (rec-base)/base', upper(metric));
    nG = numel(groups); nc = min(4,nG); nr = ceil(nG/nc);
    fig = figure('Color','w','Name',sprintf('%s norm box-violin — %s', metric, chLabel), ...
        'Position',[80 80 480*nc 400*nr]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('%s  |  %s  |  baseline-normalized', upper(metric), chLabel), 'Interpreter','none');

    for g = 1:nG
        ax = nexttile;
        conds = groups{g};
        data = cell(1, numel(conds));
        for i = 1:numel(conds)
            data{i} = get_dist(normRows, conds{i}, chLabel, metric);
        end
        boxViolinPlot(data, 'Parent', ax, 'GroupLabels', conds, ...
            'YLabel', ylab, 'Title', sprintf('group %d', g));
        yline(ax, 0, 'k:', 'HandleVisibility','off');   % baseline reference
        try, ax.XAxis.TickLabelRotation = 30; catch, end   % keep condition labels legible
        if ~isempty(yMax); ylim(ax, [-yMax yMax]); end   % symmetric y-axis cap about 0
    end
    whiten_figure(fig);   % uniform style: white bg, black text, font 20
end

% ======================================================================
function vals = get_dist(normRows, cond, chLabel, metric)
    vals = []; mm = lower(metric);
    for r = 1:numel(normRows)
        R = normRows(r);
        if strcmpi(R.condition,cond) && strcmpi(R.label,chLabel) && isfield(R.dist,mm)
            vals = [vals; R.dist.(mm)(:)]; %#ok<AGROW>
        end
    end
    vals = vals(isfinite(vals));
end
