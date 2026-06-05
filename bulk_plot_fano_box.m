function fig = bulk_plot_fano_box(normRows, groups, chLabel)
% BULK_PLOT_FANO_BOX  Box + scatter of the baseline-normalized Fano slope, one
% channel, laid out by group, using YOUR boxScatterPlot.m.
%
% The Fano slope is one scalar per (animal, condition), so each box holds the
% per-animal normalized values (rec-base)/base for that condition; the scatter
% shows the individual animals.

    assert(exist('boxScatterPlot','file')==2, ['bulk_plot_fano_box:noDep ' ...
        'boxScatterPlot.m not on path. addpath(<your processing_new folder>) first.']);

    nG = numel(groups); nc = min(4,nG); nr = ceil(nG/nc);
    fig = figure('Color','w','Name',sprintf('Fano slope norm box — %s', chLabel), ...
        'Position',[120 120 380*nc 320*nr]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('Fano slope  |  %s  |  baseline-normalized', chLabel), 'Interpreter','none');

    for g = 1:nG
        ax = nexttile;
        conds = groups{g};
        data = cell(1, numel(conds));
        for i = 1:numel(conds)
            data{i} = get_scalar(normRows, conds{i}, chLabel);
        end
        boxScatterPlot(data, 'Parent', ax, 'GroupLabels', conds, ...
            'YLabel', 'Fano slope  (rec-base)/base', 'Title', sprintf('group %d', g));
        yline(ax, 0, 'k:', 'HandleVisibility','off');
    end
    whiten_figure(fig);   % uniform style: white bg, black text, font 20
end

% ======================================================================
function vals = get_scalar(normRows, cond, chLabel)
    vals = [];
    for r = 1:numel(normRows)
        R = normRows(r);
        if strcmpi(R.condition,cond) && strcmpi(R.label,chLabel) && isfield(R.scalar,'fano')
            vals(end+1,1) = R.scalar.fano; %#ok<AGROW>
        end
    end
    vals = vals(isfinite(vals));
end
