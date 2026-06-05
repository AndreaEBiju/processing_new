function fig = bulk_plot_windowed_total(rawRows, groups)
% BULK_PLOT_WINDOWED_TOTAL  Box+violin of TOTAL firing rate x = (LVN+RVN)/s, per
% group, with five recovery windows (1, 2, 5, 10 min, full) per condition.
%
% Per (animal, condition): the per-bin RVN and LVN rates are summed on common
% bins to give x(t); baseline-normalized as (x_rec - mean(x_base))/mean(x_base),
% where mean(x_base) is the mean of the summed baseline rate. Pooled across animals.
%
% Uses YOUR boxViolinPlot.m. Needs rawRows with per-bin rate + bin times
% (reprocess once with forceRefresh if your cache predates timing).

    assert(exist('boxViolinPlot','file')==2, ['bulk_plot_windowed_total:noDep ' ...
        'boxViolinPlot.m not on path. addpath(<your processing_new folder>) first.']);

    win  = [60 120 300 600 Inf];
    wlab = {'1min','2min','5min','10min','full'};

    nG = numel(groups); nc = min(4,nG); nr = ceil(nG/nc);
    fig = figure('Color','w','Name','Total (LVN+RVN) rate windowed', ...
        'Position',[60 60 560*nc 430*nr]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact','TileSpacing','compact');
    title(tl, 'Total firing rate (LVN+RVN)/s  |  recovery windows, baseline-normalized', ...
        'Interpreter','none');

    for g = 1:nG
        ax = nexttile;
        conds = groups{g};
        data = cell(numel(conds), numel(win));
        for i = 1:numel(conds)
            for s = 1:numel(win)
                data{i,s} = gather_total(rawRows, conds{i}, win(s));
            end
        end
        boxViolinPlot(data, 'Parent', ax, 'GroupLabels', conds, ...
            'SubgroupLabels', wlab, 'YLabel', '(LVN+RVN)/s  (rec-base)/base', ...
            'Title', sprintf('group %d', g));
        yline(ax, 0, 'k:', 'HandleVisibility','off');
        try, ax.XAxis.TickLabelRotation = 30; catch, end   % keep condition labels legible
    end
    % keep ONE shared legend, placed outside the grid (declutters the panels)
    lgs = findall(fig,'Type','legend');
    for q = 2:numel(lgs); delete(lgs(q)); end
    if ~isempty(lgs)
        try, lgs(1).Layout.Tile = 'east'; catch
            try, lgs(1).Location = 'northeastoutside'; catch, end
        end
    end
    whiten_figure(fig);   % uniform style: white bg, black text, font 20
end

% ======================================================================
function vals = gather_total(rawRows, cond, win)
    vals = [];
    sel = strcmpi({rawRows.condition},cond);
    animals = unique({rawRows(sel).animal});
    for a = 1:numel(animals)
        rR = find_row(rawRows, animals{a}, cond, 'RVN', 'recovery');
        rL = find_row(rawRows, animals{a}, cond, 'LVN', 'recovery');
        bR = find_row(rawRows, animals{a}, cond, 'RVN', 'baseline');
        bL = find_row(rawRows, animals{a}, cond, 'LVN', 'baseline');
        if isempty(rR)||isempty(rL)||isempty(bR)||isempty(bL); continue; end
        [tRec, totRec]  = sum_rate(rawRows(rR), rawRows(rL));
        [~,    totBase] = sum_rate(rawRows(bR), rawRows(bL));
        mu = mean(totBase(isfinite(totBase)), 'omitnan');
        if ~isfinite(mu) || mu==0 || isempty(totRec); continue; end
        keep = isfinite(totRec) & isfinite(tRec) & (tRec <= win);
        vals = [vals; (totRec(keep) - mu)/mu]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
end

function [t, tot] = sum_rate(rowA, rowB)
    a = rowA.dist.rate(:); b = rowB.dist.rate(:);
    if isfield(rowA.dist,'rate_t'); t = rowA.dist.rate_t(:); else; t = (0:numel(a)-1)'; end
    n = min([numel(a), numel(b), numel(t)]);
    if n == 0; t = []; tot = []; return; end
    tot = a(1:n) + b(1:n); t = t(1:n);
end

function i = find_row(rawRows, animal, cond, chLabel, phase)
    i = find(strcmpi({rawRows.animal},animal) & strcmpi({rawRows.condition},cond) ...
           & strcmpi({rawRows.label},chLabel) & strcmpi({rawRows.phase},phase), 1);
end
