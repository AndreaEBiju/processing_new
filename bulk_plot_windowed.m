function fig = bulk_plot_windowed(rawRows, groups, metric, chLabel)
% BULK_PLOT_WINDOWED  Box+violin of one metric, one channel, per group, where
% each condition shows FIVE sub-distributions taken from increasing windows of
% the recovery recording: first 1, 2, 5, 10 min and the full recovery.
%
% Every recovery sample x is baseline-normalized as (x - mean(baseline))/mean(baseline),
% where mean(baseline) is the mean of that (animal,condition,channel) metric over
% the FULL baseline recording. Distributions are pooled across animals.
%
%   rawRows : out.raw from run_pipeline_bulk (needs per-sample values AND times --
%             reprocess once with forceRefresh if your cache predates timing)
%   metric  : 'rate' | 'vpp' | 'fwhm' | 'cv2'
%   chLabel : 'RVN' | 'LVN'
%
% Uses YOUR boxViolinPlot.m (grouped G x 5 mode).

    assert(exist('boxViolinPlot','file')==2, ['bulk_plot_windowed:noDep ' ...
        'boxViolinPlot.m not on path. addpath(<your processing_new folder>) first.']);

    win  = [60 120 300 600 Inf];
    wlab = {'1min','2min','5min','10min','full'};
    ylab = sprintf('%s  (rec-base)/base', upper(metric));

    nG = numel(groups); nc = min(4,nG); nr = ceil(nG/nc);
    fig = figure('Color','w','Name',sprintf('%s windowed — %s', metric, chLabel), ...
        'Position',[80 80 430*nc 340*nr]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('%s  |  %s  |  recovery windows, baseline-normalized', upper(metric), chLabel), ...
        'Interpreter','none');

    for g = 1:nG
        ax = nexttile;
        conds = groups{g};
        data = cell(numel(conds), numel(win));
        for i = 1:numel(conds)
            for s = 1:numel(win)
                data{i,s} = gather_windowed(rawRows, conds{i}, chLabel, metric, win(s));
            end
        end
        boxViolinPlot(data, 'Parent', ax, 'GroupLabels', conds, ...
            'SubgroupLabels', wlab, 'YLabel', ylab, 'Title', sprintf('group %d', g));
        yline(ax, 0, 'k:', 'HandleVisibility','off');
    end
    whiten_figure(fig);   % uniform style: white bg, black text, font 20
end

% ======================================================================
function vals = gather_windowed(rawRows, cond, chLabel, metric, win)
    vals = [];
    sel = strcmpi({rawRows.condition},cond) & strcmpi({rawRows.label},chLabel);
    animals = unique({rawRows(sel).animal});
    for a = 1:numel(animals)
        rec  = find_row(rawRows, animals{a}, cond, chLabel, 'recovery');
        base = find_row(rawRows, animals{a}, cond, chLabel, 'baseline');
        if isempty(rec) || isempty(base); continue; end
        [v, t] = metric_vt(rawRows(rec),  metric);
        bv     = metric_v (rawRows(base), metric);
        mu = mean(bv(isfinite(bv)), 'omitnan');                 % mean of FULL baseline
        if ~isfinite(mu) || mu==0 || isempty(v); continue; end
        keep = isfinite(v) & isfinite(t) & (t <= win);
        vals = [vals; (v(keep) - mu)/mu]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
end

function i = find_row(rawRows, animal, cond, chLabel, phase)
    i = find(strcmpi({rawRows.animal},animal) & strcmpi({rawRows.condition},cond) ...
           & strcmpi({rawRows.label},chLabel) & strcmpi({rawRows.phase},phase), 1);
end

function [v, t] = metric_vt(row, metric)
    switch lower(metric)
        case 'rate'; v = row.dist.rate;  t = getfield_or(row.dist,'rate_t',[]);
        case 'vpp';  v = row.dist.vpp;   t = getfield_or(row.dist,'spk_t',[]);
        case 'fwhm'; v = row.dist.fwhm;  t = getfield_or(row.dist,'spk_t',[]);
        case 'cv2';  v = row.dist.cv2;   t = getfield_or(row.dist,'cv2_t',[]);
        otherwise;   v = []; t = [];
    end
    v = v(:); t = t(:);
    n = min(numel(v), numel(t)); v = v(1:n); t = t(1:n);   % align values<->times
end

function v = metric_v(row, metric)
    switch lower(metric)
        case 'rate'; v = row.dist.rate;
        case 'vpp';  v = row.dist.vpp;
        case 'fwhm'; v = row.dist.fwhm;
        case 'cv2';  v = row.dist.cv2;
        otherwise;   v = [];
    end
    v = v(:);
end

function v = getfield_or(s, f, d)
    if isfield(s,f); v = s.(f); else; v = d; end
end
