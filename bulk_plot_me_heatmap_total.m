function fig = bulk_plot_me_heatmap_total(rawRows)
% BULK_PLOT_ME_HEATMAP_TOTAL  M x E mean-%-change heatmap of the TOTAL firing
% rate x = (LVN + RVN)/s. Per animal, the RVN and LVN mean rates are summed for
% baseline and for recovery, the per-animal percent change is taken, then
% averaged over animals. Marginals: E-alone top row, M-alone left column.

    conds = unique({rawRows.condition});
    [Mlev, Elev] = parse_ME(conds);
    eL = sort(unique(Elev(Elev>0 & Mlev==0 & ~isnan(Elev))));
    mL = sort(unique(Mlev(Mlev>0 & Elev==0 & ~isnan(Mlev))));
    if isempty(eL) || isempty(mL)
        warning('bulk_plot_me_heatmap_total:levels','Need M-alone and E-alone conditions.'); fig = []; return;
    end
    nM = numel(mL); nE = numel(eL);

    Z = nan(nM+1, nE+1);
    for j = 1:nE; Z(1,1+j) = cond_pct_total(rawRows, sprintf('E%d', eL(j))); end
    for i = 1:nM; Z(1+i,1) = cond_pct_total(rawRows, sprintf('M%d', mL(i))); end
    for i = 1:nM
        for j = 1:nE
            Z(1+i,1+j) = cond_pct_total(rawRows, sprintf('M%dE%d', mL(i), eL(j)));
        end
    end

    fig = me_heatmap_render(Z, mL, eL, ...
        'total rate (LVN+RVN)/s  :  mean % change from baseline (avg over animals)', ...
        'mean % change from baseline');
end

% ======================================================================
function v = cond_pct_total(rawRows, cond)
    vals = [];
    animals = unique({rawRows(strcmpi({rawRows.condition},cond)).animal});
    for a = 1:numel(animals)
        bR = row_idx(rawRows, animals{a}, cond, 'RVN', 'baseline');
        bL = row_idx(rawRows, animals{a}, cond, 'LVN', 'baseline');
        rR = row_idx(rawRows, animals{a}, cond, 'RVN', 'recovery');
        rL = row_idx(rawRows, animals{a}, cond, 'LVN', 'recovery');
        if isempty(bR)||isempty(bL)||isempty(rR)||isempty(rL); continue; end
        bm = mean_rate(rawRows(bR)) + mean_rate(rawRows(bL));
        rm = mean_rate(rawRows(rR)) + mean_rate(rawRows(rL));
        if ~isfinite(bm) || bm==0 || ~isfinite(rm); continue; end
        vals(end+1) = 100*(rm-bm)/bm; %#ok<AGROW>
    end
    if isempty(vals); v = NaN; else; v = mean(vals); end
end

function m = mean_rate(r)
    x = r.dist.rate; m = mean(x(isfinite(x)));
end

function i = row_idx(rawRows, animal, cond, chLabel, phase)
    i = find(strcmpi({rawRows.animal},animal) & strcmpi({rawRows.condition},cond) ...
           & strcmpi({rawRows.label},chLabel) & strcmpi({rawRows.phase},phase), 1);
end

function [Mlev, Elev] = parse_ME(conds)
    n = numel(conds); Mlev = nan(n,1); Elev = nan(n,1);
    for k = 1:n
        c = upper(strtrim(conds{k}));
        t = regexp(c,'^M(\d+)E(\d+)$','tokens','once');
        if ~isempty(t); Mlev(k)=str2double(t{1}); Elev(k)=str2double(t{2}); continue; end
        t = regexp(c,'^E(\d+)M(\d+)$','tokens','once');
        if ~isempty(t); Elev(k)=str2double(t{1}); Mlev(k)=str2double(t{2}); continue; end
        t = regexp(c,'^M(\d+)$','tokens','once');
        if ~isempty(t); Mlev(k)=str2double(t{1}); Elev(k)=0; continue; end
        t = regexp(c,'^E(\d+)$','tokens','once');
        if ~isempty(t); Elev(k)=str2double(t{1}); Mlev(k)=0; end
    end
end
