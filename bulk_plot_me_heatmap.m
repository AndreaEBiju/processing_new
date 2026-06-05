function fig = bulk_plot_me_heatmap(rawRows, metric, chLabel)
% BULK_PLOT_ME_HEATMAP  M x E heatmap of the mean percent change from baseline.
%
%   fig = bulk_plot_me_heatmap(rawRows, metric, chLabel)
%
% Each cell = average over the trials (animals) of that condition of the
% per-trial percent change
%       100 * ( mean(recovery) - mean(baseline) ) / mean(baseline)
% computed on the chosen metric and channel. The interior 3x3 is the M x E
% combinations; the marginals are placed OUTSIDE the axis labels:
%   * E-alone (E10/E100/E1000 alone) as the TOP row, above the E labels
%   * M-alone (M10/M50/M100 alone)   as the LEFT column, beside the M labels
%
%   metric  : 'rate' | 'vpp' | 'fwhm' | 'cv2'
%   chLabel : 'RVN' | 'LVN'

    conds = unique({rawRows.condition});
    [Mlev, Elev] = parse_ME(conds);
    eL = sort(unique(Elev(Elev>0 & Mlev==0 & ~isnan(Elev))));
    mL = sort(unique(Mlev(Mlev>0 & Elev==0 & ~isnan(Mlev))));
    if isempty(eL) || isempty(mL)
        warning('bulk_plot_me_heatmap:levels','Need M-alone and E-alone conditions.'); fig = []; return;
    end
    nM = numel(mL); nE = numel(eL);

    % layout: row1 = E-alone, rows 2..nM+1 = M levels; col1 = M-alone, cols 2..nE+1 = E levels
    Z = nan(nM+1, nE+1);
    for j = 1:nE; Z(1,1+j)   = cond_pct(rawRows, sprintf('E%d', eL(j)), chLabel, metric); end
    for i = 1:nM; Z(1+i,1)   = cond_pct(rawRows, sprintf('M%d', mL(i)), chLabel, metric); end
    for i = 1:nM
        for j = 1:nE
            Z(1+i,1+j) = cond_pct(rawRows, sprintf('M%dE%d', mL(i), eL(j)), chLabel, metric);
        end
    end

    fig = me_heatmap_render(Z, mL, eL, ...
        sprintf('%s  |  %s  :  mean %% change from baseline (avg over animals)', upper(metric), chLabel), ...
        'mean % change from baseline');
end

% ======================================================================
function v = cond_pct(rawRows, cond, chLabel, metric)
% mean over trials (animals) of 100*(mean(rec)-mean(base))/mean(base)
    mm = lower(metric); vals = [];
    sel = strcmpi({rawRows.condition},cond) & strcmpi({rawRows.label},chLabel);
    animals = unique({rawRows(sel).animal});
    for a = 1:numel(animals)
        bi = find(strcmpi({rawRows.animal},animals{a}) & strcmpi({rawRows.condition},cond) ...
                & strcmpi({rawRows.label},chLabel) & strcmpi({rawRows.phase},'baseline'), 1);
        ri = find(strcmpi({rawRows.animal},animals{a}) & strcmpi({rawRows.condition},cond) ...
                & strcmpi({rawRows.label},chLabel) & strcmpi({rawRows.phase},'recovery'), 1);
        if isempty(bi) || isempty(ri); continue; end
        bv = rawRows(bi).dist.(mm); rv = rawRows(ri).dist.(mm);
        bm = mean(bv(isfinite(bv))); rm = mean(rv(isfinite(rv)));
        if ~isfinite(bm) || bm==0 || ~isfinite(rm); continue; end
        vals(end+1) = 100*(rm-bm)/bm; %#ok<AGROW>
    end
    if isempty(vals); v = NaN; else; v = mean(vals); end
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
