function f = plot_total_rate_mixed(out, saveDir)
% PLOT_TOTAL_RATE_MIXED  One relative-change-from-baseline heatmap for the TOTAL
% firing rate (LVN + RVN), using the same mixed-effects model + styling as
% bulk_mixed_models. Standalone -- does NOT regenerate the other metrics.
%
%   plot_total_rate_mixed(out)            % out = output of run_pipeline_bulk
%   plot_total_rate_mixed(out, saveDir)   % also saves png/svg/fig
%
% Total rate per recording = sum of the per-channel mean firing rate. Each
% recovery is normalized to its matched baseline: y = (total_rec - total_base)/
% total_base, one value per (animal, condition) trial. Model: no-intercept M/E
% presence indicators + interaction + (1|animal); each cell is that condition's
% modeled % change from baseline, FDR-corrected vs baseline.
%
% Needs me_heatmap_render.m (and save_one_figure.m for saving) on the path and
% the Statistics & ML Toolbox.

    if nargin < 2; saveDir = ''; end
    formats = {'png','svg','fig'};
    assert(exist('me_heatmap_render','file')==2, 'me_heatmap_render.m not on path.');

    T = build_total_rate_table(out.raw);
    if height(T) < 6 || numel(unique(T.animal)) < 3
        error('plot_total_rate_mixed: too few total-rate trials (%d).', height(T));
    end

    mainAll = {'m10','m50','m100','e10','e100','e1000'};
    main = mainAll(cellfun(@(v) any(T.(v)~=0), mainAll));
    interAll = {'m10:e10','m10:e100','m10:e1000','m50:e10','m50:e100','m50:e1000', ...
                'm100:e10','m100:e100','m100:e1000'};
    inter = interAll(cellfun(@(v) any(prodcol(T,v)~=0), interAll));
    fixedFull = ['y ~ -1 + ' strjoin([main inter],' + ')];

    obsPer = grpcount(T.animal); useRand = max(obsPer) >= 2;
    if useRand
        mdl = fitlme(T, [fixedFull ' + (1|animal)'], 'FitMethod','ML');
    else
        mdl = fitlm(T, fixedFull);
    end
    CoefName = mdl.CoefficientNames; Est = mdl.Coefficients.Estimate;

    mLevAll = [10 50 100]; eLevAll = [10 100 1000];
    has = @(nm) any(strcmp(CoefName, nm));
    mLev = mLevAll(arrayfun(@(L) has(sprintf('m%d',L)), mLevAll));
    eLev = eLevAll(arrayfun(@(L) has(sprintf('e%d',L)), eLevAll));
    nM = numel(mLev); nE = numel(eLev);
    Z = nan(nM+1, nE+1); P = nan(nM+1, nE+1);
    for j = 1:nE; [e,p] = csum(mdl,CoefName,Est,{sprintf('e%d',eLev(j))}); Z(1,1+j)=100*e; P(1,1+j)=p; end
    for i = 1:nM; [e,p] = csum(mdl,CoefName,Est,{sprintf('m%d',mLev(i))}); Z(1+i,1)=100*e; P(1+i,1)=p; end
    for i = 1:nM
        for j = 1:nE
            [e,p] = csum(mdl,CoefName,Est,{sprintf('m%d',mLev(i)),sprintf('e%d',eLev(j)),sprintf('m%d:e%d',mLev(i),eLev(j))});
            Z(1+i,1+j) = 100*e; P(1+i,1+j) = p;
        end
    end
    fin = isfinite(P); q = nan(size(P)); if any(fin(:)); q(fin) = bhfdr(P(fin)); end
    Sg = isfinite(q) & q < 0.05;

    f = me_heatmap_render(Z, mLev, eLev, ...
        'Total firing rate (LVN+RVN) : relative change from baseline (%)', ...
        'relative change from baseline (%)   (* FDR<0.05 vs baseline)', Sg);
    if ~isempty(saveDir) && exist('save_one_figure','file')==2
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        save_one_figure(f, saveDir, 'mixedmodel_Total_firing_rate', formats);
    end
    fprintf('Total firing rate: %d trials, %d animals, %s.\n', ...
        height(T), numel(unique(T.animal)), class(mdl));
end

% ======================================================================
function T = build_total_rate_table(rows)
    n = numel(rows);
    if isfield(rows,'stem') && ~all(cellfun(@isempty,{rows.stem}))
        rk = arrayfun(@(r) sprintf('%s|%s', char(rows(r).stem), rows(r).phase), 1:n, 'uni', 0);
    else
        rk = arrayfun(@(r) sprintf('%s|%s|%s', rows(r).animal, rows(r).condition, rows(r).phase), 1:n, 'uni', 0);
    end
    uk = unique(rk, 'stable');
    rec = struct('phase',{},'animal',{},'condition',{},'total',{});
    for k = 1:numel(uk)
        idx = find(strcmp(rk, uk{k})); tot = 0; cnt = 0;
        for ii = idx
            if isfield(rows(ii),'mean') && isfield(rows(ii).mean,'rate') && isfinite(rows(ii).mean.rate)
                tot = tot + rows(ii).mean.rate; cnt = cnt + 1;   % sum channels in this recording
            end
        end
        if cnt == 0; continue; end
        rec(end+1) = struct('phase',rows(idx(1)).phase,'animal',rows(idx(1)).animal, ...
            'condition',rows(idx(1)).condition,'total',tot); %#ok<AGROW>
    end
    bi = find(strcmpi({rec.phase},'baseline'));
    y=[]; an={}; Mv=[]; Ev=[];
    for r = find(strcmpi({rec.phase},'recovery'))
        a = rec(r).animal; c = rec(r).condition;
        cand = bi(strcmpi({rec(bi).animal},a) & strcmpi({rec(bi).condition},c));
        if isempty(cand); continue; end
        bt = rec(cand(1)).total; rt = rec(r).total;
        if ~isfinite(bt) || bt == 0; continue; end
        [M,E] = parse_me(c);
        y(end+1,1)=(rt-bt)/bt; an{end+1,1}=a; Mv(end+1,1)=M; Ev(end+1,1)=E; %#ok<AGROW>
    end
    T = table(y(:), categorical(cellstr(string(an(:)))), 'VariableNames', {'y','animal'});
    T.m10=double(Mv(:)==10); T.m50=double(Mv(:)==50); T.m100=double(Mv(:)==100);
    T.e10=double(Ev(:)==10); T.e100=double(Ev(:)==100); T.e1000=double(Ev(:)==1000);
end

function c = prodcol(T,name); p = strsplit(name,':'); c = T.(p{1}) .* T.(p{2}); end

function [e,p] = csum(mdl, CoefName, Est, names)
    L = zeros(1, numel(CoefName)); e = 0;
    for nn = names
        idx = find(strcmp(CoefName, nn{1}), 1);
        if isempty(idx) && contains(nn{1}, ':')
            pp = strsplit(nn{1},':'); idx = find(strcmp(CoefName, [pp{2} ':' pp{1}]), 1);
        end
        if isempty(idx); e = NaN; p = NaN; return; end
        L(idx) = 1; e = e + Est(idx);
    end
    try, p = coefTest(mdl, L); catch, p = NaN; end %#ok<CTCH>
end

function [M,E] = parse_me(cond)
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end

function nO = grpcount(g); u = unique(g); nO = zeros(numel(u),1); for i=1:numel(u); nO(i)=sum(g==u(i)); end; end

function q = bhfdr(p)
    p = p(:); m = numel(p); [ps,ix] = sort(p); q = nan(m,1);
    qs = ps .* m ./ (1:m)'; for i=m-1:-1:1; qs(i)=min(qs(i),qs(i+1)); end
    q(ix) = min(qs,1);
end
