function R = bulk_mmc_models(queueFile, saveDir, opts)
% BULK_MMC_MODELS  Mechanical x Electrical mixed-effects test (animal random) for
% the gastric MMC metrics, with the SAME model, contrasts, FDR and heatmap
% rendering as bulk_mixed_models -- only the metrics/source differ.
%
%   R = bulk_mmc_models                       % gemsplots_queue.mat in pwd, display only
%   R = bulk_mmc_models(queueFile, saveDir)   % also writes CSVs + saves heatmaps
%   R = bulk_mmc_models(queueFile, saveDir, opts)
%
% METRICS (baseline-normalized % change, from compile_mmc):
%   rate  = gastric individual-firing rate
%   burst = gastric burst-episode rate
%   amp   = gastric firing peak amplitude
%
% MODEL (identical to bulk_mixed_models)
%   y = (mean recovery - mean baseline)/mean baseline, one value per
%   (animal, condition[, channel]). M/E entered as presence indicators with NO
%   intercept, so each M:E coefficient is that cell's synergy and the missing
%   (0,0) cell is fixed at 0. animal is a random intercept (repeated measures);
%   omnibus M x E via LRT (LME) or joint F (LM fallback). Every heatmap cell is
%   that condition's modeled relative change from baseline (%), per-cell p vs
%   baseline, BH-FDR within metric.
%
% opts.channelMode : 'mean' (default; average the 3 gastric channels -> one model
%   per metric) | 'separate' (model G1,G2,G3 individually -> 3 models per metric).
%
% Needs the Statistics & ML Toolbox (fitlme) and, on the path, me_heatmap_render
% + whiten_figure (+ save_one_figure to save). compile_mmc must be on the path.

    if nargin < 1 || isempty(queueFile); queueFile = fullfile(pwd,'gemsplots_queue.mat'); end
    if nargin < 2; saveDir = ''; end
    if nargin < 3; opts = struct(); end
    chMode = 'mean'; if isfield(opts,'channelMode')&&~isempty(opts.channelMode); chMode = lower(opts.channelMode); end
    formats = {'png','svg','fig'};
    assert(license('test','Statistics_Toolbox') || exist('fitlme','file')==2, ...
        'bulk_mmc_models: needs the Statistics & Machine Learning Toolbox (fitlme).');
    assert(exist('compile_mmc','file')==2, 'bulk_mmc_models: compile_mmc.m not on path.');
    assert(exist('me_heatmap_render','file')==2, 'bulk_mmc_models: me_heatmap_render.m not on path.');

    [~, norm] = compile_mmc(queueFile, opts);

    mets = {'rate','burst','amp'};
    labs = {'Gastric firing rate','Gastric burst rate','Gastric peak amplitude'};
    switch chMode
        case 'separate'; chans = {'G1','G2','G3'};
        otherwise;       chans = {'G'};   % channel-mean
    end

    summ = {}; cellRows = {};
    for k = 1:numel(mets)
        for ci = 1:numel(chans)
            T = build_mmc_table(norm, mets{k}, chans{ci});
            if strcmp(chMode,'separate'); lab = sprintf('%s (%s)', labs{k}, chans{ci});
            else;                          lab = labs{k}; end
            [s, c, hm] = fit_metric(T, lab, 'mmc');
            summ{end+1} = s; cellRows = [cellRows; c]; %#ok<AGROW>
            maybe_plot(hm, s.metric, saveDir, formats);
        end
    end

    S = vertcat(summ{:});
    Summary = struct2table(S);
    Cells   = struct2table(vertcat(cellRows{:}));

    fprintf('\n============= MMC  M x E interaction (animal random) =============\n');
    for i = 1:height(Summary)
        fprintf('%-26s [%s, n=%d, %d animals]  interaction %s = %.3g, df=%d, p=%.3g\n', ...
            Summary.metric{i}, Summary.method{i}, Summary.nObs(i), Summary.nAnimals(i), ...
            Summary.statName{i}, Summary.stat(i), Summary.df(i), Summary.pInteraction(i));
    end
    sig = Cells(Cells.pFDR < 0.05, :);
    if ~isempty(sig)
        fprintf('\nSignificant synergy cells (FDR<0.05):\n');
        for i = 1:height(sig)
            fprintf('   %-26s M%d x E%d : %+.1f%%  (p=%.3g, pFDR=%.3g)\n', ...
                sig.metric{i}, sig.M(i), sig.E(i), 100*sig.synergy(i), sig.p(i), sig.pFDR(i));
        end
    else
        fprintf('\nNo individual synergy cell survives FDR<0.05.\n');
    end

    R = struct('summary',Summary,'cells',Cells);

    if ~isempty(saveDir)
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        writetable(Summary, fullfile(saveDir,'mmc_mixedmodel_interaction_summary.csv'));
        writetable(Cells,   fullfile(saveDir,'mmc_mixedmodel_synergy_cells.csv'));
        fprintf('\n[saved] %s\n', fullfile(saveDir,'mmc_mixedmodel_interaction_summary.csv'));
        fprintf('[saved] %s\n',   fullfile(saveDir,'mmc_mixedmodel_synergy_cells.csv'));
    end
end

% ======================================================================
function T = build_mmc_table(norm, metric, channel)
% One baseline-normalized scalar per (animal, condition) for the chosen channel
% label; M/E parsed from the condition string. Mirrors build_nerve_table.
    sel = find(strcmpi({norm.label}, channel));
    y=[]; an={}; Mv=[]; Ev=[];
    for i = sel(:)'
        if ~isfield(norm(i).scalar, metric); continue; end
        v = norm(i).scalar.(metric); if ~isfinite(v); continue; end
        [M,E] = parse_me(norm(i).condition);
        y(end+1,1)=v; an{end+1,1}=norm(i).animal; Mv(end+1,1)=M; Ev(end+1,1)=E; %#ok<AGROW>
    end
    T = make_table(y, an, Mv, Ev);
end

% ======================================================================
% ---- everything below is copied verbatim from bulk_mixed_models so the MMC
%      models/heatmaps behave identically to the nerve ones (no useLog path;
%      MMC metrics are rates/amplitudes, not the skewed HRV ratios) -----------
function [s, cells, hm] = fit_metric(T, label, source)
    cells = {}; hm = [];
    s = struct('metric',label,'source',source,'method','-','nObs',height(T), ...
        'nAnimals',numel(unique(T.animal)),'statName','-','stat',NaN,'df',NaN,'pInteraction',NaN);
    if height(T) < 6 || numel(unique(T.animal)) < 3
        warning('bulk_mmc_models:tooFew','%s: too few observations (%d) -- skipped.', label, height(T));
        return;
    end

    mainAll = {'m10','m50','m100','e10','e100','e1000'};
    main = mainAll(cellfun(@(v) any(T.(v)~=0), mainAll));
    interAll = {'m10:e10','m10:e100','m10:e1000','m50:e10','m50:e100','m50:e1000', ...
                'm100:e10','m100:e100','m100:e1000'};
    inter = interAll(cellfun(@(v) any(prod_col(T,v)~=0), interAll));
    if isempty(inter)
        warning('bulk_mmc_models:noInter','%s: no combined (MxE) cells -- interaction not testable.', label);
        return;
    end

    fixedRed  = ['y ~ -1 + ' strjoin(main,' + ')];
    fixedFull = [fixedRed ' + ' strjoin(inter,' + ')];

    obsPer = groupcounts_local(T.animal);
    useRand = max(obsPer) >= 2;

    try
        if useRand
            lmeFull = fitlme(T, [fixedFull ' + (1|animal)'], 'FitMethod','ML');
            lmeRed  = fitlme(T, [fixedRed  ' + (1|animal)'], 'FitMethod','ML');
            cmp = compare(lmeRed, lmeFull);
            s.method='LME'; s.statName='chi2'; s.stat=cmp.LRStat(2);
            s.df=cmp.deltaDF(2); s.pInteraction=cmp.pValue(2);
            mdl = lmeFull;
            CoefName = lmeFull.CoefficientNames; Est = lmeFull.Coefficients.Estimate;
            SE = lmeFull.Coefficients.SE; PV = lmeFull.Coefficients.pValue;
        else
            lmFull = fitlm(T, fixedFull);
            cn = lmFull.CoefficientNames; isI = contains(cn, ':');
            L = zeros(sum(isI), numel(cn)); ii = find(isI);
            for q=1:numel(ii); L(q, ii(q)) = 1; end
            [p,F,df1,~] = coefTest(lmFull, L);
            s.method='LM'; s.statName='F'; s.stat=F; s.df=df1; s.pInteraction=p;
            mdl = lmFull;
            CoefName = lmFull.CoefficientNames; Est = lmFull.Coefficients.Estimate;
            SE = lmFull.Coefficients.SE; PV = lmFull.Coefficients.pValue;
        end
    catch ME
        warning('bulk_mmc_models:fit','%s: fit failed (%s)', label, ME.message); return;
    end

    isI = contains(CoefName, ':');
    idx = find(isI); pvals = PV(idx); pFDR = bh_fdr(pvals);
    for q = 1:numel(idx)
        [M,E] = parse_inter(CoefName{idx(q)});
        cells{end+1,1} = struct('metric',label,'source',source,'M',M,'E',E, ...
            'synergy',Est(idx(q)),'SE',SE(idx(q)),'p',pvals(q),'pFDR',pFDR(q)); %#ok<AGROW>
    end

    mLevAll = [10 50 100]; eLevAll = [10 100 1000];
    has = @(nm) any(strcmp(CoefName, nm));
    mLev = mLevAll(arrayfun(@(L) has(sprintf('m%d',L)), mLevAll));
    eLev = eLevAll(arrayfun(@(L) has(sprintf('e%d',L)), eLevAll));
    nMl = numel(mLev); nEl = numel(eLev);
    toPct = @(e) 100 * e;
    Z = nan(nMl+1, nEl+1); P = nan(nMl+1, nEl+1);
    for j = 1:nEl
        [e,p] = contrast_sum(mdl, CoefName, Est, {sprintf('e%d',eLev(j))});
        Z(1,1+j) = toPct(e); P(1,1+j) = p;
    end
    for i = 1:nMl
        [e,p] = contrast_sum(mdl, CoefName, Est, {sprintf('m%d',mLev(i))});
        Z(1+i,1) = toPct(e); P(1+i,1) = p;
    end
    for i = 1:nMl
        for j = 1:nEl
            [e,p] = contrast_sum(mdl, CoefName, Est, ...
                {sprintf('m%d',mLev(i)), sprintf('e%d',eLev(j)), sprintf('m%d:e%d',mLev(i),eLev(j))});
            Z(1+i,1+j) = toPct(e); P(1+i,1+j) = p;
        end
    end
    fin = isfinite(P); q = nan(size(P)); if any(fin(:)); q(fin) = bh_fdr(P(fin)); end
    Sg = isfinite(q) & q < 0.05;
    hm = struct('Mlev',mLev,'Elev',eLev,'Z',Z,'sig',Sg);
end

% ======================================================================
function [est, p] = contrast_sum(mdl, CoefName, Est, names)
    L = zeros(1, numel(CoefName)); est = 0;
    for n = names
        idx = findcoef(CoefName, n{1});
        if isempty(idx); est = NaN; p = NaN; return; end
        L(idx) = 1; est = est + Est(idx);
    end
    try, p = coefTest(mdl, L); catch, p = NaN; end %#ok<CTCH>
end

function idx = findcoef(CoefName, n)
    idx = find(strcmp(CoefName, n), 1);
    if isempty(idx) && contains(n, ':')
        pp = strsplit(n, ':'); idx = find(strcmp(CoefName, [pp{2} ':' pp{1}]), 1);
    end
end

% ======================================================================
function maybe_plot(hm, label, saveDir, formats)
% EVERY cell = that condition's relative change from baseline (%); '*' + bold
% border = FDR<0.05 vs baseline. Identical rendering to the nerve heatmaps.
    if isempty(hm); return; end
    interior = hm.Z(2:end, 2:end);
    if ~any(isfinite(interior(:))); return; end
    ttl = sprintf('%s : relative change from baseline (%%)', label);
    cbl = 'relative change from baseline (%)   (* FDR<0.05 vs baseline)';
    f = me_heatmap_render(hm.Z, hm.Mlev, hm.Elev, ttl, cbl, hm.sig);
    if ~isempty(saveDir) && exist('save_one_figure','file')==2
        save_one_figure(f, saveDir, ['mmc_mixedmodel_' regexprep(label,'[^\w]+','_')], formats);
    end
end

% ======================================================================
function T = make_table(y, an, Mv, Ev)
    T = table(y(:), categorical(cellstr(string(an(:)))), 'VariableNames', {'y','animal'});
    T.m10  = double(Mv(:)==10);  T.m50  = double(Mv(:)==50);  T.m100  = double(Mv(:)==100);
    T.e10  = double(Ev(:)==10);  T.e100 = double(Ev(:)==100); T.e1000 = double(Ev(:)==1000);
end

function c = prod_col(T, name)
    p = strsplit(name, ':'); c = T.(p{1}) .* T.(p{2});
end

function [M,E] = parse_inter(name)
    p = strsplit(name, ':');
    M = str2double(regexprep(p{1},'\D','')); E = str2double(regexprep(p{2},'\D',''));
end

function [M,E] = parse_me(cond)
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end

function n = groupcounts_local(g)
    u = unique(g); n = zeros(numel(u),1);
    for i=1:numel(u); n(i)=sum(g==u(i)); end
end

function q = bh_fdr(p)
    p = p(:); m = numel(p); [ps,ix] = sort(p); q = nan(m,1);
    qs = ps .* m ./ (1:m)'; for i=m-1:-1:1; qs(i)=min(qs(i),qs(i+1)); end
    q(ix) = min(qs,1);
end
