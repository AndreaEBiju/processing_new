function R = bulk_mixed_models(out, saveDir)
% BULK_MIXED_MODELS  Mixed-effects test of the Mechanical x Electrical interaction
% (synergy) for every metric, with animal as a random factor.
%
%   R = bulk_mixed_models(out)            % out = output of run_pipeline_bulk
%   R = bulk_mixed_models(out, saveDir)   % also writes CSV summaries
%
% MODEL
%   Response y = per-recording baseline-normalized change of the metric,
%               y = (mean(recovery) - mean(baseline)) / mean(baseline),
%               one value per (animal, condition[, channel]).
%   M and E are entered as PRESENCE indicators (m10,m50,m100 / e10,e100,e1000)
%   with NO intercept, so the missing (0,0) cell is fixed at 0 (true because y
%   is a % change). Each M:E coefficient is then the synergy of that combined
%   cell:  beta(MiEj) = mean(MiEj) - mean(Mi alone) - mean(Ej alone).
%
%   RVN and LVN are distinct nerves, so each nerve metric is modelled per
%   channel (RVN and LVN separately) -- NOT pooled. Each animal is measured
%   under multiple conditions (repeated measures), so the animal random
%   intercept (1|animal) is identifiable and the model is a genuine LME fit by
%   ML, with a likelihood-ratio test for the M x E interaction. (If a metric
%   happens to have one obs/animal, the code falls back to a fixed-effects LM
%   with a joint F-test -- auto-detected from obs-per-animal.)
%
%   Nerve responses are taken straight from out.norm (bulk_compile): one row per
%   recovery trial x channel, each paired to its own baseline (by trial stem),
%   repeat trials kept as SEPARATE observations, animal parsed by parse_stem --
%   the same numbers the box/violin and heatmap plots use. Systemic responses are
%   paired by (animal,condition) via the same convention as compute_norm.
%
%   Omnibus interaction: LME -> likelihood-ratio test (full vs additive, ML);
%   LM -> joint F-test on the interaction coefficients. Per-cell synergy
%   p-values are BH-FDR corrected within each metric.
%
% Requires the Statistics & ML Toolbox. Systemic metrics also need your
% buildFileQueue / loadAllSeriesCached on the path (same loaders as
% bulk_hrv_heatmaps); if absent they are skipped with a warning.

    if nargin < 2; saveDir = ''; end
    formats = {'png','svg','fig'};
    assert(license('test','Statistics_Toolbox') || exist('fitlme','file')==2, ...
        'bulk_mixed_models: needs the Statistics & Machine Learning Toolbox (fitlme).');

    summ = {};   % summary rows
    cellRows = {};   % per-cell synergy rows

    % ---------------- nerve metrics (from out.norm) ------------------------
    nMet = {'rate','vpp','fwhm','cv2'};
    nLab = {'Nerve firing rate','Vpp','FWHM','CV2'};
    nChan = {'RVN','LVN'};
    for k = 1:numel(nMet)
        for ci = 1:numel(nChan)   % RVN and LVN are distinct nerves -> separate models
            T = build_nerve_table(out.norm, nMet{k}, nChan{ci});
            [s, c, hm] = fit_metric(T, sprintf('%s (%s)', nLab{k}, nChan{ci}), 'nerve');
            summ{end+1} = s; cellRows = [cellRows; c]; %#ok<AGROW>
            maybe_plot(hm, s.metric, saveDir, formats);
        end
    end

    % ---------------- systemic metrics (HR / HRV loaders) ------------------
    haveLoaders = exist('buildFileQueue','file')==2 && exist('loadAllSeriesCached','file')==2;
    if haveLoaders
        try
            Tsys = build_systemic_tables();   % struct array .label, .T
            for k = 1:numel(Tsys)
                [s, c, hm] = fit_metric(Tsys(k).T, Tsys(k).label, 'systemic');
                summ{end+1} = s; cellRows = [cellRows; c]; %#ok<AGROW>
                maybe_plot(hm, s.metric, saveDir, formats);
            end
        catch ME
            warning('bulk_mixed_models:sys','systemic metrics skipped (%s)', ME.message);
        end
    else
        warning('bulk_mixed_models:sys', ...
            'buildFileQueue/loadAllSeriesCached not on path -- systemic metrics skipped.');
    end

    % ---------------- assemble + report ------------------------------------
    S = vertcat(summ{:});
    Summary = struct2table(S);
    Cells   = struct2table(vertcat(cellRows{:}));

    fprintf('\n================ M x E interaction (animal random) ================\n');
    for i = 1:height(Summary)
        fprintf('%-20s [%s, n=%d, %d animals]  interaction %s = %.3g, df=%d, p=%.3g\n', ...
            Summary.metric{i}, Summary.method{i}, Summary.nObs(i), Summary.nAnimals(i), ...
            Summary.statName{i}, Summary.stat(i), Summary.df(i), Summary.pInteraction(i));
    end
    sig = Cells(Cells.pFDR < 0.05, :);
    if ~isempty(sig)
        fprintf('\nSignificant synergy cells (FDR<0.05):\n');
        for i = 1:height(sig)
            fprintf('   %-20s M%d x E%d : %+.1f%%  (p=%.3g, pFDR=%.3g)\n', ...
                sig.metric{i}, sig.M(i), sig.E(i), 100*sig.synergy(i), sig.p(i), sig.pFDR(i));
        end
    else
        fprintf('\nNo individual synergy cell survives FDR<0.05.\n');
    end

    R = struct('summary',Summary,'cells',Cells);

    if ~isempty(saveDir)
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        writetable(Summary, fullfile(saveDir,'mixedmodel_interaction_summary.csv'));
        writetable(Cells,   fullfile(saveDir,'mixedmodel_synergy_cells.csv'));
        fprintf('\n[saved] %s\n', fullfile(saveDir,'mixedmodel_interaction_summary.csv'));
        fprintf('[saved] %s\n',   fullfile(saveDir,'mixedmodel_synergy_cells.csv'));
    end
end

% ======================================================================
function [s, cells, hm] = fit_metric(T, label, source)
    cells = {}; hm = [];
    s = struct('metric',label,'source',source,'method','-','nObs',height(T), ...
        'nAnimals',numel(unique(T.animal)),'statName','-','stat',NaN,'df',NaN,'pInteraction',NaN);
    if height(T) < 6 || numel(unique(T.animal)) < 3
        warning('bulk_mixed_models:tooFew','%s: too few observations (%d) -- skipped.', label, height(T));
        return;
    end

    hasCh = ismember('chLVN', T.Properties.VariableNames) && numel(unique(T.chLVN)) > 1;

    % which presence terms actually occur (avoid all-zero / rank-deficient cols)
    mainAll = {'m10','m50','m100','e10','e100','e1000'};
    main = mainAll(cellfun(@(v) any(T.(v)~=0), mainAll));
    if hasCh; main = [main {'chLVN'}]; end
    interAll = {'m10:e10','m10:e100','m10:e1000','m50:e10','m50:e100','m50:e1000', ...
                'm100:e10','m100:e100','m100:e1000'};
    inter = interAll(cellfun(@(v) any(prod_col(T,v)~=0), interAll));
    if isempty(inter)
        warning('bulk_mixed_models:noInter','%s: no combined (MxE) cells -- interaction not testable.', label);
        return;
    end

    fixedRed  = ['y ~ -1 + ' strjoin(main,' + ')];
    fixedFull = [fixedRed ' + ' strjoin(inter,' + ')];

    % random effect identifiable only with >=2 obs per animal
    obsPer = groupcounts_local(T.animal);
    useRand = max(obsPer) >= 2;

    try
        if useRand
            lmeFull = fitlme(T, [fixedFull ' + (1|animal)'], 'FitMethod','ML');
            lmeRed  = fitlme(T, [fixedRed  ' + (1|animal)'], 'FitMethod','ML');
            cmp = compare(lmeRed, lmeFull);   % theoretical LRT
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
        warning('bulk_mixed_models:fit','%s: fit failed (%s)', label, ME.message); return;
    end

    % collect per-cell synergy (interaction coefficients)
    isI = contains(CoefName, ':');
    idx = find(isI); pvals = PV(idx); pFDR = bh_fdr(pvals);
    for q = 1:numel(idx)
        [M,E] = parse_inter(CoefName{idx(q)});
        cells{end+1,1} = struct('metric',label,'source',source,'M',M,'E',E, ...
            'synergy',Est(idx(q)),'SE',SE(idx(q)),'p',pvals(q),'pFDR',pFDR(q)); %#ok<AGROW>
    end

    % ---- heatmap data: EVERY cell = that condition's relative change from
    %      baseline (%), with p vs baseline (0). Each cell's modeled mean is a
    %      contrast of the fitted coefficients:
    %        M alone (k)      = alpha_k
    %        E alone (l)      = beta_l
    %        combined (k,l)   = alpha_k + beta_l + gamma_kl
    %      p comes from coefTest on that contrast; FDR across all cells.
    mLevAll = [10 50 100]; eLevAll = [10 100 1000];
    has = @(nm) any(strcmp(CoefName, nm));
    mLev = mLevAll(arrayfun(@(L) has(sprintf('m%d',L)), mLevAll));
    eLev = eLevAll(arrayfun(@(L) has(sprintf('e%d',L)), eLevAll));
    nMl = numel(mLev); nEl = numel(eLev);
    Z = nan(nMl+1, nEl+1); P = nan(nMl+1, nEl+1);
    for j = 1:nEl
        [e,p] = contrast_sum(mdl, CoefName, Est, {sprintf('e%d',eLev(j))});
        Z(1,1+j) = 100*e; P(1,1+j) = p;
    end
    for i = 1:nMl
        [e,p] = contrast_sum(mdl, CoefName, Est, {sprintf('m%d',mLev(i))});
        Z(1+i,1) = 100*e; P(1+i,1) = p;
    end
    for i = 1:nMl
        for j = 1:nEl
            [e,p] = contrast_sum(mdl, CoefName, Est, ...
                {sprintf('m%d',mLev(i)), sprintf('e%d',eLev(j)), sprintf('m%d:e%d',mLev(i),eLev(j))});
            Z(1+i,1+j) = 100*e; P(1+i,1+j) = p;
        end
    end
    fin = isfinite(P); q = nan(size(P)); if any(fin(:)); q(fin) = bh_fdr(P(fin)); end
    Sg = isfinite(q) & q < 0.05;
    hm = struct('Mlev',mLev,'Elev',eLev,'Z',Z,'sig',Sg);
end

% ======================================================================
function [est, p] = contrast_sum(mdl, CoefName, Est, names)
% Estimate and p-value for the linear contrast that SUMS the named coefficients
% (tests sum == 0). Used to read each condition's fitted mean change from the
% fitted M/E/interaction coefficients. Works for fitlme and fitlm.
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
    if isempty(idx) && contains(n, ':')   % interaction stored in the other order
        pp = strsplit(n, ':'); idx = find(strcmp(CoefName, [pp{2} ':' pp{1}]), 1);
    end
end

% ======================================================================
function maybe_plot(hm, label, saveDir, formats)
% Render (and save) the heatmap for one metric: EVERY cell is that condition's
% relative change from baseline (%); '*' + bold border = change significantly
% different from baseline (FDR<0.05 across the cells).
    if isempty(hm); return; end
    interior = hm.Z(2:end, 2:end);
    if ~any(isfinite(interior(:))); return; end
    ttl = sprintf('%s : relative change from baseline (%%)', label);
    cbl = 'relative change from baseline (%)   (* FDR<0.05 vs baseline)';
    f = me_heatmap_render(hm.Z, hm.Mlev, hm.Elev, ttl, cbl, hm.sig);
    if ~isempty(saveDir) && exist('save_one_figure','file')==2
        save_one_figure(f, saveDir, ['mixedmodel_' regexprep(label,'[^\w]+','_')], formats);
    end
end

% ======================================================================
function T = build_nerve_table(normRows, metric, channel)
% Reuse the pipeline's OWN baseline-normalized scalars (out.norm from
% bulk_compile): one value per (animal, condition, channel), already paired
% baseline<->recovery by (animal,condition,channel), repeat trials pooled, and
% animal parsed by parse_stem -- exactly the numbers every box/violin and
% heatmap uses. y = scalar.(metric) = (mean recovery - mean baseline)/mean baseline.
    sel = find(strcmpi({normRows.label}, channel));
    y=[]; an={}; Mv=[]; Ev=[];
    for i = sel(:)'
        if ~isfield(normRows(i).scalar, metric); continue; end
        v = normRows(i).scalar.(metric); if ~isfinite(v); continue; end
        [M,E] = parse_me(normRows(i).condition);
        y(end+1,1)=v; an{end+1,1}=normRows(i).animal; Mv(end+1,1)=M; Ev(end+1,1)=E; %#ok<AGROW>
    end
    T = make_table(y, an, Mv, Ev);
end

% ======================================================================
function Tsys = build_systemic_tables()
    state = buildFileQueue(fullfile(pwd,'gemsplots_queue.mat'));
    assert(~isempty(state.files), 'no files queued for systemic metrics');
    specs = hrv_specs();
    seriesByFile = loadAllSeriesCached(state, specs, fullfile(pwd,'gemsplots_series_cache.mat'));
    nF = numel(state.files); nM = numel(specs);
    fmean = nan(nF,nM);
    for i=1:nF
        for k=1:nM
            sObj = seriesByFile{i,k};
            if ~isempty(sObj) && isfield(sObj,'y') && ~isempty(sObj.y); fmean(i,k)=mean(sObj.y,'omitnan'); end
        end
    end
    isBase = strcmpi(state.phase,'baseline'); isRec = strcmpi(state.phase,'recovery');
    Tsys = struct('label',{},'T',{});
    for k = 1:nM
        y=[]; an={}; Mv=[]; Ev=[];
        for r = find(isRec(:))'
            a = state.animal{r}; c = state.condition{r};
            % pair baseline<->recovery by (animal,condition), same convention as
            % compute_norm in plotSynergyHeatmaps / bulk_hrv_heatmaps
            cand = find(isBase(:)' & strcmpi(state.animal,a) & strcmpi(state.condition,c));
            if isempty(cand); continue; end
            bm = fmean(cand(1),k); rm = fmean(r,k);
            if ~isfinite(bm)||~isfinite(rm)||bm==0; continue; end
            [M,E]=parse_me(c);
            y(end+1,1)=(rm-bm)/bm; an{end+1,1}=a; Mv(end+1,1)=M; Ev(end+1,1)=E; %#ok<AGROW>
        end
        Tsys(end+1) = struct('label',specs(k).label,'T',make_table(y,an,Mv,Ev)); %#ok<AGROW>
    end
end

% ======================================================================
function T = make_table(y, an, Mv, Ev)
    % animal already parsed/canonicalised upstream (parse_stem -> lowercase letter)
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

% ---- HR/HRV metric specs (mirror of bulk_hrv_heatmaps) ----
function specs = hrv_specs()
    specs = [ ...
      ms('HR',             '_HRBR.mat',        'avgHeartRate','heartRateSeries'); ...
      ms('Breathing rate', '_HRBR.mat',        'avgBreathRate','breathRateSeries'); ...
      ms('HRV',            '_HRVMeasures.mat', 'hrv',   'hrv_series'); ...
      ms('RMSSD',          '_HRVMeasures.mat', 'rmssd', 'rmssd_series'); ...
      ms('pNN5',           '_HRVMeasures.mat', 'pnn5',  'pnn5_series'); ...
      ms('Sample entropy', '_HRVMeasures.mat', 'sampEn','sampEn_series'); ...
      ms('SD1',            '_HRVMeasures.mat', 'sd1',   'sd1_series'); ...
      ms('SD2',            '_HRVMeasures.mat', 'sd2',   'sd2_series') ];
end

function s = ms(label, suffix, field, seriesField)
    s = struct('label',label,'unitsIn','','unitsOut','','suffix',suffix, ...
               'field',field,'seriesField',seriesField,'timeField','metrics_t', ...
               'aggregator','auto','channel',NaN);
end
