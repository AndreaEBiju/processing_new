function bulk_hrv_heatmaps(saveDir, formats)
% BULK_HRV_HEATMAPS  M x E mean-%-change heatmaps for the HR / HRV metrics,
% reusing YOUR loaders so field/series names match your pipeline exactly.
%
%   bulk_hrv_heatmaps(saveDir)
%
% Metrics: HR, breathing rate, HRV, RMSSD, pNN5, sample entropy, SD1, SD2
% (loaded from _HRBR.mat / _HRVMeasures.mat via loadAllSeriesCached, same metric
% specs as your plotSynergyHeatmaps). Each cell =
%       mean over animals of (mean(rec series) - mean(base series)) / mean(base series) * 100
% paired by (animal, condition). Marginals: E-alone top row, M-alone left column.
%
% Requires processing_new on the path (buildFileQueue, loadAllSeriesCached) and
% the data accessible. saveDir optional (png/svg/fig).

    if nargin < 1; saveDir = ''; end
    if nargin < 2 || isempty(formats); formats = {'png','svg','fig'}; end
    for fn = {'buildFileQueue','loadAllSeriesCached'}
        assert(exist(fn{1},'file')==2, ...
            'bulk_hrv_heatmaps: %s not on path -- addpath your processing_new folder.', fn{1});
    end

    qCache = fullfile(pwd,'gemsplots_queue.mat');
    sCache = fullfile(pwd,'gemsplots_series_cache.mat');

    state = buildFileQueue(qCache);
    if isempty(state.files); fprintf('No files queued.\n'); return; end

    specs = hrv_specs();
    [seriesByFile, ~, ~] = loadAllSeriesCached(state, specs, sCache);
    normByCond = compute_norm(state, specs, seriesByFile);   % nMetrics x nConds

    uconds = unique_stable(state.condition);
    [Mlev, Elev] = parse_ME(uconds);
    eL = sort(unique(Elev(Elev>0 & Mlev==0 & ~isnan(Elev))));
    mL = sort(unique(Mlev(Mlev>0 & Elev==0 & ~isnan(Mlev))));
    if isempty(eL) || isempty(mL)
        warning('bulk_hrv_heatmaps:levels','Need M-alone and E-alone conditions in the queue.'); return;
    end
    nM = numel(mL); nE = numel(eL);
    idxFor = @(cs) find(strcmpi(uconds, cs), 1);

    for k = 1:numel(specs)
        Z = nan(nM+1, nE+1);
        for j = 1:nE
            ix = idxFor(sprintf('E%d',eL(j))); if ~isempty(ix); Z(1,1+j) = 100*normByCond(k,ix); end
        end
        for i = 1:nM
            ix = idxFor(sprintf('M%d',mL(i))); if ~isempty(ix); Z(1+i,1) = 100*normByCond(k,ix); end
        end
        for i = 1:nM
            for j = 1:nE
                ix = idxFor(sprintf('M%dE%d',mL(i),eL(j)));
                if ~isempty(ix); Z(1+i,1+j) = 100*normByCond(k,ix); end
            end
        end
        nm = specs(k).label;
        f = me_heatmap_render(Z, mL, eL, sprintf('%s : mean %% change from baseline', nm), ...
            'mean % change from baseline');
        if ~isempty(saveDir) && exist('save_one_figure','file')==2
            save_one_figure(f, saveDir, sprintf('heatmap_%s', regexprep(nm,'\s+','_')), formats);
        end
    end
end

% ======================================================================
function normByCond = compute_norm(state, specs, seriesByFile)
% (rec_mean - base_mean)/base_mean per (animal, condition), averaged over animals
% (mirrors computeNormalisedByCondition in your plotSynergyHeatmaps).
    nFiles = numel(state.files); nMet = numel(specs);
    conds = state.condition; phases = state.phase; animals = state.animal;
    fileMean = nan(nFiles, nMet);
    for i = 1:nFiles
        for k = 1:nMet
            s = seriesByFile{i,k};
            if isempty(s) || ~isfield(s,'y') || isempty(s.y); continue; end
            fileMean(i,k) = mean(s.y, 'omitnan');
        end
    end
    uconds = unique_stable(conds);
    normByCond = nan(nMet, numel(uconds));
    for c = 1:numel(uconds)
        cond = uconds{c};
        baseRows = find(strcmpi(conds,cond) & strcmpi(phases,'baseline'));
        recRows  = find(strcmpi(conds,cond) & strcmpi(phases,'recovery'));
        if isempty(baseRows) || isempty(recRows); continue; end
        for k = 1:nMet
            vals = [];
            for rr = recRows(:)'
                iBase = baseRows(find(strcmpi(animals(baseRows), animals{rr}), 1));
                if isempty(iBase); continue; end
                bm = fileMean(iBase,k); rm = fileMean(rr,k);
                if isnan(bm) || isnan(rm) || bm==0; continue; end
                vals(end+1) = (rm - bm)/bm; %#ok<AGROW>
            end
            if ~isempty(vals); normByCond(k,c) = mean(vals); end
        end
    end
end

function specs = hrv_specs()
    specs = [ ...
      ms('HR',             'bpm','bpm','_HRBR.mat',        'avgHeartRate','heartRateSeries','metrics_t'); ...
      ms('Breathing rate', 'bpm','bpm','_HRBR.mat',        'avgBreathRate','breathRateSeries','metrics_t'); ...
      ms('HRV',            's',  'ms', '_HRVMeasures.mat', 'hrv',   'hrv_series',   'metrics_t'); ...
      ms('RMSSD',          's',  'ms', '_HRVMeasures.mat', 'rmssd', 'rmssd_series', 'metrics_t'); ...
      ms('pNN5',           '%',  '%',  '_HRVMeasures.mat', 'pnn5',  'pnn5_series',  'metrics_t'); ...
      ms('Sample entropy', '',   '',   '_HRVMeasures.mat', 'sampEn','sampEn_series','metrics_t'); ...
      ms('SD1',            's',  'ms', '_HRVMeasures.mat', 'sd1',   'sd1_series',   'metrics_t'); ...
      ms('SD2',            's',  'ms', '_HRVMeasures.mat', 'sd2',   'sd2_series',   'metrics_t') ];
end

function s = ms(label, unitsIn, unitsOut, suffix, field, seriesField, timeField)
    s = struct('label',label,'unitsIn',unitsIn,'unitsOut',unitsOut,'suffix',suffix, ...
               'field',field,'seriesField',seriesField,'timeField',timeField, ...
               'aggregator','auto','channel',NaN);
end

function [Mlev, Elev] = parse_ME(conds)
    n = numel(conds); Mlev = nan(n,1); Elev = nan(n,1);
    for k = 1:n
        c = upper(strtrim(char(conds{k})));
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

function u = unique_stable(c)
    if isempty(c); u = {}; return; end
    if ~iscell(c); c = cellstr(c); end
    [~, ia] = unique(c, 'stable'); u = c(sort(ia));
end
