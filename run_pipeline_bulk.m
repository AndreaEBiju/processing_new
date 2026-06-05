function out = run_pipeline_bulk(P, opts)
% RUN_PIPELINE_BULK  Cached batch driver for the cross-condition figures.
%
%   out = run_pipeline_bulk(P, opts)
%
% Workflow (per our spec):
%   1. bulk_conventions_ui  -> 4 file suffixes, RVN/LVN columns, heartbeat var
%      (typed once, applied to every file; persisted).
%   2. pick folder(s) -> bulk_scan_files pairs each neural file with its HRBR.
%   3. For each file: if it is already in the cache (same path + mod-time and
%      not forced) reuse it; otherwise load -> process_dataset (FULL per-file
%      diagnostic figures) -> step5f Fano slope -> save_all_figures into a
%      'figures' folder next to the file -> harvest -> write to the cache.
%      Editing a plot later reprocesses NOTHING (recompiles from cache).
%   4. defineGroupsUI -> groups.
%   5. bulk_compile baseline-normalizes every recovery sample as
%      (x - mean_baseline)/mean_baseline, then renders, per channel:
%      box+violin (rate/Vpp/FWHM/CV2), synergy heatmaps, and a Fano box.
%
% opts (optional): .cacheFile (default ./bulk_cache.mat), .forceRefresh
% (default false), .metricsBox {'rate','vpp','fwhm','cv2'}, .metricsSyn
% {'rate','excess','vpp','fwhm','cv2'}.

    if nargin < 1 || isempty(P);    P = pipeline_params(); end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts,'cacheFile');    opts.cacheFile = fullfile(pwd,'bulk_cache.mat'); end
    if ~isfield(opts,'forceRefresh'); opts.forceRefresh = false; end
    if ~isfield(opts,'metricsBox');   opts.metricsBox = {'rate','vpp','fwhm','cv2'}; end
    if ~isfield(opts,'metricsSyn');   opts.metricsSyn = {'rate','excess','vpp','fwhm','cv2'}; end
    if ~isfield(opts,'perFileFigures'); opts.perFileFigures = true; end   % false = fast metrics-only
    if ~isfield(opts,'saveFigsDir');  opts.saveFigsDir = ''; end          % save combined plots here (png/svg/fig)
    cf = opts.cacheFile;

    conv = bulk_conventions_ui(cf);
    if isempty(conv); fprintf('Cancelled at conventions.\n'); out = []; return; end

    files = bulk_queue_ui(conv, cf);   % table: Add folder -> recursive scan -> edit/delete -> Continue
    if isempty(files); fprintf('No files queued.\n'); out = []; return; end

    harvest = load_harvest(cf);
    loadOpts = struct('neuralCols',[conv.rvnCol conv.lvnCol], 'labels',{{'RVN','LVN'}}, 'rVar',conv.heartVar);
    fprintf('[cache] file: %s\n[cache] %d entr(y/ies) loaded.\n', cf, numel(harvest));

    % which files need (re)processing? central cache -> per-file .mat -> reprocess
    SCHEMA = veng_schema();
    mtimes = nan(numel(files),1); need = false(numel(files),1);
    for f = 1:numel(files)
        di = dir(files(f).neural); if ~isempty(di); mtimes(f) = di.datenum; end
        if opts.forceRefresh; need(f) = true; continue; end
        ci = find_cache(harvest, files(f).neural);
        if ~isempty(ci) && harvest(ci).mtime == mtimes(f) && harvest(ci).schema == SCHEMA
            need(f) = false; continue;                      % in-memory/central cache hit
        end
        pf = load_perfile(files(f).neural, mtimes(f), SCHEMA);   % <base>_vengmetrics.mat next to file
        if ~isempty(pf)
            harvest = put_cache(harvest, files(f).neural, mtimes(f), pf);
            need(f) = false; continue;
        end
        need(f) = true;
    end
    proc = find(need);
    fprintf('[bulk] %d cached (skipped), %d to process.\n', numel(files)-numel(proc), numel(proc));

    % parallel if available (and requested); cache writes stay serial afterwards
    usePar = (~isfield(opts,'parallel') || opts.parallel) && ~isempty(proc) && ...
             license('test','Distrib_Computing_Toolbox');
    if usePar
        try; if isempty(gcp('nocreate')); parpool; end
        catch ME; warning('bulk:pool','parpool failed (%s); serial.', ME.message); usePar = false; end
    end

    filesProc = files(proc); newRows = cell(numel(proc),1);
    makeFigs = opts.perFileFigures;
    if usePar
        parfor j = 1:numel(proc)
            newRows{j} = try_process(filesProc(j), P, loadOpts, makeFigs, j, numel(proc));
        end
    else
        for j = 1:numel(proc)
            newRows{j} = try_process(filesProc(j), P, loadOpts, makeFigs, j, numel(proc));
        end
    end
    for j = 1:numel(proc)
        if ~isempty(newRows{j})
            harvest = put_cache(harvest, files(proc(j)).neural, mtimes(proc(j)), newRows{j});
        end
    end
    save_harvest(cf, harvest);

    % assemble all rawRows from cache (cached + freshly processed)
    rawRows = repmat(empty_row(), 0, 1);
    for f = 1:numel(files)
        ci = find_cache(harvest, files(f).neural);
        if isempty(ci); continue; end
        rws = harvest(ci).rows;
        for k = 1:numel(rws)
            if ~isfield(rws(k).dist,'rate_t'); continue; end   % skip rows from a pre-timing cache
            rawRows(end+1) = rws(k); %#ok<AGROW>
        end
    end
    if isempty(rawRows); fprintf('Nothing harvested.\n'); out = []; return; end

    groups = [];
    if exist('defineGroupsUI','file')==2
        prior = bulk_cache_get(cf, 'groups'); if ~iscell(prior); prior = {}; end
        [groups, ok] = defineGroupsUI(unique({rawRows.condition}), prior);   % prefilled from cache
        if ~ok; groups = []; end
        if ok && ~isempty(groups); bulk_cache_set(cf, 'groups', groups); end % remember for next run
    end

    normRows = bulk_compile(rawRows);
    out = struct('raw', rawRows, 'norm', normRows, 'conv', conv, 'groups', {groups});
    try; save(fullfile(fileparts(cf),'BULK.mat'), 'out'); catch; end

    if ~isempty(groups)
        for ci2 = 1:2
            chL = loadOpts.labels{ci2};
            for m = 1:numel(opts.metricsBox)
                bulk_plot_boxviolin(normRows, groups, opts.metricsBox{m}, chL);  % full-duration
                bulk_plot_windowed(rawRows,  groups, opts.metricsBox{m}, chL);   % 1/2/5/10min/full
            end
            bulk_plot_fano_box(normRows, groups, chL);
            for m = 1:numel(opts.metricsSyn)
                bulk_plot_synergy(normRows, opts.metricsSyn{m}, chL);
            end
        end
        bulk_plot_windowed_total(rawRows, groups);   % combined (LVN+RVN) rate, windowed
        if ~isempty(opts.saveFigsDir)
            if ~exist(opts.saveFigsDir,'dir'); mkdir(opts.saveFigsDir); end
            save_all_figures(opts.saveFigsDir, 'bulk', {'png','svg','fig'});
            fprintf('[bulk] saved combined figures to %s\n', opts.saveFigsDir);
        end
        % all M x E mean-%-change heatmaps: nerve (rate/Vpp/FWHM/CV2 x RVN/LVN),
        % total (LVN+RVN) rate, and HR / breathing / HRV / RMSSD / pNN5 /
        % sample entropy / SD1 / SD2. These self-save via save_one_figure, so run
        % AFTER save_all_figures to avoid double-saving the same figures.
        try
            bulk_all_heatmaps(out, opts.saveFigsDir, {'png','svg','fig'});
        catch ME
            warning('bulk:heatmaps','heatmap generation failed (%s)', ME.message);
        end
    else
        fprintf('[bulk] groups not set. Render later, e.g.:\n');
        fprintf('   bulk_plot_boxviolin(out.norm, groups, ''rate'', ''RVN'');\n');
    end
end

% ======================================================================
function rows = try_process(F, P, loadOpts, makeFigs, j, n)
    fprintf('[bulk %d/%d] processing: %s  (%s / %s)\n', j, n, F.stem, F.condition, F.phase);
    try
        rows = process_and_harvest(F, P, loadOpts, makeFigs);
    catch ME
        warning('bulk:fail','  skipped %s (%s)', F.stem, ME.message); rows = [];
    end
end

function rows = process_and_harvest(F, P, loadOpts, makeFigs)
    D = bulk_load_one(F.neural, F.hrbr, loadOpts);
    if makeFigs
        prev = get(0,'DefaultFigureVisible'); set(0,'DefaultFigureVisible','off');
        cu = onCleanup(@() set(0,'DefaultFigureVisible',prev)); %#ok<NASGU>
        D  = process_dataset(D, P, true);        % FULL per-file diagnostic figures
        Rf = step5f_fano_slope(D, P, true);
        [nfDir, nfBase] = fileparts(F.neural);
        try
            save_all_figures(fullfile(nfDir,'figures'), nfBase, {'png','fig','svg'});
        catch ME
            warning('bulk:figsave','  figure save failed (%s)', ME.message);
        end
        close all force; %#ok<CLALL>
    else
        D  = process_dataset(D, P, false);       % fast: no per-file figures
        Rf = step5f_fano_slope(D, P, false);
    end
    rows = repmat(empty_row(), numel(D.neuralChannels), 1);
    for k = 1:numel(D.neuralChannels)
        rows(k) = harvest_channel(D, Rf, k, F.animal, F.condition, F.phase);
    end
    % persist per-file metrics next to the neural file, so later runs/plots read
    % them without recomputing -- even if the central cache is deleted.
    try
        [nfDir, nfBase] = fileparts(F.neural);
        rows_ = rows; srcMtime = file_mtime(F.neural); schema = veng_schema(); %#ok<NASGU>
        save(fullfile(nfDir, [nfBase '_vengmetrics.mat']), 'rows_', 'srcMtime', 'schema');
    catch ME
        warning('bulk:perfile','per-file metrics save failed for %s (%s)', F.stem, ME.message);
    end
end

% ---- per-channel harvest (raw distributions + means + fano slope) ----
function r = empty_row()
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'rate_t',[],'vpp',[],'fwhm',[],'spk_t',[],'cv2',[],'cv2_t',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN);
end

function r = harvest_channel(D, Rf, k, animal, condition, phase)
    r = empty_row();
    r.animal = animal; r.condition = condition; r.phase = phase; r.label = D.channelLabels{k};
    fs = D.fs; valid = D.validMask(:,k);
    if isfield(D,'metrics') && numel(D.metrics)>=k && ~isempty(D.metrics(k).fr_hz)
        r.dist.rate = D.metrics(k).fr_hz(:);
        if isfield(D.metrics,'fr_t'); r.dist.rate_t = D.metrics(k).fr_t(:); end
    end
    if isfield(D.spikes,'Vpp_uv');   r.dist.vpp  = D.spikes(k).Vpp_uv(:);   end
    if isfield(D.spikes,'width_ms'); r.dist.fwhm = D.spikes(k).width_ms(:); end
    if isfield(D.spikes,'alignedTimes') && ~isempty(D.spikes(k).alignedTimes)
        r.dist.spk_t = D.spikes(k).alignedTimes(:);          % spike times (s)
        cen = round(D.spikes(k).alignedTimes(:)*fs)+1;
    elseif isfield(D.spikes,'times'); cen = round(D.spikes(k).times(:)*fs)+1;
    else; cen = []; end
    [r.dist.cv2, r.dist.cv2_t] = cv2_clean(cen, valid, fs);
    if isfield(D,'metrics') && numel(D.metrics)>=k;   r.mean.rate   = D.metrics(k).meanRate_hz; end
    if isfield(D,'envelope') && numel(D.envelope)>=k; r.mean.excess = D.envelope(k).meanExcess_uv; end
    r.mean.vpp  = median(r.dist.vpp,'omitnan');
    r.mean.fwhm = median(r.dist.fwhm,'omitnan');
    r.mean.cv2  = mean(r.dist.cv2,'omitnan');
    if ~isempty(Rf) && numel(Rf)>=k && isfield(Rf,'slopeNorm'); r.fanoSlope = Rf(k).slopeNorm; end
end

function [cv2, cv2t] = cv2_clean(cen, valid, fs)
% per-pair CV2 from gap-clean consecutive ISIs + the time (s) of each pair's
% middle spike (used for within-recovery windowing).
    cv2 = []; cv2t = [];
    if nargin < 3 || isempty(fs); fs = 1; end
    cen = sort(cen(:)); cen = cen(cen>=1 & cen<=numel(valid));
    if numel(cen) < 3; return; end
    cvv = [0; cumsum(double(valid(:)))];
    span = diff(cen); vbet = cvv(cen(2:end)) - cvv(cen(1:end-1));
    isClean = (vbet == span); isi = span;
    for i = 1:numel(isi)-1
        if isClean(i) && isClean(i+1)
            a = isi(i); b = isi(i+1);
            if (a+b)>0
                cv2(end+1,1)  = 2*abs(b-a)/(a+b);    %#ok<AGROW>
                cv2t(end+1,1) = (cen(i+1)-1)/fs;     %#ok<AGROW>
            end
        end
    end
end

% ---- cache helpers ----
function s = veng_schema(); s = 2; end   % bump when the harvested row schema changes

function harvest = load_harvest(cf)
    harvest = bulk_cache_get(cf, 'harvest');
    if isempty(harvest) || ~isstruct(harvest)
        harvest = struct('neuralPath',{},'mtime',{},'rows',{},'schema',{});
        return;
    end
    if ~isfield(harvest,'schema'); [harvest.schema] = deal(-1); end   % old entries -> stale
end
function ci = find_cache(harvest, neuralPath)
    ci = [];
    if isempty(harvest); return; end
    ci = find(strcmp({harvest.neuralPath}, neuralPath), 1);
end
function harvest = put_cache(harvest, neuralPath, mt, rows)
    e = struct('neuralPath',neuralPath,'mtime',mt,'rows',rows,'schema',veng_schema());
    ci = find_cache(harvest, neuralPath);
    if isempty(ci); harvest(end+1) = e; else; harvest(ci) = e; end
end

function rows = load_perfile(neuralPath, mt, schema)
% read a per-file <base>_vengmetrics.mat if it is current (matching mtime+schema)
    rows = [];
    [d, b] = fileparts(neuralPath);
    pf = fullfile(d, [b '_vengmetrics.mat']);
    if ~exist(pf,'file'); return; end
    try
        S = load(pf);
        if isfield(S,'schema') && isequal(S.schema,schema) && isfield(S,'srcMtime') ...
                && isequal(S.srcMtime,mt) && isfield(S,'rows_')
            rows = S.rows_;
        end
    catch
    end
end

function mt = file_mtime(p)
    di = dir(p); if isempty(di); mt = NaN; else; mt = di.datenum; end
end
function save_harvest(cf, harvest)
    bulk_cache_set(cf, 'harvest', harvest);
end
