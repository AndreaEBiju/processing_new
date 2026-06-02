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
    cf = opts.cacheFile;

    conv = bulk_conventions_ui(cf);
    if isempty(conv); fprintf('Cancelled at conventions.\n'); out = []; return; end

    files = bulk_queue_ui(conv, cf);   % table: Add folder -> recursive scan -> edit/delete -> Continue
    if isempty(files); fprintf('No files queued.\n'); out = []; return; end

    harvest = load_harvest(cf);
    loadOpts = struct('neuralCols',[conv.rvnCol conv.lvnCol], 'labels',{{'RVN','LVN'}}, 'rVar',conv.heartVar);
    fprintf('[cache] file: %s\n[cache] %d entr(y/ies) loaded.\n', cf, numel(harvest));

    % which files need (re)processing?
    mtimes = nan(numel(files),1); need = false(numel(files),1);
    for f = 1:numel(files)
        di = dir(files(f).neural); if ~isempty(di); mtimes(f) = di.datenum; end
        ci = find_cache(harvest, files(f).neural);
        if isempty(ci)
            need(f) = true;
        elseif opts.forceRefresh
            need(f) = true;
        elseif harvest(ci).mtime ~= mtimes(f)
            need(f) = true;
            fprintf('[cache] mtime changed -> reprocess: %s\n', files(f).stem);
        else
            need(f) = false;
        end
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
        for k = 1:numel(rws); rawRows(end+1) = rws(k); end %#ok<AGROW>
    end
    if isempty(rawRows); fprintf('Nothing harvested.\n'); out = []; return; end

    groups = [];
    if exist('defineGroupsUI','file')==2
        [groups, ok] = defineGroupsUI(unique({rawRows.condition}), {});
        if ~ok; groups = []; end
    end

    normRows = bulk_compile(rawRows);
    out = struct('raw', rawRows, 'norm', normRows, 'conv', conv, 'groups', {groups});
    try; save(fullfile(fileparts(cf),'BULK.mat'), 'out'); catch; end

    if ~isempty(groups)
        for ci2 = 1:2
            chL = loadOpts.labels{ci2};
            for m = 1:numel(opts.metricsBox)
                bulk_plot_boxviolin(normRows, groups, opts.metricsBox{m}, chL);
            end
            bulk_plot_fano_box(normRows, groups, chL);
            for m = 1:numel(opts.metricsSyn)
                bulk_plot_synergy(normRows, opts.metricsSyn{m}, chL);
            end
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
end

% ---- per-channel harvest (raw distributions + means + fano slope) ----
function r = empty_row()
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'vpp',[],'fwhm',[],'cv2',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN);
end

function r = harvest_channel(D, Rf, k, animal, condition, phase)
    r = empty_row();
    r.animal = animal; r.condition = condition; r.phase = phase; r.label = D.channelLabels{k};
    fs = D.fs; valid = D.validMask(:,k);
    if isfield(D,'metrics') && numel(D.metrics)>=k && ~isempty(D.metrics(k).fr_hz)
        r.dist.rate = D.metrics(k).fr_hz(:);
    end
    if isfield(D.spikes,'Vpp_uv');   r.dist.vpp  = D.spikes(k).Vpp_uv(:);   end
    if isfield(D.spikes,'width_ms'); r.dist.fwhm = D.spikes(k).width_ms(:); end
    if isfield(D.spikes,'alignedTimes') && ~isempty(D.spikes(k).alignedTimes)
        cen = round(D.spikes(k).alignedTimes(:)*fs)+1;
    elseif isfield(D.spikes,'times'); cen = round(D.spikes(k).times(:)*fs)+1;
    else; cen = []; end
    r.dist.cv2 = cv2_clean(cen, valid);
    if isfield(D,'metrics') && numel(D.metrics)>=k;   r.mean.rate   = D.metrics(k).meanRate_hz; end
    if isfield(D,'envelope') && numel(D.envelope)>=k; r.mean.excess = D.envelope(k).meanExcess_uv; end
    r.mean.vpp  = median(r.dist.vpp,'omitnan');
    r.mean.fwhm = median(r.dist.fwhm,'omitnan');
    r.mean.cv2  = mean(r.dist.cv2,'omitnan');
    if ~isempty(Rf) && numel(Rf)>=k && isfield(Rf,'slopeNorm'); r.fanoSlope = Rf(k).slopeNorm; end
end

function cv2 = cv2_clean(cen, valid)
    cv2 = [];
    cen = sort(cen(:)); cen = cen(cen>=1 & cen<=numel(valid));
    if numel(cen) < 3; return; end
    cvv = [0; cumsum(double(valid(:)))];
    span = diff(cen); vbet = cvv(cen(2:end)) - cvv(cen(1:end-1));
    isClean = (vbet == span); isi = span;
    for i = 1:numel(isi)-1
        if isClean(i) && isClean(i+1)
            a = isi(i); b = isi(i+1);
            if (a+b)>0; cv2(end+1,1) = 2*abs(b-a)/(a+b); end %#ok<AGROW>
        end
    end
end

% ---- cache helpers ----
function harvest = load_harvest(cf)
    harvest = bulk_cache_get(cf, 'harvest');
    if isempty(harvest) || ~isstruct(harvest)
        harvest = struct('neuralPath',{},'mtime',{},'rows',{});
    end
end
function ci = find_cache(harvest, neuralPath)
    ci = [];
    if isempty(harvest); return; end
    ci = find(strcmp({harvest.neuralPath}, neuralPath), 1);
end
function harvest = put_cache(harvest, neuralPath, mt, rows)
    e = struct('neuralPath',neuralPath,'mtime',mt,'rows',rows);
    ci = find_cache(harvest, neuralPath);
    if isempty(ci); harvest(end+1) = e; else; harvest(ci) = e; end
end
function save_harvest(cf, harvest)
    bulk_cache_set(cf, 'harvest', harvest);
end
