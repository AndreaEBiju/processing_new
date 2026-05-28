function [values, hitCount, loadCount] = loadAllMetricsCached(state, metricSpecs, cacheFile, opts)
%LOADALLMETRICSCACHED  Pre-load every (source, metric) scalar with a disk
% cache + optional parallel loading. First run is parallelised across
% files (one task per source); subsequent runs read from the local cache
% and skip parallelisation entirely.
%
%   [values, hitCount, loadCount] = loadAllMetricsCached(state, metricSpecs, cacheFile)
%   [values, hitCount, loadCount] = loadAllMetricsCached(state, metricSpecs, cacheFile, opts)
%
%   opts (struct, all optional)
%     .useParallel  true (default) — use parfeval over files when the
%                   Parallel Computing Toolbox is available
%     .nWorkers     numeric — desired pool size; default = the existing
%                   pool if running, else MATLAB's default
%
%   Cache invalidation by source-file mtime — re-running batch_process
%   on a file invalidates its cache entries. Delete the cache file to
%   force a full reload.

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'useParallel'), opts.useParallel = true; end
    if ~isfield(opts,'nWorkers'),    opts.nWorkers    = []; end

    nFiles   = numel(state.files);
    nMetrics = numel(metricSpecs);
    values   = nan(nFiles, nMetrics);

    % ----- restore cache -----
    cacheMap = containers.Map('KeyType','char','ValueType','any');
    if ~isempty(cacheFile) && exist(cacheFile, 'file')
        try
            c = load(cacheFile);
            if isfield(c, 'cacheMap')
                cacheMap = c.cacheMap;
            end
        catch ME
            warning('loadAllMetricsCached:restore', ...
                'Cache restore failed (%s); starting fresh.', ME.message);
        end
    end

    % per-metric stable keys
    specKeys = cell(nMetrics,1);
    for k = 1:nMetrics
        specKeys{k} = makeSpecKey(metricSpecs(k));
    end

    % ----- pass 1: cache hits + plan loads -----
    needsLoad  = false(nFiles, nMetrics);
    fileMtimes = zeros(nFiles, 1);
    hitCount = 0;

    for i = 1:nFiles
        info = dir(state.files{i});
        if ~isempty(info), fileMtimes(i) = info(1).datenum; end
        for k = 1:nMetrics
            ck = [state.files{i} '||' specKeys{k}];
            if isKey(cacheMap, ck) && cacheMap(ck).mtime == fileMtimes(i)
                values(i,k) = cacheMap(ck).value;
                hitCount = hitCount + 1;
            else
                needsLoad(i,k) = true;
            end
        end
    end
    loadCount = sum(needsLoad(:));

    if loadCount == 0
        return;
    end

    % ----- group by file (one task per file) -----
    fileIdxToLoad = find(any(needsLoad, 2));
    nF = numel(fileIdxToLoad);
    perFileSpecIdx = cell(nF,1);
    for j = 1:nF
        perFileSpecIdx{j} = find(needsLoad(fileIdxToLoad(j), :));
    end

    % ----- progress UI -----
    wb = waitbar(0, sprintf('Loading %d file(s)  (%d values)...', nF, loadCount), ...
        'Name','plotAveragedMetrics: loading');
    try, set(wb,'WindowState','normal'); catch, end %#ok<CTCH>
    setWaitbarInterpNone(wb);                 % file paths contain \ and _
    cleanup = onCleanup(@() safeCloseWb(wb)); %#ok<NASGU>

    % ----- start parpool if requested + available -----
    canParallel = opts.useParallel && hasParallelToolbox();
    if canParallel && isempty(gcp('nocreate'))
        if ishandle(wb)
            waitbar(0, wb, 'Starting parallel pool...');
        end
        try
            if isempty(opts.nWorkers)
                parpool('local');
            else
                parpool('local', opts.nWorkers);
            end
        catch ME
            warning('loadAllMetricsCached:parpool', ...
                'Could not start parpool: %s. Falling back to serial.', ME.message);
            canParallel = false;
        end
    end

    % ----- dispatch loads -----
    nDone = 0;
    if canParallel
        % one parfeval per file; each worker reads that file's metrics
        futures = parallel.FevalFuture.empty(0,1);
        for j = 1:nF
            i = fileIdxToLoad(j);
            specsForFile = metricSpecs(perFileSpecIdx{j});
            futures(j,1) = parfeval(@loadFileMetrics, 1, ...
                state.files{i}, specsForFile);
        end

        for j = 1:nF
            try
                [completedIdx, vals] = fetchNext(futures);
            catch ME
                warning('loadAllMetricsCached:parfeval', ...
                    'parfeval failed at task %d: %s', j, ME.message);
                cancel(futures);
                break;
            end
            i = fileIdxToLoad(completedIdx);
            kList = perFileSpecIdx{completedIdx};
            for kk = 1:numel(kList)
                k = kList(kk);
                values(i,k) = vals(kk);
                ck = [state.files{i} '||' specKeys{k}];
                cacheMap(ck) = struct('mtime', fileMtimes(i), 'value', vals(kk));
            end
            nDone = nDone + 1;
            updateBar(wb, nDone, nF, state.files{i});
        end
    else
        % serial fallback
        for j = 1:nF
            i = fileIdxToLoad(j);
            specsForFile = metricSpecs(perFileSpecIdx{j});
            vals = loadFileMetrics(state.files{i}, specsForFile);
            kList = perFileSpecIdx{j};
            for kk = 1:numel(kList)
                k = kList(kk);
                values(i,k) = vals(kk);
                ck = [state.files{i} '||' specKeys{k}];
                cacheMap(ck) = struct('mtime', fileMtimes(i), 'value', vals(kk));
            end
            nDone = nDone + 1;
            updateBar(wb, nDone, nF, state.files{i});
        end
    end

    % ----- persist cache -----
    if ~isempty(cacheFile)
        try
            save(cacheFile, 'cacheMap', '-v7.3');
        catch ME
            warning('loadAllMetricsCached:save', ...
                'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
function vals = loadFileMetrics(file, specs)
% Run on a worker (parfeval) or on main thread (serial fallback). Loads
% every requested metric from a single source file's result .mat files.
    nS = numel(specs);
    vals = nan(nS,1);
    for s = 1:nS
        vals(s) = loadMetric(file, specs(s));
    end
end

function k = makeSpecKey(spec)
    if isfield(spec,'channel') && isnumeric(spec.channel) && ~isnan(spec.channel)
        ch = num2str(spec.channel);
    else
        ch = '-';
    end
    k = sprintf('%s|%s|%s|%s', ...
        char(spec.suffix), char(spec.field), ...
        char(spec.aggregator), ch);
end

function tf = hasParallelToolbox()
    tf = ~isempty(ver('parallel'));
end

function safeCloseWb(wb)
    try
        if ishandle(wb), close(wb); end
    catch
    end
end

function setWaitbarInterpNone(wb)
% Force the waitbar's message text to render literally — file paths
% contain \ and _ which the LaTeX interpreter (pubfig_setup default)
% otherwise parses as escape sequences / subscripts.
    try
        ax  = findall(wb, 'Type','axes');
        txt = findall(ax, 'Type','text');
        if ~isempty(txt), set(txt, 'Interpreter','none'); end
    catch
    end
end

function updateBar(wb, nDone, nTotal, currentFile)
    if ishandle(wb)
        waitbar(nDone/nTotal, wb, sprintf( ...
            '[%d / %d files]   %s', nDone, nTotal, shortName(currentFile)));
    end
end

function s = shortName(fp)
    [d, b, e] = fileparts(char(fp));
    [~, parent] = fileparts(d);
    s = [parent filesep b e];
    if numel(s) > 70, s = ['...' s(end-66:end)]; end
end
