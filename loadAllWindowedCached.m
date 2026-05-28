function [baselineByFile, windowByFile, hitCount, loadCount] = ...
        loadAllWindowedCached(state, metricSpecs, winSecs, cacheFile, opts)
%LOADALLWINDOWEDCACHED  Time-series-aware version of loadAllMetricsCached.
%
% For every queued source file × every metric, compute
%   baselineByFile(i,k)   = mean of the FULL time series
%   windowByFile(i,k,w)   = mean of the time series over t in [t0, t0+W]
%                            for each window W in winSecs (Inf = full)
%
% Results are cached on disk (cacheFile) per (file, spec, windowKey) so
% subsequent runs skip the load. First run uses parfeval (one task per
% file) when the Parallel Computing Toolbox is available.
%
%   [bl, win, hits, loads] = loadAllWindowedCached(state, metricSpecs, winSecs, cacheFile)
%   ...                    = loadAllWindowedCached(..., opts)
%
% opts.useParallel  default true
% opts.nWorkers     default = MATLAB's local default pool size

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'useParallel'), opts.useParallel = true; end
    if ~isfield(opts,'nWorkers'),    opts.nWorkers    = []; end

    nFiles   = numel(state.files);
    nMetrics = numel(metricSpecs);
    nW       = numel(winSecs);
    baselineByFile = nan(nFiles, nMetrics);
    windowByFile   = nan(nFiles, nMetrics, nW);

    % ----- restore cache -----
    cacheMap = containers.Map('KeyType','char','ValueType','any');
    if ~isempty(cacheFile) && exist(cacheFile,'file')
        try
            c = load(cacheFile);
            if isfield(c, 'cacheMap'), cacheMap = c.cacheMap; end
        catch ME
            warning('loadAllWindowedCached:restore', ...
                'Cache restore failed (%s); starting fresh.', ME.message);
        end
    end

    specKeys = cell(nMetrics,1);
    for k = 1:nMetrics
        specKeys{k} = makeSpecKey(metricSpecs(k));
    end
    winKeys = cell(nW,1);
    for w = 1:nW
        if isinf(winSecs(w)), winKeys{w} = 'win_full';
        else,                 winKeys{w} = sprintf('win_%g', winSecs(w));
        end
    end

    % ----- pass 1: cache hits + plan loads -----
    hitCount = 0;
    needsLoad  = false(nFiles, 1);
    fileMtimes = zeros(nFiles, 1);
    for i = 1:nFiles
        info = dir(state.files{i});
        if ~isempty(info), fileMtimes(i) = info(1).datenum; end
        for k = 1:nMetrics
            ck_bl = [state.files{i} '||' specKeys{k} '||baseline'];
            if isKey(cacheMap, ck_bl) && cacheMap(ck_bl).mtime == fileMtimes(i) ...
                    && ~isnan(cacheMap(ck_bl).value)
                baselineByFile(i,k) = cacheMap(ck_bl).value;
                hitCount = hitCount + 1;
            else
                needsLoad(i) = true;
            end
            for w = 1:nW
                ck_w = [state.files{i} '||' specKeys{k} '||' winKeys{w}];
                if isKey(cacheMap, ck_w) && cacheMap(ck_w).mtime == fileMtimes(i) ...
                        && ~isnan(cacheMap(ck_w).value)
                    windowByFile(i,k,w) = cacheMap(ck_w).value;
                    hitCount = hitCount + 1;
                else
                    needsLoad(i) = true;
                end
            end
        end
    end

    filesToLoad = find(needsLoad);
    nLoad = numel(filesToLoad);
    loadCount = 0;
    if nLoad == 0, return; end

    % ----- progress UI -----
    wb = waitbar(0, sprintf('Loading time series from %d file(s)...', nLoad), ...
        'Name','plotWindowedMetrics: loading');
    try, set(wb,'WindowState','normal'); catch, end %#ok<CTCH>
    try
        ax_ = findall(wb,'Type','axes');
        set(findall(ax_,'Type','text'),'Interpreter','none');
    catch, end %#ok<CTCH>
    cleanup = onCleanup(@() safeCloseWb(wb)); %#ok<NASGU>

    % ----- start parpool if available -----
    canParallel = opts.useParallel && hasParallelToolbox();
    if canParallel && isempty(gcp('nocreate'))
        if ishandle(wb), waitbar(0, wb, 'Starting parallel pool...'); end
        try
            if isempty(opts.nWorkers), parpool('local');
            else,                      parpool('local', opts.nWorkers); end
        catch ME
            warning('loadAllWindowedCached:parpool', ...
                'Could not start parpool: %s. Falling back to serial.', ME.message);
            canParallel = false;
        end
    end

    % ----- dispatch -----
    nDone = 0;
    if canParallel
        futures = parallel.FevalFuture.empty(0,1);
        for j = 1:nLoad
            i = filesToLoad(j);
            futures(j,1) = parfeval(@computeFileValues, 2, ...
                state.files{i}, metricSpecs, winSecs);
        end
        for j = 1:nLoad
            try
                [completedIdx, blOne, winOne] = fetchNext(futures);
            catch ME
                warning('loadAllWindowedCached:parfeval', ...
                    'parfeval failed at task %d: %s', j, ME.message);
                cancel(futures);
                break;
            end
            i = filesToLoad(completedIdx);
            for k = 1:nMetrics
                baselineByFile(i,k) = blOne(k);
                if ~isnan(blOne(k))
                    cacheMap([state.files{i} '||' specKeys{k} '||baseline']) = ...
                        struct('mtime',fileMtimes(i),'value',blOne(k));
                    loadCount = loadCount + 1;
                end
                for w = 1:nW
                    windowByFile(i,k,w) = winOne(k,w);
                    if ~isnan(winOne(k,w))
                        cacheMap([state.files{i} '||' specKeys{k} '||' winKeys{w}]) = ...
                            struct('mtime',fileMtimes(i),'value',winOne(k,w));
                        loadCount = loadCount + 1;
                    end
                end
            end
            nDone = nDone + 1;
            if ishandle(wb)
                [~, b, e] = fileparts(state.files{i});
                waitbar(nDone/nLoad, wb, sprintf('[%d / %d] %s', nDone, nLoad, [b e]));
            end
        end
    else
        for j = 1:nLoad
            i = filesToLoad(j);
            [blOne, winOne] = computeFileValues(state.files{i}, metricSpecs, winSecs);
            for k = 1:nMetrics
                baselineByFile(i,k) = blOne(k);
                if ~isnan(blOne(k))
                    cacheMap([state.files{i} '||' specKeys{k} '||baseline']) = ...
                        struct('mtime',fileMtimes(i),'value',blOne(k));
                    loadCount = loadCount + 1;
                end
                for w = 1:nW
                    windowByFile(i,k,w) = winOne(k,w);
                    if ~isnan(winOne(k,w))
                        cacheMap([state.files{i} '||' specKeys{k} '||' winKeys{w}]) = ...
                            struct('mtime',fileMtimes(i),'value',winOne(k,w));
                        loadCount = loadCount + 1;
                    end
                end
            end
            nDone = nDone + 1;
            if ishandle(wb)
                [~, b, e] = fileparts(state.files{i});
                waitbar(nDone/nLoad, wb, sprintf('[%d / %d] %s', nDone, nLoad, [b e]));
            end
        end
    end

    % ----- persist cache -----
    if ~isempty(cacheFile)
        try
            save(cacheFile, 'cacheMap', '-v7.3');
        catch ME
            warning('loadAllWindowedCached:save', ...
                'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
function [baseline, windows] = computeFileValues(srcFile, metricSpecs, winSecs)
% Runs on a parfeval worker (or main thread in serial fallback). For one
% source file, loads each metric's time series and computes:
%   baseline(k)    : mean of the whole series (omitnan)
%   windows(k,w)   : mean of samples with t-t0 <= winSecs(w)
% Returns NaN for any (k) whose series cannot be loaded.
    nMetrics = numel(metricSpecs);
    nW = numel(winSecs);
    baseline = nan(nMetrics, 1);
    windows  = nan(nMetrics, nW);
    for k = 1:nMetrics
        [y, t] = loadMetricSeries(srcFile, metricSpecs(k));
        if isempty(y) || isempty(t), continue; end
        baseline(k) = mean(y, 'omitnan');
        t0 = t(1);
        for w = 1:nW
            W = winSecs(w);
            if isinf(W), sel = true(size(t));
            else,        sel = (t - t0) <= W;
            end
            if any(sel)
                windows(k, w) = mean(y(sel), 'omitnan');
            end
        end
    end
end


function k = makeSpecKey(spec)
    if isfield(spec,'channel') && isnumeric(spec.channel) && ~isnan(spec.channel)
        ch = num2str(spec.channel);
    else
        ch = '-';
    end
    k = sprintf('%s|%s|%s|%s', ...
        char(spec.suffix), char(spec.seriesField), ...
        char(spec.timeField), ch);
end

function tf = hasParallelToolbox()
    tf = ~isempty(ver('parallel'));
end

function safeCloseWb(wb)
    try, if ishandle(wb), close(wb); end, catch, end
end
