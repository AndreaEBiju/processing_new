function [values, hitCount, loadCount] = loadAllMetricsCached(state, metricSpecs, cacheFile)
%LOADALLMETRICSCACHED  Pre-load every (source, metric) scalar with a disk
% cache so subsequent runs are instant — essential when source files live
% on a virtual filesystem like Google Drive where every load() is slow.
%
%   [values, hitCount, loadCount] = loadAllMetricsCached(state, metricSpecs, cacheFile)
%
%   state.files    cell of source-file paths
%   metricSpecs    struct array (from defineMetricsUI / defaultMetricSpecs)
%   cacheFile      path to a small .mat file (kept LOCAL, e.g. next to
%                  the queue cache on the SSD; never put it on Google
%                  Drive — defeats the purpose).
%
%   Returns
%     values     nFiles x nMetrics matrix of scalar metric values (NaN
%                where the source or its output file is missing).
%     hitCount   number of values read from cache (no disk I/O on source)
%     loadCount  number that had to be loaded fresh from disk
%
%   Cache invalidation: each entry stores the source file's modification
%   time; if mtime changes (re-ran batch_process), the entry is reloaded.
%
%   To force a full reload, delete cacheFile.

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

    % per-metric key strings so we don't recompute every iteration
    specKeys = cell(nMetrics,1);
    for k = 1:nMetrics
        specKeys{k} = makeSpecKey(metricSpecs(k));
    end

    % ----- progress UI -----
    wb = waitbar(0, 'Loading metric values...', ...
        'Name','plotAveragedMetrics: loading');
    try, set(wb,'WindowState','normal'); catch, end %#ok<CTCH>
    cleanup = onCleanup(@() safeCloseWb(wb)); %#ok<NASGU>

    hitCount  = 0;
    loadCount = 0;
    nDirty    = 0;

    % ----- main loop -----
    for i = 1:nFiles
        srcFile = state.files{i};

        % source-file mtime drives cache validity
        srcMtime = 0;
        info = dir(srcFile);
        if ~isempty(info), srcMtime = info(1).datenum; end

        if ishandle(wb)
            waitbar(i/nFiles, wb, sprintf( ...
                '[%d / %d]  %s\nCache hits: %d   Fresh loads: %d', ...
                i, nFiles, shortName(srcFile), hitCount, loadCount));
        end

        for k = 1:nMetrics
            ck = [srcFile '||' specKeys{k}];
            useCache = isKey(cacheMap, ck);
            if useCache
                entry = cacheMap(ck);
                useCache = (entry.mtime == srcMtime);
            end
            if useCache
                values(i,k) = entry.value;
                hitCount = hitCount + 1;
            else
                v = loadMetric(srcFile, metricSpecs(k));
                values(i,k) = v;
                cacheMap(ck) = struct('mtime', srcMtime, 'value', v);
                loadCount = loadCount + 1;
                nDirty    = nDirty + 1;
            end
        end
    end

    % ----- persist cache (only if anything changed) -----
    if nDirty > 0 && ~isempty(cacheFile)
        try
            save(cacheFile, 'cacheMap', '-v7.3');
        catch ME
            warning('loadAllMetricsCached:save', ...
                'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
function k = makeSpecKey(spec)
% Stable key string that uniquely identifies a metric definition.
    if isfield(spec,'channel') && isnumeric(spec.channel) && ~isnan(spec.channel)
        ch = num2str(spec.channel);
    else
        ch = '-';
    end
    k = sprintf('%s|%s|%s|%s', ...
        char(spec.suffix), char(spec.field), ...
        char(spec.aggregator), ch);
end

function safeCloseWb(wb)
    try
        if ishandle(wb), close(wb); end
    catch
        % already gone
    end
end

function s = shortName(fp)
    [d, b, e] = fileparts(char(fp));
    [~, parent] = fileparts(d);
    s = [parent filesep b e];
    if numel(s) > 70, s = ['...' s(end-66:end)]; end
end
