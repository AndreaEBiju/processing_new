function [seriesByFile, hitCount, loadCount] = ...
        loadAllSeriesCached(state, metricSpecs, cacheFile, opts)
%LOADALLSERIESCACHED  Full time series (Y, T) per (file, metric) with disk
% cache and parfeval parallel load. Mirrors loadAllWindowedCached's
% pattern but stores the actual vectors (not just scalars) so downstream
% plotting can do anything with the trace.
%
%   [seriesByFile, hits, loads] = loadAllSeriesCached(state, metricSpecs, cacheFile)
%   ...                          = loadAllSeriesCached(..., opts)
%
% Returns
%   seriesByFile : nFiles x nMetrics cell array.  Each cell is either []
%                  (couldn't load) or struct with fields .y .t (column
%                  vectors in the file's own units).
%
% Cache: keyed by (filepath, specKey), mtime-invalidated. Stored at the
% caller-supplied path, separate from the scalar/windowed cache.

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'useParallel'), opts.useParallel = true; end
    if ~isfield(opts,'nWorkers'),    opts.nWorkers    = []; end

    nFiles   = numel(state.files);
    nMetrics = numel(metricSpecs);
    seriesByFile = cell(nFiles, nMetrics);

    % ----- restore cache -----
    cacheMap = containers.Map('KeyType','char','ValueType','any');
    if ~isempty(cacheFile) && exist(cacheFile,'file')
        try
            c = load(cacheFile);
            if isfield(c, 'cacheMap'), cacheMap = c.cacheMap; end
        catch ME
            warning('loadAllSeriesCached:restore', ...
                'Cache restore failed (%s); starting fresh.', ME.message);
        end
    end

    specKeys = cell(nMetrics,1);
    for k = 1:nMetrics
        specKeys{k} = makeSpecKey(metricSpecs(k));
    end

    % ----- pass 1: cache hits + decide what to load -----
    hitCount = 0;
    needsLoad  = false(nFiles, 1);
    fileMtimes = zeros(nFiles, 1);
    for i = 1:nFiles
        info = dir(state.files{i});
        if ~isempty(info), fileMtimes(i) = info(1).datenum; end
        for k = 1:nMetrics
            ck = [state.files{i} '||' specKeys{k}];
            if isKey(cacheMap, ck) && cacheMap(ck).mtime == fileMtimes(i) ...
                    && ~isempty(cacheMap(ck).y) && ~isempty(cacheMap(ck).t)
                seriesByFile{i,k} = struct('y', cacheMap(ck).y, 't', cacheMap(ck).t);
                hitCount = hitCount + 1;
            else
                needsLoad(i) = true;
            end
        end
    end

    filesToLoad = find(needsLoad);
    nLoad = numel(filesToLoad);
    loadCount = 0;
    if nLoad == 0, return; end

    % ----- progress UI -----
    wb = waitbar(0, sprintf('Loading time series from %d file(s)...', nLoad), ...
        'Name','plotTimeTraces: loading');
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
            warning('loadAllSeriesCached:parpool', ...
                'Could not start parpool: %s. Falling back to serial.', ME.message);
            canParallel = false;
        end
    end

    nDone = 0;
    if canParallel
        futures = parallel.FevalFuture.empty(0,1);
        for j = 1:nLoad
            i = filesToLoad(j);
            futures(j,1) = parfeval(@computeFileSeries, 1, ...
                state.files{i}, metricSpecs);
        end
        for j = 1:nLoad
            try
                [completedIdx, fileSeries] = fetchNext(futures);
            catch ME
                warning('loadAllSeriesCached:parfeval', ...
                    'parfeval failed at task %d: %s', j, ME.message);
                cancel(futures);
                break;
            end
            i = filesToLoad(completedIdx);
            for k = 1:nMetrics
                if ~isempty(fileSeries{k})
                    seriesByFile{i,k} = fileSeries{k};
                    cacheMap([state.files{i} '||' specKeys{k}]) = ...
                        struct('mtime', fileMtimes(i), ...
                               'y',     fileSeries{k}.y, ...
                               't',     fileSeries{k}.t);
                    loadCount = loadCount + 1;
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
            fileSeries = computeFileSeries(state.files{i}, metricSpecs);
            for k = 1:nMetrics
                if ~isempty(fileSeries{k})
                    seriesByFile{i,k} = fileSeries{k};
                    cacheMap([state.files{i} '||' specKeys{k}]) = ...
                        struct('mtime', fileMtimes(i), ...
                               'y',     fileSeries{k}.y, ...
                               't',     fileSeries{k}.t);
                    loadCount = loadCount + 1;
                end
            end
            nDone = nDone + 1;
            if ishandle(wb)
                [~, b, e] = fileparts(state.files{i});
                waitbar(nDone/nLoad, wb, sprintf('[%d / %d] %s', nDone, nLoad, [b e]));
            end
        end
    end

    if ~isempty(cacheFile)
        try
            save(cacheFile, 'cacheMap', '-v7.3');
        catch ME
            warning('loadAllSeriesCached:save', 'Cache save failed: %s', ME.message);
        end
    end
end


% ==========================================================================
function fileSeries = computeFileSeries(srcFile, metricSpecs)
% Same batched-load pattern as loadAllWindowedCached/computeFileValues:
% one load() per unique result file with all needed variables.
    nMetrics = numel(metricSpecs);
    fileSeries = cell(nMetrics, 1);

    [d, base, ~] = fileparts(srcFile);
    stem = regexprep(base, '_blankmotion$', '');
    targetPath = cell(nMetrics,1);
    seriesName = cell(nMetrics,1);
    timeName   = cell(nMetrics,1);
    channelIdx = nan(nMetrics,1);
    for k = 1:nMetrics
        spec = metricSpecs(k);
        if isempty(spec.seriesField) || isempty(spec.timeField), continue; end
        suffix = char(spec.suffix);
        if ~endsWith(suffix, '.mat'), suffix = [suffix '.mat']; end
        candidates = {fullfile(d, [base suffix]), fullfile(d, [stem suffix])};
        if contains(stem, 'stim_rec', 'IgnoreCase', true) && ...
           ~endsWith(stem, '_recovery', 'IgnoreCase', true) && ...
           ~endsWith(stem, '_stim',     'IgnoreCase', true)
            candidates{end+1} = fullfile(d, [stem '_recovery' suffix]); %#ok<AGROW>
        end
        for ci = 1:numel(candidates)
            if exist(candidates{ci}, 'file')
                targetPath{k} = candidates{ci}; break;
            end
        end
        seriesName{k} = char(spec.seriesField);
        timeName{k}   = char(spec.timeField);
        if isfield(spec,'channel') && isnumeric(spec.channel) && ~isnan(spec.channel)
            channelIdx(k) = round(double(spec.channel));
        end
    end

    valid = ~cellfun(@isempty, targetPath);
    if ~any(valid), return; end
    paths = unique(targetPath(valid));

    wState = warning('off', 'MATLAB:load:variableNotFound');
    cleanupW = onCleanup(@() warning(wState)); %#ok<NASGU>

    for p = 1:numel(paths)
        fp = paths{p};
        idx = find(strcmp(targetPath, fp));
        wantFields = unique([seriesName(idx); timeName(idx)]);
        wantFields = wantFields(~cellfun(@isempty, wantFields));

        try
            S = load(fp, wantFields{:});
        catch
            continue;
        end

        actualName = containers.Map('KeyType','char','ValueType','char');
        fns = fieldnames(S);
        for ff = 1:numel(wantFields)
            w = wantFields{ff};
            if any(strcmp(fns, w))
                actualName(w) = w;
            else
                hit = strcmpi(fns, w);
                if sum(hit) == 1, actualName(w) = fns{find(hit,1)}; end
            end
        end

        for kk = 1:numel(idx)
            k = idx(kk);
            sN = ''; tN = '';
            if isKey(actualName, seriesName{k}), sN = actualName(seriesName{k}); end
            if isKey(actualName, timeName{k}),   tN = actualName(timeName{k});   end
            if isempty(sN) || isempty(tN), continue; end

            Y = S.(sN);
            T = S.(tN);
            if ~isnumeric(Y) || ~isnumeric(T) || isempty(Y) || isempty(T), continue; end

            if ~isnan(channelIdx(k))
                Y = selectChannelLocal(Y, channelIdx(k));
                if isempty(Y), continue; end
            end

            Y = Y(:); T = T(:);
            n = min(numel(Y), numel(T));
            if n < 2, continue; end
            fileSeries{k} = struct('y', Y(1:n), 't', T(1:n));
        end
    end
end


function y = selectChannelLocal(x, ch)
    y = [];
    if ~ismatrix(x), return; end
    [r, c] = size(x);
    if c == 1 && r > 1, if ch == 1, y = x; end, return; end
    if r == 1 && c > 1, if ch == 1, y = x.'; end, return; end
    if c <= 8 && c <= r
        if ch <= c, y = x(:, ch); end
    elseif r <= 8 && r < c
        if ch <= r, y = x(ch, :).'; end
    else
        if ch <= c, y = x(:, ch); end
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
