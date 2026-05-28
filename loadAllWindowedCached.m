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
% source file, groups metrics by the result file they target, then issues
% ONE load() per unique result file pulling every needed (seriesField,
% timeField) at once. v7 .mat files have to be scanned end-to-end, so
% batching brings load count from 11-per-source down to ~3-per-source.

    nMetrics = numel(metricSpecs);
    nW = numel(winSecs);
    baseline = nan(nMetrics, 1);
    windows  = nan(nMetrics, nW);

    % ----- per-metric: resolve target path + field names -----
    [d, base, ~] = fileparts(srcFile);
    stem = regexprep(base, '_blankmotion$', '');
    targetPath  = cell(nMetrics,1);
    seriesName  = cell(nMetrics,1);
    timeName    = cell(nMetrics,1);
    channelIdx  = nan(nMetrics,1);
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

    % ----- group metrics by target file path -----
    valid = ~cellfun(@isempty, targetPath);
    if ~any(valid), return; end
    paths = unique(targetPath(valid));

    wState = warning('off', 'MATLAB:load:variableNotFound');
    cleanupW = onCleanup(@() warning(wState));

    for p = 1:numel(paths)
        fp = paths{p};
        idx = find(strcmp(targetPath, fp));
        wantFields = unique([seriesName(idx); timeName(idx)]);
        wantFields = wantFields(~cellfun(@isempty, wantFields));

        % single load() with every needed variable
        try
            S = load(fp, wantFields{:});
        catch
            continue;
        end

        % resolve actual field names (case-insensitive fallback) once
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

        % distribute to each metric using this file
        for kk = 1:numel(idx)
            k  = idx(kk);
            sN = '';
            tN = '';
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
            Y = Y(1:n); T = T(1:n);

            baseline(k) = mean(Y, 'omitnan');
            t0 = T(1);
            for w = 1:nW
                W = winSecs(w);
                if isinf(W), sel = true(size(T));
                else,        sel = (T - t0) <= W;
                end
                if any(sel)
                    windows(k, w) = mean(Y(sel), 'omitnan');
                end
            end
        end
    end
end

function y = selectChannelLocal(x, ch)
% Same channel-picker as loadMetric / loadMetricSeries. Inlined so the
% parfeval worker doesn't need an extra function on the path.
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
