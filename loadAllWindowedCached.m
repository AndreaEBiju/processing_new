function [baselineByFile, windowByFile, hitCount, loadCount] = ...
        loadAllWindowedCached(state, metricSpecs, winSecs, cacheFile, opts)
%LOADALLWINDOWEDCACHED  Thin wrapper around loadAllSeriesCached: pulls
% every (file, metric) time series through the shared series cache,
% then derives baseline-mean + per-window means in memory.
%
%   [baselineByFile, windowByFile, hitCount, loadCount] = ...
%       loadAllWindowedCached(state, metricSpecs, winSecs, cacheFile)
%   ...                          = loadAllWindowedCached(..., opts)
%
% Returns
%   baselineByFile : nFiles x nMetrics       (mean of the full series)
%   windowByFile   : nFiles x nMetrics x nW  (mean over t-t0 <= winSecs(w))
%   hitCount, loadCount : cache statistics from loadAllSeriesCached
%
% Important: `cacheFile` should point at the SHARED series cache
% (gemsplots_series_cache.mat), not a separate windowed cache. This way
% plotWindowedMetrics, plotWindowedViolins, and plotTimeTraces all read
% from / write to the same cache, so whichever runs first populates it
% for the others — no duplicated file loads.

    if nargin < 5, opts = struct(); end

    % ----- load the actual time series (cached + parallel) -----
    [seriesByFile, hitCount, loadCount] = ...
        loadAllSeriesCached(state, metricSpecs, cacheFile, opts);

    % ----- derive the aggregated scalars in memory -----
    [nFiles, nMetrics] = size(seriesByFile);
    nW = numel(winSecs);
    baselineByFile = nan(nFiles, nMetrics);
    windowByFile   = nan(nFiles, nMetrics, nW);

    for i = 1:nFiles
        for k = 1:nMetrics
            s = seriesByFile{i,k};
            if isempty(s) || isempty(s.y) || isempty(s.t), continue; end
            Y = s.y;  T = s.t;
            baselineByFile(i,k) = mean(Y, 'omitnan');
            t0 = T(1);
            for w = 1:nW
                W = winSecs(w);
                if isinf(W), sel = true(size(T));
                else,        sel = (T - t0) <= W;
                end
                if any(sel)
                    windowByFile(i,k,w) = mean(Y(sel), 'omitnan');
                end
            end
        end
    end
end
