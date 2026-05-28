function [y, t] = loadMetricSeries(srcFile, spec)
%LOADMETRICSERIES  Pull a time-series + matching time vector from a
% result .mat file. Mirrors loadMetric's path resolution exactly so
% queue entries pointing at any of the supported source-file shapes
% (baseline, post-split recovery, original stim_rec) all work.
%
%   [y, t] = loadMetricSeries(srcFile, spec)
%
%   spec must define
%     .suffix      e.g. '_HRBR.mat'
%     .seriesField e.g. 'heartRateSeries'
%     .timeField   e.g. 'metrics_t'
%     .channel     optional column index for multi-channel matrices;
%                  NaN/empty = take whole series.
%
%   Returns
%     y : column vector of values  (size N x 1)  or [] if nothing found
%     t : column vector of times   (size N x 1)  in the file's own units
%         (typically seconds; loadMetricSeries does NOT unit-convert).
%
%   Both are returned empty on any failure (missing file, missing field,
%   wrong shape) so callers can NaN-pad safely.

    y = []; t = [];
    if ~isstruct(spec) ...
       || ~isfield(spec,'suffix') ...
       || ~isfield(spec,'seriesField') || isempty(spec.seriesField) ...
       || ~isfield(spec,'timeField')   || isempty(spec.timeField)
        return;
    end

    % --- resolve result file path (same logic as loadMetric) ---
    [d, base, ~] = fileparts(srcFile);
    stem = regexprep(base, '_blankmotion$', '');
    suffix = char(spec.suffix);
    if ~endsWith(suffix, '.mat'), suffix = [suffix '.mat']; end

    candidates = {fullfile(d, [base suffix]), ...
                  fullfile(d, [stem suffix])};
    if contains(stem, 'stim_rec', 'IgnoreCase', true) && ...
       ~endsWith(stem, '_recovery', 'IgnoreCase', true) && ...
       ~endsWith(stem, '_stim',     'IgnoreCase', true)
        candidates{end+1} = fullfile(d, [stem '_recovery' suffix]);
    end

    fp = '';
    for ci = 1:numel(candidates)
        if exist(candidates{ci}, 'file')
            fp = candidates{ci};
            break;
        end
    end
    if isempty(fp), return; end

    % --- load both fields ---
    try
        S = load(fp, char(spec.seriesField), char(spec.timeField));
    catch
        return;
    end
    if ~isfield(S, spec.seriesField) || ~isfield(S, spec.timeField)
        return;
    end
    Y = S.(spec.seriesField);
    T = S.(spec.timeField);
    if ~isnumeric(Y) || ~isnumeric(T) || isempty(Y) || isempty(T)
        return;
    end

    % --- channel selection on Y if needed ---
    if isfield(spec,'channel') && ~isempty(spec.channel) && ...
       isnumeric(spec.channel) && ~isnan(spec.channel)
        ch = round(double(spec.channel));
        Y = selectChannel(Y, ch);
        if isempty(Y), return; end
    end

    % --- align lengths (sometimes series and t differ by 1) ---
    Y = Y(:);
    T = T(:);
    n = min(numel(Y), numel(T));
    if n < 2, return; end
    y = Y(1:n);
    t = T(1:n);
end


% ==========================================================================
function y = selectChannel(x, ch)
% Same channel-picker as loadMetric — choose a column when ncols is the
% smaller, narrow dim; otherwise try rows.
    y = [];
    if ~ismatrix(x), return; end
    [r, c] = size(x);
    if c == 1 && r > 1
        if ch == 1, y = x; end
        return;
    end
    if r == 1 && c > 1
        if ch == 1, y = x.'; end
        return;
    end
    if c <= 8 && c <= r
        if ch <= c, y = x(:, ch); end
    elseif r <= 8 && r < c
        if ch <= r, y = x(ch, :).'; end
    else
        if ch <= c, y = x(:, ch); end
    end
end
