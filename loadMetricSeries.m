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

    % --- load both fields (suppress MATLAB's per-call "not found" warning;
    %     we'll emit a single, more useful diagnostic if anything is missing)
    wState = warning('off', 'MATLAB:load:variableNotFound');
    try
        S = load(fp, char(spec.seriesField), char(spec.timeField));
    catch
        warning(wState);
        fprintf('[loadMetricSeries] FAIL load %s\n', fp);
        return;
    end
    warning(wState);

    % --- resolve actual variable names (case-insensitive fallback) ---
    sName = resolveField(S, spec.seriesField, fp);
    tName = resolveField(S, spec.timeField,   fp);

    if isempty(sName) || isempty(tName)
        % final attempt: whos the whole file to find a case/punct-insensitive
        % match. Slower (reads the file header) but useful for diagnosis.
        try, w = whos('-file', fp); catch, w = []; end
        if ~isempty(w)
            if isempty(sName), sName = matchByCanonical(spec.seriesField, w, fp); end
            if isempty(tName), tName = matchByCanonical(spec.timeField,   w, fp); end
            % If we found names via whos, do a targeted re-load
            if ~isempty(sName) && ~isempty(tName)
                try
                    S = load(fp, sName, tName);
                catch
                    S = struct();
                end
            end
        end
    end

    if isempty(sName) || ~isfield(S, sName) || ...
       isempty(tName) || ~isfield(S, tName)
        avail = '';
        if ~isempty(S)
            avail = strjoin(fieldnames(S),', ');
        end
        fprintf(['[loadMetricSeries] MISSING in %s\n' ...
                 '   wanted:   series=''%s''  time=''%s''\n' ...
                 '   loaded:   %s\n'], ...
            fp, char(spec.seriesField), char(spec.timeField), avail);
        return;
    end

    Y = S.(sName);
    T = S.(tName);
    if ~isnumeric(Y) || ~isnumeric(T) || isempty(Y) || isempty(T)
        fprintf('[loadMetricSeries] NON-NUMERIC %s in %s\n', sName, fp);
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
function name = resolveField(S, wanted, ~)
% Return the actual fieldname in struct S matching `wanted` exactly,
% or case-insensitively if there's exactly one such match. Otherwise ''.
    name = '';
    if isempty(S) || ~isstruct(S), return; end
    fns = fieldnames(S);
    if isempty(fns), return; end
    w = char(wanted);
    if any(strcmp(fns, w))
        name = w; return;
    end
    hit = strcmpi(fns, w);
    if sum(hit) == 1
        name = fns{hit};
    end
end

function name = matchByCanonical(wanted, w, ~)
% Find a variable in whos-listing w whose canonical form (lowercase,
% non-alphanumeric stripped) matches `wanted`. Returns '' if 0 or >1.
    name = '';
    if isempty(w), return; end
    canon = @(s) lower(regexprep(char(s), '[^a-zA-Z0-9]', ''));
    target = canon(wanted);
    hits = false(numel(w),1);
    for k = 1:numel(w)
        hits(k) = strcmp(canon(w(k).name), target);
    end
    if sum(hits) == 1
        idx = find(hits,1);
        name = w(idx).name;
    end
end

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
