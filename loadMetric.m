function v = loadMetric(srcFile, spec)
%LOADMETRIC  Pull one aggregated scalar metric from a result .mat file.
%
%   v = loadMetric(srcFile, spec)
%     srcFile : full path to the SOURCE recording file
%               (typically <stem>_blankmotion.mat). The result file is
%               located by stripping any trailing '_blankmotion' from
%               the basename and appending spec.suffix.
%     spec    : struct with fields
%        .suffix      string  e.g. '_HRBR.mat' (auto-appends '.mat')
%        .field       string  variable name inside the .mat
%        .aggregator  string  {auto,mean,median,max,min,first,last,sum,scalar}
%                             'auto' = take scalar as-is, else omitnan mean
%        .channel     numeric optional column index (NaN/[] = take whole field)
%
%   Returns NaN on any failure (missing file, missing field, wrong shape,
%   non-numeric data, etc.) so callers can safely build NaN-padded
%   matrices for box plots.

    v = NaN;
    if ~isstruct(spec) || ~isfield(spec,'suffix') || ~isfield(spec,'field')
        return;
    end

    % --- resolve result file path ---
    [d, base, ~] = fileparts(srcFile);
    stem = regexprep(base, '_blankmotion$', '');
    suffix = char(spec.suffix);
    if ~endsWith(suffix, '.mat'), suffix = [suffix '.mat']; end

    % Try multiple naming variants because batch_process saves outputs
    % differently depending on the source kind:
    %   <stem>_<suffix>             for baseline / explicit recovery files
    %   <stem>_recovery_<suffix>    for stim_rec sources that got split
    %                               (batch_process appends '_recovery' to
    %                               the condition before saving)
    candidates = {fullfile(d, [stem suffix])};
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
    if isempty(fp)
        return;
    end

    % --- load + extract field ---
    try
        S = load(fp, spec.field);
    catch
        return;
    end
    if ~isfield(S, spec.field), return; end
    x = S.(spec.field);
    if ~isnumeric(x) || isempty(x), return; end

    % --- optional channel selection ---
    if isfield(spec,'channel') && ~isempty(spec.channel) && ...
       isnumeric(spec.channel) && ~isnan(spec.channel)
        ch = round(double(spec.channel));
        x = selectChannel(x, ch);
        if isempty(x), return; end
    end

    % --- aggregate ---
    aggregator = lower(strtrim(char(getf(spec,'aggregator'))));
    if isempty(aggregator), aggregator = 'auto'; end
    xv = double(x(:));
    switch aggregator
        case 'auto'
            if isscalar(xv), v = xv;
            else,            v = mean(xv, 'omitnan');
            end
        case {'scalar','first'},  v = xv(1);
        case 'last',              v = xv(end);
        case 'mean',              v = mean   (xv, 'omitnan');
        case 'median',            v = median (xv, 'omitnan');
        case 'max',               v = max    (xv, [], 'omitnan');
        case 'min',               v = min    (xv, [], 'omitnan');
        case 'sum',               v = sum    (xv, 'omitnan');
        otherwise,                v = mean   (xv, 'omitnan');
    end
end


% ==========================================================================
function y = selectChannel(x, ch)
% Pick a single channel from a 2D matrix. Channels are columns when ncols
% is small (<=8) or smaller than nrows; otherwise rows. Returns [] if the
% requested channel doesn't exist.
    y = [];
    if ~ismatrix(x), return; end
    [r, c] = size(x);
    if c == 1 && r > 1                              % column vector
        if ch == 1, y = x; end
        return;
    end
    if r == 1 && c > 1                              % row vector
        if ch == 1, y = x.'; end
        return;
    end
    if c <= 8 && c <= r                             % [time x channels]
        if ch <= c, y = x(:, ch); end
    elseif r <= 8 && r < c                          % [channels x time]
        if ch <= r, y = x(ch, :).'; end
    else                                            % ambiguous — try cols
        if ch <= c, y = x(:, ch); end
    end
end

function v = getf(s, name)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = []; end
end
