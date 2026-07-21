function out = dfaGapAware(x, validMask, P)
% DFAGAPAWARE Gap-aware Detrended Fluctuation Analysis.
%
% out = dfaGapAware(x, validMask, P)
%
% x         - series to analyze, e.g. RR_intervals (beat-indexed, cardiac)
%             or fr_hz (1 Hz-binned firing rate, nerve). One value per
%             index; index spacing need not be uniform in real time (RR
%             series) as long as x(i) is the i-th valid measurement.
% validMask - logical, same length as x. true = usable sample/beat/bin.
%             false = inside a blanked/excluded region (motion artifact,
%             QC floor failure, etc). Pass true(size(x)) if x has already
%             had invalid points removed and every remaining point is
%             known-contiguous in real time (NOT recommended for motion-
%             artifact data -- see fullSegmentsToLocalSegments.m for
%             translating sample-index blank segments into this mask).
% P         - (optional) pipeline_params() struct. Reads P.dfa* fields
%             (see pipeline_params.m Step 7 block). Falls back to
%             built-in defaults if omitted or fields are missing.
%             P.dfaRunsOverride (optional): an Nx2 [start end] index list
%             to use INSTEAD of deriving runs from validMask. Use this
%             when gap structure isn't naturally per-sample (e.g. cardiac
%             RR intervals, where every beat is individually valid but a
%             dropped beat invalidates the TRANSITION between two
%             surviving beats -- validMask can't express that, a runs
%             list can). When set, validMask's only remaining job is
%             sizing/demeaning; pass true(size(x)) for it in that case.
%
% out fields:
%   nVals        - box sizes used (after truncation), in units of x's index
%   F            - fluctuation function F(n), same length as nVals
%   nWindows     - number of valid detrending windows pooled at each n
%                  (QC diagnostic -- NOT a physiological metric)
%   alpha1       - short-term scaling exponent (P.dfaShortRange)
%   alpha2       - long-term scaling exponent (n > P.dfaLongRangeMinScale)
%   alphaFull    - single exponent over the full retained nVals range
%   R2_1, R2_2, R2_full - fit quality (R^2) for each of the above
%   nCross       - crossover scale where the two local fits meet (interp),
%                  NaN if alpha1/alpha2 could not both be estimated
%   localAlpha   - struct with .n, .slope: scale-dependent local exponent
%                  (3-point moving slope of log F vs log n), for plotting
%                  a temporal-spectrum-of-scale-exponents style curve
%   runs         - Nruns x 2 index list of contiguous valid runs used
%   excludedScales - box sizes dropped for having too few pooled windows
%
% Method: pooled valid-window F^2(n) accumulation across contiguous clean
% sub-segments (Ma, Bartsch, Bernaola-Galvan, Yoneyama & Ivanov, Phys. Rev.
% E 81, 031101, 2010; "blocks adjustment", arXiv:0708.1628). Each
% contiguous valid run is treated as its own finite series; non-overlapping
% detrending windows of size n are tiled from BOTH ends of each run (as in
% the classical DFA convention) so leftover samples at one end aren't
% silently dropped. Windows are pooled across ALL runs before averaging,
% so short and long clean stretches all contribute rather than requiring
% one unbroken record. No interpolation is performed across gaps -- see
% conversation rationale: interpolation across gaps of more than a few
% seconds risks fabricating the long-range structure DFA is meant to
% measure.
%
% Global mean/demeaning uses only valid samples, computed once over the
% whole series (not per-run), matching the standard DFA convention of a
% single global profile origin.

if nargin < 3 || isempty(P)
    P = struct();
end
P = fillDefaultDfaParams(P);

x = x(:);
validMask = logical(validMask(:));
N = numel(x);
if numel(validMask) ~= N
    error('dfaGapAware:sizeMismatch', 'x and validMask must be the same length.');
end

% ---- global demeaning over valid samples only ----
xValid = x(validMask);
if numel(xValid) < P.dfaMinRunLen * 4
    warning('dfaGapAware:tooFewValid', ...
        'Only %d valid points (%d required); DFA output will be unreliable or NaN.', ...
        numel(xValid), P.dfaMinRunLen * 4);
end
mu = mean(xValid, 'omitnan');
xd = x - mu; % demeaned; invalid entries never read past this point

% ---- contiguous valid runs ----
% Two ways in: (a) a validMask, from which runs are found via the same
% diff-based idiom as step2_noise_sigma.m's flat_run_mask, or (b) a
% pre-built Nx2 runs list via P.dfaRunsOverride, for callers whose gap
% structure isn't naturally a per-sample validMask -- e.g. the cardiac RR
% case, where every individual beat is "valid" but the TRANSITION between
% two valid beats can be invalid (a dropped beat), which validMask alone
% can't represent. See DFA_IMPLEMENTATION_PLAN.md S5 for the cardiac
% wrapper's construction of this override.
if isfield(P, 'dfaRunsOverride') && ~isempty(P.dfaRunsOverride)
    runs = P.dfaRunsOverride;
    runs = runs(runs(:,2) - runs(:,1) + 1 >= P.dfaMinRunLen, :);
else
    d = diff([0; validMask(:); 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    runs = [starts, ends];
    runs = runs(runs(:,2) - runs(:,1) + 1 >= P.dfaMinRunLen, :);
end

if isempty(runs)
    error('dfaGapAware:noValidRuns', ...
        'No contiguous valid run reaches P.dfaMinRunLen (%d).', P.dfaMinRunLen);
end
runLens = runs(:,2) - runs(:,1) + 1;
totalValidLen = sum(runLens);

fprintf('[dfaGapAware] %d valid run(s), total %d samples/beats (%.1f%% of %d).\n', ...
    size(runs,1), totalValidLen, 100*totalValidLen/N, N);

% ---- box sizes ----
nMin = P.dfaMinScale;
nMax = max(nMin + 1, floor(P.dfaMaxScaleFrac * max(runLens)));
nVals = unique(round(logspace(log10(nMin), log10(nMax), P.dfaNScales)));

F2sum = zeros(numel(nVals), 1);
F2n   = zeros(numel(nVals), 1); % pooled window count per scale

for k = 1:numel(nVals)
    n = nVals(k);
    t = (1:n)';
    for r = 1:size(runs,1)
        i0 = runs(r,1); i1 = runs(r,2);
        L = i1 - i0 + 1;
        if L < n; continue; end
        y = cumsum(xd(i0:i1)); % local profile, restarts at each run
        Ns = floor(L / n);
        % forward tiling
        for v = 1:Ns
            seg = y((v-1)*n+1 : v*n);
            p = polyfit(t, seg, P.dfaDetrendOrder);
            resid = seg - polyval(p, t);
            F2sum(k) = F2sum(k) + mean(resid.^2);
            F2n(k) = F2n(k) + 1;
        end
        % reverse tiling to use the leftover tail (classic DFA convention)
        rem = mod(L, n);
        if rem > 0 && Ns >= 1
            yRev = flipud(y);
            for v = 1:Ns
                seg = yRev((v-1)*n+1 : v*n);
                p = polyfit(t, seg, P.dfaDetrendOrder);
                resid = seg - polyval(p, t);
                F2sum(k) = F2sum(k) + mean(resid.^2);
                F2n(k) = F2n(k) + 1;
            end
        end
    end
end

F = sqrt(F2sum ./ max(F2n, 1));
F(F2n < P.dfaMinWindowsPerScale) = NaN;

excludedScales = nVals(F2n < P.dfaMinWindowsPerScale | ~isfinite(F));
keep = isfinite(F) & F > 0;
nValsKept = nVals(keep);
Fkept = F(keep);
nWindowsKept = F2n(keep);

if numel(nValsKept) < 3
    warning('dfaGapAware:tooFewScales', ...
        'Only %d usable scale(s) survived the window-count floor; check gap burden.', ...
        numel(nValsKept));
end

logn = log10(nValsKept);
logF = log10(Fkept);

[alpha1, R2_1] = fitRange(logn, logF, nValsKept, P.dfaShortRange);
[alpha2, R2_2] = fitRange(logn, logF, nValsKept, [P.dfaLongRangeMinScale, inf]);
[alphaFull, R2_full] = fitRange(logn, logF, nValsKept, [nValsKept(1), nValsKept(end)]);

% crossover: intersection of the two local linear fits (if both exist)
nCross = NaN;
if ~isnan(alpha1) && ~isnan(alpha2) && numel(nValsKept) >= 4
    idxShort = nValsKept >= P.dfaShortRange(1) & nValsKept <= P.dfaShortRange(2);
    idxLong  = nValsKept > P.dfaLongRangeMinScale;
    if any(idxShort) && any(idxLong)
        p1 = polyfit(logn(idxShort), logF(idxShort), 1);
        p2 = polyfit(logn(idxLong),  logF(idxLong),  1);
        if abs(p1(1) - p2(1)) > eps
            logNx = (p2(2) - p1(2)) / (p1(1) - p2(1));
            nCross = 10^logNx;
        end
    end
end

% local (scale-dependent) slope: centered 3-point finite difference in log-log space
localN = nValsKept(2:end-1);
localSlope = nan(numel(localN), 1);
for i = 2:numel(nValsKept)-1
    localSlope(i-1) = (logF(i+1) - logF(i-1)) / (logn(i+1) - logn(i-1));
end

out = struct();
out.nVals = nValsKept;
out.F = Fkept;
out.nWindows = nWindowsKept;
out.alpha1 = alpha1;
out.alpha2 = alpha2;
out.alphaFull = alphaFull;
out.R2_1 = R2_1;
out.R2_2 = R2_2;
out.R2_full = R2_full;
out.nCross = nCross;
out.localAlpha = struct('n', localN, 'slope', localSlope);
out.runs = runs;
out.excludedScales = excludedScales;
out.mu = mu;
out.P = P;

fprintf(['[dfaGapAware] alpha1=%.3f (R2=%.2f) | alpha2=%.3f (R2=%.2f) | ' ...
    'alphaFull=%.3f (R2=%.2f) | crossover n~%.1f | %d/%d scales excluded (window floor)\n'], ...
    alpha1, R2_1, alpha2, R2_2, alphaFull, R2_full, nCross, ...
    numel(excludedScales), numel(nVals));

end

% ========================================================================
function [a, R2] = fitRange(logn, logF, nVals, range)
% logn/logF/nVals can arrive with mismatched row/column orientation (nVals
% traces back to a row-shaped logspace(), F to a column-shaped zeros()) --
% force a consistent orientation here so ssRes below is a scalar elementwise
% sum, not an implicit-expansion (row-vs-column) broadcast into a matrix.
logn = logn(:); logF = logF(:); nVals = nVals(:);
idx = nVals >= range(1) & nVals <= range(2);
if sum(idx) < 3
    a = NaN; R2 = NaN;
    return;
end
p = polyfit(logn(idx), logF(idx), 1);
a = p(1);
fit = polyval(p, logn(idx));
ssRes = sum((logF(idx) - fit).^2);
ssTot = sum((logF(idx) - mean(logF(idx))).^2);
if ssTot <= 0
    R2 = NaN;
else
    R2 = 1 - ssRes/ssTot;
end
end

% ========================================================================
function P = fillDefaultDfaParams(P)
% Mirrors the pipeline_params.m Step 7 block; used standalone if a P
% struct/field is not supplied (e.g. calling dfaGapAware outside the repo).
d = struct( ...
    'dfaMinScale',           4, ...
    'dfaMaxScaleFrac',       0.25, ...
    'dfaNScales',            20, ...
    'dfaShortRange',         [4 16], ...
    'dfaLongRangeMinScale',  16, ...
    'dfaMinWindowsPerScale', 8, ...
    'dfaDetrendOrder',       1, ...
    'dfaMinRunLen',          4);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(P, fn{i}) || isempty(P.(fn{i}))
        P.(fn{i}) = d.(fn{i});
    end
end
end
