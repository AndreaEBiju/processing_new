function out = hrvRunsAware(RR_intervals, RR_times, P)
%HRVRUNSAWARE  HRV statistics that never difference across a blanked gap.
%
%   out = HRVRUNSAWARE(RR_intervals, RR_times) returns a struct with fields
%   hrv, rmssd, pnn5, sd1, sd2, sampEn, appxEn, plus the bookkeeping the caller
%   needs to report what was used: nRuns, nUsedDiffs, nSplicedDiffs.
%
%   WHAT CHANGES, AND WHY IT MATTERS
%   Successive-difference statistics (RMSSD, pNN5, SD1, and SD2 through SD1) are
%   computed over adjacent pairs *within* a run only. Previously diff() ran over the
%   whole concatenated series, so every blanked gap contributed one spurious
%   difference between two intervals that are not adjacent in time.
%
%   The error is one-directional. Blanking removes beats, so the interval spanning a
%   gap is long and its differences against its neighbours are large: RMSSD, SD1 and
%   pNN5 all come out inflated, and inflated in proportion to how much was blanked.
%   A better detector blanks more precisely, produces more gaps, and therefore
%   reports worse HRV - which is the instrument, not the physiology.
%
%   SDNN (hrv) is NOT a successive-difference statistic: it is the spread of the
%   intervals themselves, and every retained interval is a real measurement whatever
%   its neighbours. It is computed over all intervals, unchanged.
%
%   SampEn and ApEn need a contiguous series and have no pooled definition, so they
%   are computed per run and combined as a length-weighted mean over runs with
%   enough intervals to support them. Runs too short are excluded rather than padded.
%
%   INPUTS
%     RR_intervals : [n x 1] seconds
%     RR_times     : [n x 1] seconds, the beat each interval starts at
%     P            : (optional) pipeline_params() struct
%
%   See also RRRUNS, DFARR_GAPAWARE.

    if nargin < 3 || isempty(P); P = pipeline_params(); end

    RR_intervals = RR_intervals(:);
    RR_times     = RR_times(:);

    out = struct('hrv', NaN, 'rmssd', NaN, 'pnn5', NaN, 'sd1', NaN, 'sd2', NaN, ...
                 'sampEn', NaN, 'appxEn', NaN, 'nRuns', 0, ...
                 'nUsedDiffs', 0, 'nSplicedDiffs', 0);

    if numel(RR_intervals) < 2
        return;
    end

    runs = rrRuns(RR_intervals, RR_times, P.dfaBeatGapTolSec);
    out.nRuns = size(runs, 1);

    % Successive differences, within runs only.
    diffRR = [];
    for r = 1:size(runs, 1)
        seg = RR_intervals(runs(r,1):runs(r,2));
        if numel(seg) >= 2
            diffRR = [diffRR; diff(seg)]; %#ok<AGROW>
        end
    end

    out.nUsedDiffs    = numel(diffRR);
    % How many differences the old code took that this one refuses: one per break.
    out.nSplicedDiffs = (numel(RR_intervals) - 1) - numel(diffRR);

    % SDNN uses every interval: each is a real measurement regardless of adjacency.
    out.hrv = std(RR_intervals, 'omitnan');

    if isempty(diffRR)
        return;
    end

    out.rmssd = sqrt(mean(diffRR.^2, 'omitnan'));

    absdiffRR = abs(diffRR .* 1000);
    out.pnn5  = (sum(absdiffRR > 5, 'omitnan') / numel(absdiffRR)) * 100;

    out.sd1 = sqrt(0.5) * std(diffRR, 'omitnan');
    out.sd2 = sqrt(max(0, 2*out.hrv^2 - out.sd1^2));

    % Entropies: per run, length-weighted. A run shorter than this cannot support
    % an m=2 template count and is excluded rather than padded into existence.
    minRunForEntropy = 10;
    sampAcc = 0; appxAcc = 0; wAcc = 0;
    for r = 1:size(runs, 1)
        seg = RR_intervals(runs(r,1):runs(r,2));
        if numel(seg) >= minRunForEntropy
            s = sampleEntropy(seg);
            a = approximateEntropy(seg);
            if isfinite(s) && isfinite(a)
                w       = numel(seg);
                sampAcc = sampAcc + w*s;
                appxAcc = appxAcc + w*a;
                wAcc    = wAcc + w;
            end
        end
    end
    if wAcc > 0
        out.sampEn = sampAcc / wAcc;
        out.appxEn = appxAcc / wAcc;
    end
end
