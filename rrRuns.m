function runs = rrRuns(RR_intervals, RR_times, tol)
%RRRUNS  Contiguous runs of RR intervals with no gap between them.
%
%   runs = RRRUNS(RR_intervals, RR_times) returns an [nRuns x 2] matrix of
%   [startIdx, endIdx] into RR_intervals. Two consecutive intervals belong to the
%   same run only when the second starts where the first ended, to within tol.
%
%   WHY THIS EXISTS AS ITS OWN FUNCTION
%   The same computation was already inside dfaRR_gapAware.m and nowhere else, so
%   every other successive-difference statistic - RMSSD, pNN5, SD1, SD2, SampEn,
%   ApEn - was differencing straight across blanked gaps. A difference taken across
%   a gap is not a successive difference: it is the difference between two intervals
%   that are not adjacent in time, and it enters the statistic as if it were.
%
%   The bug direction matters. Blanking removes beats, so a gap-spanning "interval"
%   is long, and its difference against its neighbours is large - so RMSSD, SD1 and
%   pNN5 are inflated by exactly the spans the detector removed. A better detector
%   produces more, shorter, better-placed gaps and therefore looks worse. That is
%   the instrument problem task 08 exists to fix.
%
%   INPUTS
%     RR_intervals : [n x 1] RR intervals, seconds
%     RR_times     : [n x 1] time of the beat each interval starts at, seconds
%     tol          : (optional) gap tolerance in seconds. Default comes from
%                    pipeline_params().dfaBeatGapTolSec, which is ~1e-3 s and is
%                    floating-point slop, NOT a physiological tolerance.
%
%   OUTPUT
%     runs : [nRuns x 2] of [startIdx, endIdx], inclusive, into RR_intervals.
%
%   See also DFARR_GAPAWARE, HR_BR_HRVANALYSIS_NEW.

    if nargin < 3 || isempty(tol)
        P = pipeline_params();
        tol = P.dfaBeatGapTolSec;
    end

    RR_intervals = RR_intervals(:);
    RR_times     = RR_times(:);
    n = numel(RR_intervals);

    if n == 0
        runs = zeros(0, 2);
        return;
    end
    if n == 1
        runs = [1 1];
        return;
    end

    % gapAfter(i) is how far interval i+1 starts after interval i ended.
    gapAfter   = [RR_times(2:end) - RR_times(1:end-1) - RR_intervals(1:end-1); 0];
    breakAfter = abs(gapAfter) > tol;

    starts = [1; find(breakAfter(1:end-1)) + 1];
    ends   = [find(breakAfter(1:end-1)); n];
    runs   = [starts, ends];
end
