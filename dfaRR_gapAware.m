function out = dfaRR_gapAware(RR_intervals, RR_times, P)
% DFARR_GAPAWARE  Cardiac wrapper for dfaGapAware: RR-interval DFA that is
% aware of dropped beats between otherwise-valid RR_intervals entries.
%
%   out = dfaRR_gapAware(RR_intervals, RR_times, P)
%
% RR_intervals/RR_times already exclude any interval that itself spans a
% gap (computeValidRRIntervals), but consecutive surviving entries can
% still have a real-time jump between them wherever a beat was dropped.
% Every individual RR_intervals entry is legitimately valid; it's the
% TRANSITION between two valid entries that can be broken -- a validMask
% over RR_intervals cannot represent this. So a runs list is built here
% and passed to dfaGapAware via P.dfaRunsOverride (see
% DFA_IMPLEMENTATION_PLAN.md S1).
%
% RR_intervals - vector, seconds
% RR_times     - vector, seconds, onset time of each RR interval
% P            - (optional) pipeline_params() struct

if nargin < 3 || isempty(P); P = pipeline_params(); end
tol = P.dfaBeatGapTolSec; % ~1e-3 s; floating-point slop only, NOT physiological

n = numel(RR_intervals);
gapAfter = [RR_times(2:end) - RR_times(1:end-1) - RR_intervals(1:end-1); 0];
breakAfter = abs(gapAfter) > tol; % true at index i means run breaks between i, i+1

starts = [1; find(breakAfter(1:end-1)) + 1];
ends   = [find(breakAfter(1:end-1)); n];
runs = [starts, ends];

P.dfaRunsOverride = runs;
out = dfaGapAware(RR_intervals, true(n,1), P);
end
