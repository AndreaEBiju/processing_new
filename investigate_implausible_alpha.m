function summary = investigate_implausible_alpha(datasetRoot)
% INVESTIGATE_IMPLAUSIBLE_ALPHA  Follow-up on the 31 cardiac files the real
% queue-scoped run_dfa_batch pass flagged with dfa_alpha1/alpha2 outside
% the ~0.5-1.3 plausible range for rat HRV. Per DFA_IMPLEMENTATION_PLAN.md
% T15: "flag, do not silently accept, anything far outside this as a
% likely masking bug rather than a real finding" -- this script is that
% follow-through.
%
% For each flagged file, recomputes dfaRR_gapAware on its saved
% RR_intervals/RR_times to pull the full diagnostic (nWindows per scale,
% run fragmentation) that isn't saved in the .mat's two scalar dfa_alpha1/
% dfa_alpha2 fields, and reports per file:
%   - total RR-series duration and % valid
%   - number of runs / mean run length (fragmentation)
%   - the minimum pooled-window count among scales used in the alpha2 fit
%     range -- a fit resting on scales that barely clear
%     P.dfaMinWindowsPerScale is fragile/fit-noise-prone, which is a
%     data-quantity limitation of DFA on short/gappy records, not
%     necessarily a code bug
%
% Verdict per file:
%   'fragile-fit'       - alpha2 fit range's thinnest scale is within 2x
%                         P.dfaMinWindowsPerScale (unstable fit, expected
%                         limitation, not a masking bug)
%   'short-recording'   - total valid RR duration < 120s (same root cause)
%   'UNEXPLAINED'       - neither of the above; worth a direct
%                         test_dfa_t15_crosscheck() run before trusting
%                         this file's alpha values -- this is the pattern
%                         that would indicate an actual masking bug
%
%   summary = investigate_implausible_alpha(datasetRoot)
%
% datasetRoot - root folder matching the batch's queue-derived paths,
%               e.g. 'G:\Shared drives\BIONICs Lab Workspace\Project Folders\GEMS\Survivals'
%               (default: that same path, as printed by run_dfa_batch)

if nargin < 1 || isempty(datasetRoot)
    datasetRoot = 'G:\Shared drives\BIONICs Lab Workspace\Project Folders\GEMS\Survivals';
end

P = pipeline_params();

% Relative paths of the 31 files flagged "IMPLAUSIBLE dfa_alpha" by the
% queue-scoped run_dfa_batch pass (see the run's printed summary).
relFiles = { ...
    '05062026\E1000_FRE_E1000_stim_rec_1406\E1000_FRE_E1000_stim_rec_1406_stim_HRVMeasures.mat'
    '05062026\E1000_JEL_E1000_bl_1315\E1000_JEL_E1000_bl_1315_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05062026\E1000_JEL_E1000_stim_rec_1325\E1000_JEL_E1000_stim_rec_1325_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05062026\E1000_JEL_E1000_stim_rec_1325\E1000_JEL_E1000_stim_rec_1325_recovery_HRVMeasures.mat'
    '05062026\E1000_JEL_E1000_stim_rec_1325\E1000_JEL_E1000_stim_rec_1325_stim_HRVMeasures.mat'
    '05082026\M100E10_JEL_CME2_stim_rec_2216\M100E10_JEL_CME2_stim_rec_2216_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05082026\M100E10_JEL_CME2_stim_rec_2216\M100E10_JEL_CME2_stim_rec_2216_recovery_HRVMeasures.mat'
    '05082026\M100E10_ORE_CME2_bl_2323\M100E10_ORE_CME2_bl_2323_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05082026\M100E10_ORE_CME2_stim_rec_2333\M100E10_ORE_CME2_stim_rec_2333_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05082026\M10E100_ORE_CME_stim_rec_1050\M10E100_ORE_CME_stim_rec_1050_recovery_HRVMeasures.mat'
    '05092026\M10E1000_Ore_CME_bl_1243\M10E1000_Ore_CME_bl_1243_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05092026\M10E1000_Ore_CME_stim_rec_1253\M10E1000_Ore_CME_stim_rec_1253_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05092026\M50E1000_LOL_CME_bl_2244\M50E1000_LOL_CME_bl_2244_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05092026\M50E1000_Ore_CME_bl_2120\M50E1000_Ore_CME_bl_2120_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05092026\M50E100_Ore_CME_stim_rec_1040\M50E100_Ore_CME_stim_rec_1040_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05102026\M100E1000_ORE_CME_bl_1147\M100E1000_ORE_CME_bl_1147_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05102026\M100E1000_ORE_CME_stim_rec_1157\M100E1000_ORE_CME_stim_rec_1157_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05112026\M10_ORE_M10_bl_1038\M10_ORE_M10_bl_1038_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05112026\M50_Lol_M50_bl_2202\M50_Lol_M50_bl_2202_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    '05112026\M50_Lol_M50_stim_rec_2212\M50_Lol_M50_stim_rec_2212_notched_v0.2.2_recovery_HRVMeasures.mat'
    '05112026\M50_ORE_M50_stim_rec_2045\M50_ORE_M50_stim_rec_2045_notched_v0.2.2_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05052026\E10_FRE_E10_bl_0020\E10_FRE_E10_bl_0020_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    'TDT\0505-0506-0507\05052026\E10_Jel_E10_bl_2255\E10_Jel_E10_bl_2255_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    'TDT\0505-0506-0507\05052026\E10_Jel_E10_stim_rec_2305\E10_Jel_E10_stim_rec_2305_notched_v0.2.2_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05062026\M10E10_FRE_M10E10_stim_rec_2001\M10E10_FRE_M10E10_stim_rec_2001_notched_v0.2.2_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05062026\M10E10_JEL_M10E10_stim_rec_1842\M10E10_JEL_M10E10_stim_rec_1842_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05062026\M50E10_JEL_M50E10_bl_2124\M50E10_JEL_M50E10_bl_2124_notched_v0.2.2_blankmotion_HRVMeasures.mat'
    'TDT\0505-0506-0507\05062026\M50E10_JEL_M50E10_stim_rec_2134\M50E10_JEL_M50E10_stim_rec_2134_notched_v0.2.2_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05062026\M50E10_JEL_M50E10_stim_rec_2134\M50E10_JEL_M50E10_stim_rec_2134_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05072026\M10E100_JEL_CME_stim_rec_1338\M10E100_JEL_CME_stim_rec_1338_recovery_HRVMeasures.mat'
    'TDT\0505-0506-0507\05072026\M10E100_lol_CME_stim_rec_1417\M10E100_lol_CME_stim_rec_1417_notched_v0.2.2_recovery_HRVMeasures.mat' ...
};

n = numel(relFiles);
summary = repmat(struct('file','','ok',false,'durationSec',NaN,'pctValid',NaN, ...
    'nRuns',NaN,'meanRunLenSec',NaN,'alpha1',NaN,'alpha2',NaN, ...
    'minWindowsInAlpha2Range',NaN,'verdict',''), n, 1);

for i = 1:n
    f = fullfile(datasetRoot, relFiles{i});
    summary(i).file = relFiles{i};
    if ~isfile(f)
        summary(i).verdict = 'FILE NOT FOUND';
        fprintf('[%2d/%2d] NOT FOUND: %s\n', i, n, f);
        continue;
    end
    try
        S = load(f, 'RR_intervals','RR_times','dfa_alpha1','dfa_alpha2');
        dfaOut = dfaRR_gapAware(S.RR_intervals, S.RR_times, P);

        runs = dfaOut.runs;
        runLenSec = arrayfun(@(r) sum(S.RR_intervals(runs(r,1):runs(r,2))), 1:size(runs,1));
        totalDur = S.RR_times(end) + S.RR_intervals(end) - S.RR_times(1);
        validDur = sum(runLenSec);

        inAlpha2Range = dfaOut.nVals > P.dfaLongRangeMinScale;
        if any(inAlpha2Range) && isfinite(dfaOut.alpha2)
            minW = min(dfaOut.nWindows(inAlpha2Range));
        else
            minW = NaN;
        end

        summary(i).ok = true;
        summary(i).durationSec = totalDur;
        summary(i).pctValid = 100*validDur/totalDur;
        summary(i).nRuns = size(runs,1);
        summary(i).meanRunLenSec = mean(runLenSec);
        summary(i).alpha1 = dfaOut.alpha1;
        summary(i).alpha2 = dfaOut.alpha2;
        summary(i).minWindowsInAlpha2Range = minW;

        if isfinite(minW) && minW <= 2*P.dfaMinWindowsPerScale
            summary(i).verdict = 'fragile-fit (alpha2 rests on scales barely above the window floor)';
        elseif validDur < 120
            summary(i).verdict = 'short-recording (< 120s valid RR data)';
        else
            summary(i).verdict = 'UNEXPLAINED -- worth a closer look';
        end

        fprintf(['[%2d/%2d] %-70s | dur=%6.1fs valid=%4.1f%% runs=%3d | a1=%.3f a2=%.3f | ' ...
                 'minW(a2 range)=%s | %s\n'], ...
            i, n, relFiles{i}, totalDur, summary(i).pctValid, summary(i).nRuns, ...
            summary(i).alpha1, summary(i).alpha2, num2str(minW), summary(i).verdict);
    catch ME
        summary(i).verdict = ['ERROR: ' ME.message];
        fprintf('[%2d/%2d] ERROR loading/processing %s: %s\n', i, n, relFiles{i}, ME.message);
    end
end

ok = [summary.ok];
fprintf('\n============ investigate_implausible_alpha summary ============\n');
fprintf('%d/%d files processed.\n', sum(ok), n);
if any(ok)
    verdicts = {summary(ok).verdict};
    cats = unique(verdicts);
    for c = 1:numel(cats)
        fprintf('  %3d file(s): %s\n', sum(strcmp(verdicts, cats{c})), cats{c});
    end
end
fprintf(['If most fall in fragile-fit/short-recording, this is a real DFA\n' ...
         'data-quantity limitation (expected per the literature-based caveat\n' ...
         'in DFA_IMPLEMENTATION_PLAN.md T5/T15), not a masking bug. Any\n' ...
         'UNEXPLAINED file deserves a direct test_dfa_t15_crosscheck() run\n' ...
         'before trusting its alpha values.\n']);
end
