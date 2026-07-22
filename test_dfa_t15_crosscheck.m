function report = test_dfa_t15_crosscheck(hrvFile)
% TEST_DFA_T15_CROSSCHECK  Real-data cross-check for HR_BR_HRVAnalysis_new.m's
% cardiac DFA, per DFA_IMPLEMENTATION_PLAN.md T15: confirm that every
% RR-level run break dfaRR_gapAware infers (a "beat was dropped here")
% actually lines up with a real invalidMask stretch from the SAME file's
% sample-level ground truth, rather than the two disagreeing.
%
%   report = test_dfa_t15_crosscheck(hrvFile)
%
% hrvFile - path to a real '*_HRVMeasures.mat' (must contain RR_intervals,
%           RR_times, invalidMask, t, dfa_alpha1 -- all saved by
%           HR_BR_HRVAnalysis_new.m)
%
% Example (path as printed by run_dfa_batch's real-dataset run -- one of
% the files flagged with a wide alpha1/alpha2 split, alpha1=0.291 alpha2=1.063):
%   hrvFile = 'G:\Shared drives\BIONICs Lab Workspace\Project Folders\GEMS\Survivals\TDT\0505-0506-0507\05062026\M10E10_JEL_M10E10_stim_rec_1842\M10E10_JEL_M10E10_stim_rec_1842_recovery_HRVMeasures.mat';
%   report = test_dfa_t15_crosscheck(hrvFile);
%
% report.alphaRecomputeDiff - recomputed vs saved dfa_alpha1 (sanity check
%   that the saved scalar truly came from the RR_intervals/RR_times stored
%   alongside it, not a stale value from an earlier run of this file).
% report.nBreaks/.nExplained/.nUnexplained - RR-run breaks vs invalidMask
%   overlap (nUnexplained should be 0; investigate the saved figure if not).
% Saves a comparison figure to test_figs/T15_crosscheck_<name>.png.

S = load(hrvFile, 'RR_intervals','RR_times','invalidMask','t','dfa_alpha1','dfa_alpha2');
if ~isfield(S,'RR_intervals') || numel(S.RR_intervals) < 2
    error('test_dfa_t15_crosscheck:noRR', 'No usable RR_intervals in %s', hrvFile);
end

P = pipeline_params();
dfaOut = dfaRR_gapAware(S.RR_intervals, S.RR_times, P);

alphaDiff = abs(dfaOut.alpha1 - S.dfa_alpha1);
matchTag = 'MATCH'; if ~(alphaDiff < 1e-9); matchTag = 'MISMATCH -- investigate'; end
fprintf('[T15 crosscheck] recomputed alpha1=%.4f vs saved=%.4f (diff=%.2e) -- %s\n', ...
    dfaOut.alpha1, S.dfa_alpha1, alphaDiff, matchTag);

invRanges = mask_to_ranges(S.invalidMask, S.t);

runs = dfaOut.runs;
nBreaks = max(0, size(runs,1) - 1);
explained = false(nBreaks,1);
for r = 1:nBreaks
    g0 = S.RR_times(runs(r,2));
    g1 = S.RR_times(runs(r+1,1));
    explained(r) = ~isempty(invRanges) && any(invRanges(:,1) <= g1 & invRanges(:,2) >= g0);
end
nExplained = sum(explained);
nUnexplained = nBreaks - nExplained;
if nUnexplained > 0
    fprintf('[T15 crosscheck] %d/%d RR-run break(s) overlap invalidMask -- %d UNEXPLAINED\n', ...
        nExplained, nBreaks, nUnexplained);
    warning('t15_crosscheck:unexplained', ...
        '%d RR-run break(s) do not overlap any invalidMask stretch in %s -- investigate.', ...
        nUnexplained, hrvFile);
else
    fprintf('[T15 crosscheck] %d/%d RR-run break(s) overlap invalidMask (all explained)\n', ...
        nExplained, nBreaks);
end

report = struct('file', hrvFile, 'alphaRecomputeDiff', alphaDiff, ...
    'nBreaks', nBreaks, 'nExplained', nExplained, 'nUnexplained', nUnexplained);

% ---- figure ----
[~, base] = fileparts(hrvFile);
figDir = fullfile(fileparts(mfilename('fullpath')), 'test_figs');
if ~isfolder(figDir); mkdir(figDir); end
fig = figure('Color','w','Name',sprintf('T15 crosscheck -- %s', base),'Position',[120 80 1100 650]);
tl = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
title(tl, sprintf('T15 real-data crosscheck: %s', base), 'Interpreter','none');

ax1 = nexttile(tl);
ylim(ax1, [-1 1]); hold(ax1,'on');
shade_ranges(ax1, invRanges, [0.85 0.3 0.3], 0.3);
title(ax1, sprintf('Sample-level invalidMask (%.1f%% invalid)', 100*mean(S.invalidMask)), 'Interpreter','none');
xlabel(ax1,'Time (s)'); ax1.YTick = [];
xlim(ax1, [S.t(1) S.t(end)]);

ax2 = nexttile(tl);
plot(ax2, S.RR_times, S.RR_intervals*1000, 'k.-'); hold(ax2,'on');
for r = 1:nBreaks
    g0 = S.RR_times(runs(r,2)); g1 = S.RR_times(runs(r+1,1));
    col = [0.6 0.9 0.6]; if ~explained(r); col = [1 0.6 0]; end
    yl = ylim(ax2);
    patch(ax2,[g0 g1 g1 g0],[yl(1) yl(1) yl(2) yl(2)],col,'FaceAlpha',0.4,'EdgeColor','none');
    ylim(ax2, yl);
end
title(ax2, 'RR-run breaks: green=explained by invalidMask, orange=UNEXPLAINED', 'Interpreter','none');
xlabel(ax2,'Time (s)'); ylabel(ax2,'RR interval (ms)');
linkaxes([ax1 ax2],'x');

savePath = fullfile(figDir, sprintf('T15_crosscheck_%s', base));
savefig(fig, [savePath '.fig']);
exportgraphics(fig, [savePath '.png'], 'Resolution',150);
end

% ========================================================================
function ranges = mask_to_ranges(mask, t)
    if ~any(mask); ranges = zeros(0,2); return; end
    d = diff([0; mask(:); 0]);
    s = find(d==1); e = find(d==-1) - 1;
    e = min(e, numel(t));
    ranges = [t(s) t(e)];
end

function shade_ranges(ax, ranges, color, alpha)
    if isempty(ranges); return; end
    yl = ylim(ax);
    for i = 1:size(ranges,1)
        patch(ax,[ranges(i,1) ranges(i,2) ranges(i,2) ranges(i,1)],[yl(1) yl(1) yl(2) yl(2)], ...
            color,'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
    end
    ylim(ax, yl);
end
