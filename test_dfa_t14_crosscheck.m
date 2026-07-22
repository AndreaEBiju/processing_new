function report = test_dfa_t14_crosscheck(neuralFile, hrbrFile, P)
% TEST_DFA_T14_CROSSCHECK  Real-data cross-check for step7_dfa_report, per
% DFA_IMPLEMENTATION_PLAN.md T14: confirm that the gap structure DFA acts
% on (step7's binValidMask / dfaOut.runs) is actually explained by the
% file's own artifact annotations -- D.removedSegmentIdx (the segments
% browseMotionArtifacts.m produced) and D.cardiacBlank (the automatic
% +/-15ms R-peak blank) -- rather than some unit/format mismatch between
% what was annotated and what feeds the DFA-level mask.
%
%   report = test_dfa_t14_crosscheck(neuralFile, hrbrFile, P)
%
% neuralFile - path to a real '*_blankmotion.mat' / '*_recovery_blankmotion.mat'
% hrbrFile   - matching '*_HRBR.mat' sibling (R-peak locations)
% P          - (optional) pipeline_params()
%
% Example (path as printed by run_dfa_batch on the real GEMS dataset --
% chosen because the console log showed 340.3s flagged-flat / only 32%
% valid on this file, a heavy-artifact case worth checking):
%   root = 'G:\Shared drives\BIONICs Lab Workspace\Project Folders\GEMS\Survivals\05252026\M50E100_lol_CME1_bl_2124';
%   report = test_dfa_t14_crosscheck( ...
%       fullfile(root, 'M50E100_lol_CME1_bl_2124_notched_v0.2.2_blankmotion.mat'), ...
%       fullfile(root, 'M50E100_lol_CME1_bl_2124_notched_v0.2.2_blankmotion_HRBR.mat'));
%
% report(k) per channel: .label, .nGaps (DFA run breaks), .nExplainedByMotion,
% .nExplainedByCardiacOnly, .nUnexplained (should be 0 -- investigate the
% saved figure if not, since it means step7's gap detection is disagreeing
% with the actual annotated/known blanking on this file).
% Saves a comparison figure per channel to test_figs/T14_crosscheck_<label>.png.

if nargin < 3 || isempty(P); P = pipeline_params(); end

loadOpts = struct('neuralCols',[1 2], 'labels',{{'RVN','LVN'}}, 'rVar','heartlocs');
D = bulk_load_one(neuralFile, hrbrFile, loadOpts);
D = process_dataset(D, P, false);
D = step7_dfa_report(D, P, false);

figDir = fullfile(fileparts(mfilename('fullpath')), 'test_figs');
if ~isfolder(figDir); mkdir(figDir); end

fs = D.fs;
motionRanges  = idx_to_ranges(D.removedSegmentIdx, fs);      % human-annotated (browseMotionArtifacts.m)
cardiacRanges = mask_to_ranges_fs(D.cardiacBlank, fs);       % automatic R-peak blank

nCh = numel(D.neuralChannels);
report = repmat(struct('label','','nGaps',0,'nExplainedByMotion',0, ...
    'nExplainedByCardiacOnly',0,'nUnexplained',0), 1, nCh);

for k = 1:nCh
    M = D.metrics(k);
    dfaOut = D.dfa(k);
    t = M.fr_t(:);
    runs = dfaOut.runs;

    nGaps = max(0, size(runs,1) - 1);
    explainedMotion = false(nGaps,1);
    explainedCardiacOnly = false(nGaps,1);
    for g = 1:nGaps
        g0 = t(min(runs(g,2),   numel(t)));
        g1 = t(min(runs(g+1,1), numel(t)));
        explainedMotion(g) = range_overlaps(motionRanges, g0, g1);
        if ~explainedMotion(g)
            explainedCardiacOnly(g) = range_overlaps(cardiacRanges, g0, g1);
        end
    end
    unexplained = ~explainedMotion & ~explainedCardiacOnly;

    report(k).label = D.channelLabels{k};
    report(k).nGaps = nGaps;
    report(k).nExplainedByMotion = sum(explainedMotion);
    report(k).nExplainedByCardiacOnly = sum(explainedCardiacOnly);
    report(k).nUnexplained = sum(unexplained);

    fprintf(['[T14 crosscheck] %s: %d DFA gap(s) -- %d motion-explained, ' ...
             '%d cardiac-only, %d UNEXPLAINED\n'], report(k).label, nGaps, ...
        report(k).nExplainedByMotion, report(k).nExplainedByCardiacOnly, report(k).nUnexplained);
    if report(k).nUnexplained > 0
        warning('t14_crosscheck:unexplained', ...
            '%s: %d gap(s) not explained by removedSegmentIdx or cardiacBlank -- possible mismatch, inspect figure.', ...
            report(k).label, report(k).nUnexplained);
    end

    fig = figure('Color','w','Name',sprintf('T14 crosscheck -- %s', report(k).label), ...
        'Position',[120 80 1100 700]);
    tl = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
    title(tl, sprintf('T14 real-data crosscheck: %s', report(k).label), 'Interpreter','none');

    ax1 = nexttile(tl);
    plot(ax1, D.t, D.filtered(:,k), 'Color',[0.3 0.3 0.3],'LineWidth',0.3); hold(ax1,'on');
    shade_ranges(ax1, motionRanges,  [0.85 0.3 0.3], 0.25); % red: removedSegmentIdx (browseMotionArtifacts)
    shade_ranges(ax1, cardiacRanges, [0.3 0.3 0.85], 0.15); % blue: cardiacBlank
    title(ax1, 'Annotation sources: removedSegmentIdx (red) vs cardiacBlank (blue)', 'Interpreter','none');
    xlabel(ax1,'Time (s)'); ylabel(ax1,'filtered (V)');

    ax2 = nexttile(tl);
    plot(ax2, t, M.fr_hz, 'k','LineWidth',0.5); hold(ax2,'on');
    ylim(ax2, ylim(ax2)); % lock limits before shading gap patches
    for g = 1:nGaps
        g0 = t(min(runs(g,2),   numel(t)));
        g1 = t(min(runs(g+1,1), numel(t)));
        col = [0.6 0.9 0.6]; if unexplained(g); col = [1 0.6 0]; end
        yl = ylim(ax2);
        patch(ax2,[g0 g1 g1 g0],[yl(1) yl(1) yl(2) yl(2)],col,'FaceAlpha',0.4,'EdgeColor','none');
    end
    title(ax2, 'DFA gaps: green=explained, orange=UNEXPLAINED', 'Interpreter','none');
    xlabel(ax2,'Time (s)'); ylabel(ax2,'firing rate (Hz)');
    linkaxes([ax1 ax2],'x');

    savePath = fullfile(figDir, sprintf('T14_crosscheck_%s', report(k).label));
    savefig(fig, [savePath '.fig']);
    exportgraphics(fig, [savePath '.png'], 'Resolution',150);
end
end

% ========================================================================
function ranges = idx_to_ranges(segIdx, fs)
    if isempty(segIdx); ranges = zeros(0,2); return; end
    ranges = sort(segIdx,1) / fs;
end

function ranges = mask_to_ranges_fs(mask, fs)
    if ~any(mask); ranges = zeros(0,2); return; end
    d = diff([0; mask(:); 0]);
    s = find(d==1); e = find(d==-1) - 1;
    ranges = [s(:) e(:)] / fs;
end

function tf = range_overlaps(ranges, g0, g1)
    tf = false;
    if isempty(ranges); return; end
    tf = any(ranges(:,1) <= g1 & ranges(:,2) >= g0);
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
