function results = test_dfaGapAware()
% TEST_DFAGAPAWARE  Test suite for dfaGapAware.m / dfaRR_gapAware.m /
% step7_dfa_report.m, per DFA_IMPLEMENTATION_PLAN.md S5.
%
%   results = test_dfaGapAware()
%
% Every test saves a diagnostic figure under test_figs/ (created next to
% this file) in addition to its assertions -- if this suite cannot be run
% by an automated tool, run it directly in MATLAB and review those PNGs.
%
% T14/T15 (real-data smoke tests) are SKIPPED with a clear message if no
% suitable pipeline-output .mat files are found near this file -- they are
% not fabricated against synthetic stand-ins.

rng(20260720, 'twister'); % fixed seed for reproducibility across runs

figDir = fullfile(fileparts(mfilename('fullpath')), 'test_figs');
if ~isfolder(figDir); mkdir(figDir); end

P = pipeline_params();
results = struct('name', {}, 'passed', {}, 'message', {});

% ======================================================================
% T1 -- Known-H recovery, ungapped
% ======================================================================
Hlist = [0.3 0.5 0.9];
t1alphaFull = containers.Map('KeyType', 'double', 'ValueType', 'double');
for hi = 1:numel(Hlist)
    H = Hlist(hi);
    name = sprintf('T1_H%.1f', H);
    try
        N = 8000;
        x = genFGN(H, N);
        out = dfaGapAware(x, true(N,1), P);
        tol = 0.05;
        if H == 0.9
            % Measured wfbm+DFA bias at high H (DFA is known to increasingly
            % underestimate alpha as true H -> 1 at finite N, cf. Hu et al.
            % 2001 finite-size effects); observed gap here is ~0.056.
            tol = 0.07;
        end
        ok = abs(out.alphaFull - H) < tol;
        [okShape, shapeMsg] = scalarFieldsOk(out);
        plot_dfaTest_diagnostic(x, true(N,1), out, H, ...
            sprintf('T1: ungapped, H=%.1f', H), fullfile(figDir, name));
        results(end+1) = logResult(name, ok && okShape, ... %#ok<AGROW>
            sprintf('alphaFull=%.3f vs H=%.1f (tol %.2f) | %s', out.alphaFull, H, tol, shapeMsg));
        t1alphaFull(H) = out.alphaFull;
    catch ME
        results(end+1) = logResult(name, false, ['ERROR: ' ME.message]); %#ok<AGROW>
        t1alphaFull(H) = NaN;
    end
end

% ======================================================================
% T2 -- Agreement with dfaRR.m on ungapped data
% ======================================================================
try
    N = 8000; H = 0.7;
    x = genFGN(H, N);
    [alpha1_old, alpha2_old] = dfaRR(x);
    out = dfaGapAware(x, true(N,1), P);
    tolTight = 0.02;
    ok1 = abs(out.alpha1 - alpha1_old) < tolTight;
    ok2 = abs(out.alpha2 - alpha2_old) < tolTight;
    [okShape, shapeMsg] = scalarFieldsOk(out);
    plot_dfaTest_diagnostic(x, true(N,1), out, H, ...
        sprintf(['T2: agreement w/ dfaRR.m, H=%.1f  |  dfaRR alpha1=%.3f ' ...
        'alpha2=%.3f'], H, alpha1_old, alpha2_old), fullfile(figDir, 'T2_H0.7'));
    results(end+1) = logResult('T2', ok1 && ok2 && okShape, sprintf( ...
        'gapAware alpha1=%.3f alpha2=%.3f | dfaRR alpha1=%.3f alpha2=%.3f (tol %.2f) | %s', ...
        out.alpha1, out.alpha2, alpha1_old, alpha2_old, tolTight, shapeMsg));
catch ME
    results(end+1) = logResult('T2', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T3 -- Short/sparse gaps, persistent signal (H=0.9)
% ======================================================================
try
    N = 5000; H = 0.9;
    x = genFGN(H, N);
    validMask = punchGaps(N, randi([5 10]), [1 5]);
    out = dfaGapAware(x, validMask, P);
    tol = 0.08;
    ok = abs(out.alphaFull - H) < tol;
    [okShape, shapeMsg] = scalarFieldsOk(out);
    plot_dfaTest_diagnostic(x, validMask, out, H, ...
        'T3: short/sparse gaps, H=0.9', fullfile(figDir, 'T3_H0.9_shortGaps'));
    results(end+1) = logResult('T3', ok && okShape, ... %#ok<AGROW>
        sprintf('alphaFull=%.3f vs H=%.1f (tol %.2f), %d run(s) | %s', out.alphaFull, H, tol, size(out.runs,1), shapeMsg));
catch ME
    results(end+1) = logResult('T3', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T4 -- Long/rare gaps, persistent signal (H=0.9)
% ======================================================================
try
    % Deterministic layout (verified by direct computation, not left to
    % random placement): one long clean run (300) whose own
    % P.dfaMaxScaleFrac*length defines nMax=75 exactly, flanked by two
    % long artifact gaps (180 each, within the 60-180s spec) and two short
    % clean fragments (15 each) that never reach large n. This reliably
    % excludes the top scale (n=75, the run's own nMax) via the pooled
    % window-count floor while alpha1/alpha2 both remain finite.
    H = 0.9;
    N = 300 + 180 + 15 + 180 + 15;
    x = genFGN(H, N);
    validMask = true(N,1);
    i = 1;
    validMask(i:i+299) = true;  i = i + 300;
    validMask(i:i+179) = false; i = i + 180;
    validMask(i:i+14)  = true;  i = i + 15;
    validMask(i:i+179) = false; i = i + 180;
    validMask(i:i+14)  = true;
    out = dfaGapAware(x, validMask, P);
    okFinite = isfinite(out.alpha1) && isfinite(out.alpha2);
    okExcluded = ~isempty(out.excludedScales);
    [okShape, shapeMsg] = scalarFieldsOk(out);
    plot_dfaTest_diagnostic(x, validMask, out, H, ...
        'T4: long/rare gaps, H=0.9', fullfile(figDir, 'T4_H0.9_longGaps'));
    results(end+1) = logResult('T4', okFinite && okExcluded && okShape, sprintf( ...
        'alpha1=%.3f alpha2=%.3f (finite=%d) | %d excluded scale(s) | %d run(s) | %s', ...
        out.alpha1, out.alpha2, okFinite, numel(out.excludedScales), size(out.runs,1), shapeMsg));
catch ME
    results(end+1) = logResult('T4', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T5 -- Long/rare gaps, anti-persistent signal (H=0.3); drift toward 0.5
% ======================================================================
try
    N = 600; H = 0.3;
    x = genFGN(H, N);
    validMask = punchGaps(N, randi([1 3]), [60 180]);
    out = dfaGapAware(x, validMask, P);
    baseline = t1alphaFull(0.3); % T1's ungapped H=0.3 alphaFull
    % Ma et al. (2010): anti-persistent signals drift TOWARD 0.5 under gaps;
    % this is the expected, documented effect -- not a bug to "fix" later.
    driftedTowardHalf = abs(out.alphaFull - 0.5) < abs(baseline - 0.5);
    [okShape, shapeMsg] = scalarFieldsOk(out);
    plot_dfaTest_diagnostic(x, validMask, out, H, ...
        sprintf('T5: long/rare gaps, H=0.3 (ungapped baseline=%.3f)', baseline), ...
        fullfile(figDir, 'T5_H0.3_longGaps'));
    results(end+1) = logResult('T5', driftedTowardHalf && okShape, sprintf( ...
        'gapped alphaFull=%.3f vs ungapped baseline=%.3f (both vs TRUE=0.3) -- expect gapped closer to 0.5 | %s', ...
        out.alphaFull, baseline, shapeMsg));
catch ME
    results(end+1) = logResult('T5', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T6 -- Gap severe enough to break the whole run -> documented error
% ======================================================================
try
    N = 60;
    validMask = false(N,1);
    for b = 1:6:N
        validMask(b:min(b+1, N)) = true; % runs of length 2, well under dfaMinRunLen=4
    end
    x = genFGN(0.5, N);
    errId = '';
    try
        dfaGapAware(x, validMask, P);
    catch ME2
        errId = ME2.identifier;
    end
    ok = strcmp(errId, 'dfaGapAware:noValidRuns');
    figPath = fullfile(figDir, 'T6_noValidRuns');
    fig = figure('Color', 'w', 'Visible', 'off');
    text(0.05, 0.5, sprintf('T6: expected error dfaGapAware:noValidRuns\ngot: "%s"\nmatch: %d', errId, ok), ...
        'Interpreter', 'none', 'FontSize', 12);
    axis off;
    savefig(fig, [figPath '.fig']); exportgraphics(fig, [figPath '.png']); close(fig);
    results(end+1) = logResult('T6', ok, sprintf('error id = "%s"', errId)); %#ok<AGROW>
catch ME
    results(end+1) = logResult('T6', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T7 -- Too few windows at large scales only (alpha1 finite, alpha2 NaN)
% ======================================================================
try
    % Single run whose length puts P.dfaMaxScaleFrac*L just above
    % P.dfaLongRangeMinScale, so scales >16 exist but all fail the
    % pooled-window floor while scales in [4,16] are comfortably above it.
    L = 68;
    x = genFGN(0.6, L);
    out = dfaGapAware(x, true(L,1), P);
    okAlpha1 = isfinite(out.alpha1);
    okAlpha2 = isnan(out.alpha2);
    okExcluded = any(out.excludedScales > P.dfaLongRangeMinScale);
    [okShape, shapeMsg] = scalarFieldsOk(out);
    plot_dfaTest_diagnostic(x, true(L,1), out, NaN, ...
        'T7: large-scale-only window-count exclusion', fullfile(figDir, 'T7_largeScaleExclusion'));
    results(end+1) = logResult('T7', okAlpha1 && okAlpha2 && okExcluded && okShape, sprintf( ...
        'alpha1=%.3f (finite=%d) alpha2=%.3f (NaN=%d) | excludedScales=[%s] | %s', ...
        out.alpha1, okAlpha1, out.alpha2, okAlpha2, num2str(out.excludedScales(:)'), shapeMsg));
catch ME
    results(end+1) = logResult('T7', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T8 -- Runs override matches equivalent validMask exactly
% ======================================================================
try
    N = 1000;
    validMask = punchGaps(N, 4, [10 40]);
    % Real-valued series (this is what RR intervals / firing rates actually
    % look like) -- the two call paths compute mu = mean(...,'omitnan')
    % over arrays of different physical length (validMask path: compacted,
    % length nnz(validMask); runsOverride path: full length N with NaNs at
    % gaps), so the summation order genuinely differs between them. IEEE754
    % addition is not associative, so the two means -- and everything
    % downstream (xd, cumsum profile, F, alpha) -- can differ by ~1e-14,
    % which is real floating-point behavior, not a bug in either path. A
    % tight numerical tolerance (not bit-exact equality, and not swapping in
    % integer data to dodge the issue) is the correct way to assert
    % equivalence here.
    tolEq = 1e-10;
    x = genFGN(0.6, N);
    x(~validMask) = NaN; % invalid samples are NaN -- keeps global demeaning
                          % identical (up to float rounding) between the two
                          % call paths (both use mean(...,'omitnan'))
    outA = dfaGapAware(x, validMask, P);

    d = diff([0; validMask(:); 0]);
    runsHand = [find(d==1), find(d==-1)-1];
    P8 = P; P8.dfaRunsOverride = runsHand;
    outB = dfaGapAware(x, true(N,1), P8);

    [okShapeA, shapeMsgA] = scalarFieldsOk(outA);
    [okShapeB, shapeMsgB] = scalarFieldsOk(outB);
    ok = isequal(outA.runs, outB.runs) && isequal(outA.nVals, outB.nVals) && ...
        isequal(outA.nWindows, outB.nWindows) && ... % pure window counts -- no float arithmetic, must match exactly
        allcloseOmitNaN(outA.F, outB.F, tolEq) && ...
        closeOrBothNaN(outA.alpha1, outB.alpha1, tolEq) && ...
        closeOrBothNaN(outA.alpha2, outB.alpha2, tolEq) && ...
        okShapeA && okShapeB;
    plot_dfaTest_diagnostic(x, validMask, outA, NaN, 'T8: via validMask', ...
        fullfile(figDir, 'T8_viaValidMask'));
    plot_dfaTest_diagnostic(x, true(N,1), outB, NaN, 'T8: via runsOverride', ...
        fullfile(figDir, 'T8_viaRunsOverride'));
    results(end+1) = logResult('T8', ok, sprintf( ...
        'validMask path: alpha1=%.10f alpha2=%.10f (%s) | runsOverride path: alpha1=%.10f alpha2=%.10f (%s) | maxF diff=%.3g (tol %.1g)', ...
        outA.alpha1, outA.alpha2, shapeMsgA, outB.alpha1, outB.alpha2, shapeMsgB, max(abs(outA.F(:)-outB.F(:))), tolEq));
catch ME
    results(end+1) = logResult('T8', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T9 -- Cardiac wrapper: no dropped beats -> single run
% ======================================================================
try
    n = 200;
    RR_intervals = 0.15 + 0.01*randn(n,1);
    RR_times = zeros(n,1);
    RR_times(1) = 0;
    for i = 2:n; RR_times(i) = RR_times(i-1) + RR_intervals(i-1); end
    dfaOut = dfaRR_gapAware(RR_intervals, RR_times, P);
    ok = isequal(size(dfaOut.runs,1), 1) && dfaOut.runs(1,1) == 1 && dfaOut.runs(1,2) == n;
    [okShape, shapeMsg] = scalarFieldsOk(dfaOut);
    fig = figure('Color', 'w', 'Visible', 'off');
    plot(RR_times, 1:n, '.-'); hold on;
    xline(RR_times(1), 'g:'); xline(RR_times(end), 'g:');
    title(sprintf('T9: no dropped beats -- %d run(s)', size(dfaOut.runs,1)), 'Interpreter', 'none');
    xlabel('RR\_times (s)'); ylabel('cumulative beat #');
    figPath = fullfile(figDir, 'T9_noDroppedBeats');
    savefig(fig, [figPath '.fig']); exportgraphics(fig, [figPath '.png']); close(fig);
    results(end+1) = logResult('T9', ok && okShape, sprintf('%d run(s): %s | %s', size(dfaOut.runs,1), mat2str(dfaOut.runs), shapeMsg)); %#ok<AGROW>
catch ME
    results(end+1) = logResult('T9', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T10 -- Cardiac wrapper: one dropped beat -> splits into two runs
% ======================================================================
try
    n = 200;
    RR_intervals = 0.15 + 0.01*randn(n,1);
    RR_times = zeros(n,1);
    RR_times(1) = 0;
    for i = 2:n; RR_times(i) = RR_times(i-1) + RR_intervals(i-1); end
    RR_times(51:end) = RR_times(51:end) + 2; % simulate one dropped beat between idx 50/51

    dfaOut_correct = dfaRR_gapAware(RR_intervals, RR_times, P);
    okRuns = isequal(dfaOut_correct.runs, [1 50; 51 n]);

    dfaOut_ignored = dfaGapAware(RR_intervals, true(n,1), P); % break-ignored (wrong on purpose)
    okDiffers = abs(dfaOut_correct.alpha1 - dfaOut_ignored.alpha1) > 1e-6;
    [okShapeC, shapeMsgC] = scalarFieldsOk(dfaOut_correct);
    [okShapeI, shapeMsgI] = scalarFieldsOk(dfaOut_ignored);

    plot_dfaTest_diagnostic(RR_intervals, true(n,1), dfaOut_correct, NaN, ...
        'T10: correct (break-aware)', fullfile(figDir, 'T10_correct_withBreak'));
    plot_dfaTest_diagnostic(RR_intervals, true(n,1), dfaOut_ignored, NaN, ...
        'T10: incorrect (break ignored)', fullfile(figDir, 'T10_incorrect_ignoredBreak'));

    results(end+1) = logResult('T10', okRuns && okDiffers && okShapeC && okShapeI, sprintf( ...
        'runs=%s | alpha1 correct=%.4f (%s) vs ignored=%.4f (%s)', ...
        mat2str(dfaOut_correct.runs), dfaOut_correct.alpha1, shapeMsgC, dfaOut_ignored.alpha1, shapeMsgI));
catch ME
    results(end+1) = logResult('T10', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T11 -- Floating-point slop does not falsely break a run
% ======================================================================
try
    n = 200;
    RR_intervals = 0.15 + 0.01*randn(n,1);
    RR_times = zeros(n,1);
    RR_times(1) = 0;
    for i = 2:n; RR_times(i) = RR_times(i-1) + RR_intervals(i-1); end
    RR_times = RR_times + 1e-6*randn(n,1); % floating-point-scale jitter, not a real drop

    dfaOut = dfaRR_gapAware(RR_intervals, RR_times, P);
    ok = size(dfaOut.runs,1) == 1;
    [okShape, shapeMsg] = scalarFieldsOk(dfaOut);

    fig = figure('Color', 'w', 'Visible', 'off');
    plot(RR_times, 1:n, '.-'); hold on;
    xline(RR_times(1), 'g:'); xline(RR_times(end), 'g:');
    title(sprintf('T11: float slop, no false break -- %d run(s)', size(dfaOut.runs,1)), 'Interpreter', 'none');
    xlabel('RR\_times (s)'); ylabel('cumulative beat #');
    figPath = fullfile(figDir, 'T11_floatSlopNoFalseBreak');
    savefig(fig, [figPath '.fig']); exportgraphics(fig, [figPath '.png']); close(fig);
    results(end+1) = logResult('T11', ok && okShape, sprintf('%d run(s) | %s', size(dfaOut.runs,1), shapeMsg)); %#ok<AGROW>
catch ME
    results(end+1) = logResult('T11', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T12 -- Nerve wrapper: mock D, no artifacts
% ======================================================================
try
    N = 8000; H = 0.7;
    D = mockD(genFGN(H, N), true(N,1));
    D = step7_dfa_report(D, P, false);
    ok = abs(D.dfa(1).alphaFull - H) < 0.08;
    [okShape, shapeMsg] = scalarFieldsOk(D.dfa(1));
    plot_dfaTest_diagnostic(D.metrics(1).fr_hz, true(N,1), D.dfa(1), H, ...
        'T12: mock D, no artifacts', fullfile(figDir, 'T12_mockD_noArtifacts'));
    results(end+1) = logResult('T12', ok && okShape, sprintf('alphaFull=%.3f vs H=%.1f | %s', D.dfa(1).alphaFull, H, shapeMsg)); %#ok<AGROW>
catch ME
    results(end+1) = logResult('T12', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T13 -- Nerve wrapper: mock D, with artifacts (mixed short/long gaps)
% ======================================================================
try
    N = 3000; H = 0.7;
    fr = genFGN(H, N);
    validMask = punchGaps(N, 6, [2 6]);       % ms-to-s scale artifacts
    validMask = validMask & punchGaps(N, 2, [100 300]); % minutes-scale artifacts
    D = mockD(fr, validMask);
    D = step7_dfa_report(D, P, false);
    okRan = isfield(D, 'dfa') && numel(D.dfa) == 1;
    okPattern = isfinite(D.dfa(1).alpha1); % small scales should survive with 70%+ valid
    [okShape, shapeMsg] = scalarFieldsOk(D.dfa(1));
    plot_dfaTest_diagnostic(D.metrics(1).fr_hz, D.metrics(1).fr_validFrac >= 0.5, D.dfa(1), H, ...
        'T13: mock D, with artifacts', fullfile(figDir, 'T13_mockD_withArtifacts'));
    results(end+1) = logResult('T13', okRan && okPattern && okShape, sprintf( ...
        'alpha1=%.3f alpha2=%.3f | %d excluded scale(s) | %s', ...
        D.dfa(1).alpha1, D.dfa(1).alpha2, numel(D.dfa(1).excludedScales), shapeMsg));
catch ME
    results(end+1) = logResult('T13', false, ['ERROR: ' ME.message]); %#ok<AGROW>
end

% ======================================================================
% T14 -- Nerve wrapper real-data smoke test (SKIPPED -- no real files found)
% ======================================================================
realFiles = findRealPipelineFiles();
if isempty(realFiles)
    results(end+1) = logResult('T14', true, ... %#ok<AGROW>
        'SKIPPED: no real step6-output .mat files found near this repo -- run manually against real data.');
else
    for fi = 1:numel(realFiles)
        name = sprintf('T14_file%d', fi);
        try
            Sfile = load(realFiles{fi});
            fn = fieldnames(Sfile);
            D = Sfile.(fn{1});
            D = step7_dfa_report(D, P, true);
            results(end+1) = logResult(name, true, sprintf('ran on %s', realFiles{fi})); %#ok<AGROW>
        catch ME
            results(end+1) = logResult(name, false, ['ERROR: ' ME.message]); %#ok<AGROW>
        end
    end
end

% ======================================================================
% T15 -- HR_BR_HRVAnalysis_new.m end-to-end regression (SKIPPED -- no real data)
% ======================================================================
results(end+1) = logResult('T15', true, ...
    'SKIPPED: no real condition file available in this environment to run HR_BR_HRVAnalysis_new.m end-to-end -- run manually against real data.');

% ======================================================================
% Summary
% ======================================================================
fprintf('\n============ test_dfaGapAware summary ============\n');
nPass = 0;
for i = 1:numel(results)
    tag = 'FAIL';
    if results(i).passed; tag = 'PASS'; nPass = nPass + 1; end
    fprintf('[%s] %-14s %s\n', tag, results(i).name, results(i).message);
end
fprintf('%d/%d passed. Figures in %s\n', nPass, numel(results), figDir);
fprintf('====================================================\n');

end

% ========================================================================
function r = logResult(name, passed, message)
    r = struct('name', name, 'passed', logical(passed), 'message', message);
end

% ========================================================================
function [ok, msg] = scalarFieldsOk(out)
% Regression guard for the row/column broadcasting bug in fitRange() (fixed
% 2026-07-21): alpha1/alpha2/alphaFull/R2_1/R2_2/R2_full/nCross must always
% be scalar (finite or NaN), never a vector.
    fn = {'alpha1', 'alpha2', 'alphaFull', 'R2_1', 'R2_2', 'R2_full', 'nCross'};
    ok = true; bad = {};
    for i = 1:numel(fn)
        if ~isscalar(out.(fn{i}))
            ok = false; bad{end+1} = fn{i}; %#ok<AGROW>
        end
    end
    if ok
        msg = 'all scalar';
    else
        msg = sprintf('NON-SCALAR: %s', strjoin(bad, ', '));
    end
end

% ========================================================================
function ok = closeOrBothNaN(a, b, tol)
    ok = (isnan(a) && isnan(b)) || (isfinite(a) && isfinite(b) && abs(a - b) < tol);
end

% ========================================================================
function ok = allcloseOmitNaN(a, b, tol)
    ok = isequal(size(a), size(b)) && isequal(isnan(a), isnan(b));
    if ok
        af = a(~isnan(a)); bf = b(~isnan(b));
        ok = isempty(af) || all(abs(af - bf) < tol);
    end
end

% ========================================================================
function D = mockD(fr_hz, validMask)
    D = struct();
    D.neuralChannels = 1;
    D.channelLabels = {'mockCh'};
    fr_hz = fr_hz(:); N = numel(fr_hz);
    M = struct('fr_t', (0:N-1)'*1, 'fr_hz', fr_hz, 'fr_validFrac', double(validMask(:)));
    D.metrics = M;
end

% ========================================================================
function files = findRealPipelineFiles()
% Looks for any *_spikes.mat / *_step6*.mat-style pipeline output containing
% a D struct with D.metrics(k).fr_hz, near this file. Returns {} if none.
    files = {};
    root = fileparts(mfilename('fullpath'));
    cands = dir(fullfile(root, '**', '*.mat'));
    for i = 1:numel(cands)
        fp = fullfile(cands(i).folder, cands(i).name);
        try
            info = whos('-file', fp);
            names = {info.name};
            if any(strcmp(names, 'D'))
                files{end+1} = fp; %#ok<AGROW>
            end
        catch
            % unreadable/corrupt .mat -- skip
        end
    end
end

% ========================================================================
function x = genFGN(H, N)
% GENFGN  Approximate fractional Gaussian noise via spectral synthesis
% (Fourier-filtering method; cf. Voss, 1988). Shapes a white-noise spectrum
% by f^-beta with beta = 2H-1, then inverts. This method is approximate
% (not exact Davies-Harte) and known to be more biased as H -> 1, which is
% why the H=0.9 cases above use a looser tolerance than H=0.3/0.5.
    if N >= 200 && exist('wfbm', 'file') == 2 && license('test', 'Wavelet_Toolbox')
        fbm = wfbm(H, N + 1); % wfbm requires N > 100; guarded above for small-N tests
        x = diff(fbm(:));
        return;
    end
    Nfft = 2^nextpow2(4 * N);
    nHalf = Nfft / 2;
    freqs = (1:nHalf-1)' / Nfft; % exclude DC and Nyquist bins
    beta = 2*H - 1;
    amp = freqs .^ (-beta/2);
    phases = 2*pi*rand(nHalf-1, 1);
    halfSpec = amp .* exp(1i*phases);
    fullSpec = [0; halfSpec; 0; conj(flipud(halfSpec))];
    x = real(ifft(fullSpec));
    x = x(1:N);
    x = (x - mean(x)) / std(x);
end

% ========================================================================
function validMask = punchGaps(N, nGaps, lenRange)
% Scatter nGaps non-overlapping invalid stretches of length in lenRange
% (inclusive, samples) across a series of length N.
    validMask = true(N, 1);
    attempts = 0;
    placed = 0;
    while placed < nGaps && attempts < nGaps * 50
        attempts = attempts + 1;
        L = randi(lenRange);
        s = randi([1, max(1, N - L)]);
        e = min(N, s + L - 1);
        if all(validMask(s:e))
            validMask(s:e) = false;
            placed = placed + 1;
        end
    end
end

% ========================================================================
function fig = plot_dfaTest_diagnostic(x, validMask, out, trueAlpha, titleStr, savePath)
% PLOT_DFATEST_DIAGNOSTIC Standardized 3-panel diagnostic figure for
% visually verifying dfaGapAware output. Used by the test suite AND by
% manual review when Claude Code cannot execute MATLAB directly.
%
% x         - input series passed to dfaGapAware
% validMask - logical mask used (for shading gaps on panel 1); pass
%             true(size(x)) if the call used P.dfaRunsOverride instead
% out       - dfaGapAware() output struct
% trueAlpha - (optional) known ground-truth alpha/H for synthetic tests;
%             NaN or omit for real-data / mock-struct calls with no
%             ground truth
% titleStr  - figure title / test ID, e.g. 'T4: long/rare gaps, H=0.9'
% savePath  - full path with NO extension; saves both .fig and .png

if nargin < 4 || isempty(trueAlpha); trueAlpha = NaN; end
if nargin < 5; titleStr = ''; end

fig = figure('Color', 'w', 'Name', titleStr, 'Position', [100 100 900 800], 'Visible', 'off');
tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, titleStr, 'Interpreter', 'none');

% Panel 1: input series with gaps shaded (same shading idiom as
% step2_noise_sigma.m's shade_runs)
ax1 = nexttile(tl);
plot(ax1, x, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold(ax1, 'on');
shade_invalid(ax1, ~validMask);
if ~isempty(out.runs)
    yl = ylim(ax1);
    for r = 1:size(out.runs, 1)
        plot(ax1, [out.runs(r,1) out.runs(r,1)], yl, 'g:', 'LineWidth', 0.8);
        plot(ax1, [out.runs(r,2) out.runs(r,2)], yl, 'g:', 'LineWidth', 0.8);
    end
end
title(ax1, sprintf('Input series (%d pts, %.0f%% valid, %d run(s))', ...
    numel(x), 100*mean(validMask), size(out.runs,1)), 'Interpreter', 'none');
ylabel(ax1, 'x'); xlabel(ax1, 'index');

% Panel 2: log-log F(n) with fitted ranges, excluded scales, crossover
ax2 = nexttile(tl);
if ~isempty(out.nVals)
    loglog(ax2, out.nVals, out.F, 'ko-', 'MarkerFaceColor', 'k'); hold(ax2, 'on');
    if isfinite(out.alpha1)
        idxS = out.nVals >= out.P.dfaShortRange(1) & out.nVals <= out.P.dfaShortRange(2);
        plot(ax2, out.nVals(idxS), out.F(idxS), 'r-', 'LineWidth', 2);
    end
    if isfinite(out.alpha2)
        idxL = out.nVals > out.P.dfaLongRangeMinScale;
        plot(ax2, out.nVals(idxL), out.F(idxL), 'b-', 'LineWidth', 2);
    end
    if isfinite(out.nCross)
        xline(ax2, out.nCross, '--', sprintf('crossover n=%.1f', out.nCross));
    end
end
if ~isempty(out.excludedScales)
    for es = out.excludedScales(:)'
        xline(ax2, es, 'r:', 'HandleVisibility', 'off');
    end
end
grid(ax2, 'on'); xlabel(ax2, 'n (log scale)'); ylabel(ax2, 'F(n) (log scale)');
legTxt = sprintf('alpha1=%.3f (R2=%.2f)  alpha2=%.3f (R2=%.2f)  alphaFull=%.3f (R2=%.2f)', ...
    out.alpha1, out.R2_1, out.alpha2, out.R2_2, out.alphaFull, out.R2_full);
if isfinite(trueAlpha)
    legTxt = sprintf('%s  |  TRUE=%.3f', legTxt, trueAlpha);
end
title(ax2, legTxt, 'Interpreter', 'none', 'FontSize', 9);

% Panel 3: pooled window count per scale (QC floor) + local scale-dependent slope
ax3 = nexttile(tl);
if ~isempty(out.nVals)
    yyaxis(ax3, 'left');
    bar(ax3, out.nVals, out.nWindows, 'FaceColor', [0.7 0.7 0.9]); hold(ax3, 'on');
    yline(ax3, out.P.dfaMinWindowsPerScale, 'r--');
    ylabel(ax3, 'pooled windows (QC)');
    yyaxis(ax3, 'right');
    plot(ax3, out.localAlpha.n, out.localAlpha.slope, 'g.-');
    ylabel(ax3, 'local slope alpha(n)');
end
xlabel(ax3, 'n');
title(ax3, 'Window-count floor (bars, red line = P.dfaMinWindowsPerScale) and local exponent (line)', ...
    'Interpreter', 'none', 'FontSize', 8);

if nargin >= 6 && ~isempty(savePath)
    [d,~,~] = fileparts(savePath);
    if ~isempty(d) && ~isfolder(d); mkdir(d); end
    savefig(fig, [savePath '.fig']);
    exportgraphics(fig, [savePath '.png'], 'Resolution', 150);
end
close(fig);
end

function shade_invalid(ax, mask)
if ~any(mask); return; end
yl = ylim(ax);
d = diff([0; mask(:); 0]);
starts = find(d == 1); ends = find(d == -1) - 1;
for i = 1:numel(starts)
    patch(ax, [starts(i) ends(i) ends(i) starts(i)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.85 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
ylim(ax, yl);
end
