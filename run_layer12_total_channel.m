function results = run_layer12_total_channel(rawRows, saveDir)
%RUN_LAYER12_TOTAL_CHANNEL Layer 1-2 fit to TOTAL (RVN+LVN summed) firing
%   rate, cross-checked against bulk_mixed_models' existing "Total firing
%   rate (LVN+RVN)" metric -- an ALTERNATIVE to the per-channel fit in
%   run_layer12_first_pass.m, run in PARALLEL (does not replace it; that
%   driver is unmodified and still the per-channel deliverable).
%
%   results = run_layer12_total_channel(rawRows)
%   results = run_layer12_total_channel(rawRows, saveDir)
%
%   rawRows must be the SAME per-channel (RVN/LVN) shape run_layer12_first_pass.m
%   takes -- this driver sums each matched trial's RVN+LVN raw mean rate
%   internally (mirroring bulk_mixed_models.m's build_total_rate_table
%   grouping convention exactly: match by stem+phase, fallback to
%   animal+condition+phase) rather than requiring a separately-prepared input.
%
%   WHY: bulk_mixed_models' cross-check already shows "Total firing rate
%   (LVN+RVN)" has one of the strongest M x E interaction signals of any
%   metric -- this fits the SAME r_vagus model (model_layer12_equations.m,
%   unchanged) to that combined signal directly, instead of separately to
%   RVN and LVN.
%
%   IDENTIFIABILITY CAVEAT -- READ BEFORE TRUSTING RAW VALUES: with only ONE
%   channel here, Stage 1 is back to a plain 3-parameter fit
%   [W_c_TOTAL, u_half, f_max_c] against just 3 E-alone data points. There is
%   no second channel left to pool against, so the exact (u_half, f_max_c)
%   swap symmetry PROVEN in model_fitting_plan.md L12 is NOT resolved by
%   pooling here (unlike run_layer12_first_pass.m's per-channel path, which
%   pools RVN+LVN specifically to fix this). This driver instead runs the
%   same L13-style dual-basin Stage 2 check validated useful for the
%   per-channel case: fit Stage 2 with both the normal and the
%   swapped-and-rescaled Stage 1 solution, and let Stage 2's u_M dependence
%   pick a winner (report which, and how strongly).

    if nargin < 2; saveDir = ''; end

    % ---- Step 1: unit tests gate (equations file unchanged -- same gate) --
    fprintf('==== Step 1: unit tests (test_layer12_equations) ====\n');
    testResults = runtests('test_layer12_equations');
    nFailed = nnz(~[testResults.Passed]);
    if nFailed > 0
        disp(testResults);
        error('run_layer12_total_channel:testsFailed', ...
            '%d/%d unit tests FAILED -- aborting before any fitting.', nFailed, numel(testResults));
    end
    fprintf('All %d unit tests passed.\n\n', numel(testResults));

    % ---- Step 2: build TOTAL pseudo-channel rows (sum RVN+LVN raw rate) ---
    totalRows = build_total_rows_local(rawRows);
    fprintf(['[total] %d TOTAL rows built (paired RVN+LVN sums, per bulk_mixed_models.m''s ' ...
        'build_total_rate_table grouping convention).\n\n'], numel(totalRows));

    % ---- Step 3: bulk_mixed_models cross-check -----------------------------
    % Uses the ORIGINAL per-channel rawRows, NOT totalRows -- bulk_mixed_models
    % computes "Total firing rate (LVN+RVN)" itself internally (from out.raw),
    % so it needs the per-channel rows to do that, same as run_layer12_first_pass.m.
    fprintf('==== bulk_mixed_models cross-check (for "Total firing rate (LVN+RVN)" synergy) ====\n');
    R = [];
    if exist('bulk_mixed_models','file')==2 && exist('bulk_compile','file')==2
        try
            normRows = bulk_compile(rawRows);
            out = struct('raw', rawRows, 'norm', normRows);
            R = bulk_mixed_models(out, saveDir, [], false);   % includeSystemic=false, see run_layer12_first_pass.m
        catch ME
            warning('run_layer12_total_channel:bmm', 'bulk_mixed_models cross-check failed (%s) -- proceeding without it.', ME.message);
        end
    else
        warning('run_layer12_total_channel:bmmMissing', 'bulk_mixed_models/bulk_compile not on path -- proceeding without the cross-check.');
    end
    fprintf('\n');

    % ---- Step 4: Stage 1, single-group (NO pooling possible -- one channel) --
    fprintf('==== Stage 1: TOTAL E-alone fit (3 params, single group) ====\n');
    fitRaw = fit_stage1_single_local(totalRows);
    fprintf('\n');

    % ---- Step 5: Stage 2 for BOTH candidates (L15 applied to the single-
    % channel case): the swap symmetry (L12) is unresolved by pooling here (no
    % second channel), so build both candidates exactly as
    % fit_layer12_stage1_electrical.m does, fit Stage 2 for both, and SELECT
    % the winner by margin -- do not just report the "normal" one as if it
    % were the answer (that was the pre-L15 behavior).
    rescaleFactor = fitRaw.f_max_c / fitRaw.u_half;
    candidates(1) = struct('u_half',fitRaw.u_half, 'f_max_c',fitRaw.f_max_c, ...
        'W_c_TOTAL',fitRaw.W_c_TOTAL, 'label','raw');
    candidates(2) = struct('u_half',fitRaw.f_max_c, 'f_max_c',fitRaw.u_half, ...
        'W_c_TOTAL',fitRaw.W_c_TOTAL*rescaleFactor, 'label','swapped');

    fprintf('==== Stage 2: TOTAL interaction fit, both candidates ====\n');
    thetaFullCand = cell(1,2); resnorm2Cand = nan(1,2); stage2ErrCand = cell(1,2);
    for ci = 1:2
        try
            [thetaFullCand{ci}, resnorm2Cand(ci)] = fit_layer12_stage2_interaction(totalRows, 'TOTAL', candidates(ci));
        catch ME
            warning('run_layer12_total_channel:stage2', 'candidate %d (%s) fit failed (%s).', ci, candidates(ci).label, ME.message);
            stage2ErrCand{ci} = ME.message;
        end
    end
    fprintf('\n');

    % ---- Step 6: select the winning candidate (L15, single-channel case) --
    MIN_MARGIN_PCT = 10;
    resnorm2ForSelection = resnorm2Cand; resnorm2ForSelection(isnan(resnorm2ForSelection)) = Inf;
    [bestResnorm2, winnerIdx] = min(resnorm2ForSelection);
    otherIdx = 3 - winnerIdx;
    if isfinite(resnorm2ForSelection(otherIdx)) && isfinite(bestResnorm2) && bestResnorm2 > 0
        marginPct = 100 * (resnorm2ForSelection(otherIdx) - bestResnorm2) / bestResnorm2;
    else
        marginPct = Inf;
    end
    basinAmbiguous = marginPct < MIN_MARGIN_PCT;
    fprintf('==== Basin selection (L15, TOTAL) ====\n');
    fprintf('  candidate 1 (raw)     resnorm2 = %.6g\n', resnorm2Cand(1));
    fprintf('  candidate 2 (swapped) resnorm2 = %.6g\n', resnorm2Cand(2));
    fprintf('  SELECTED: candidate %d (%s), margin = %.1f%% (threshold %.0f%%)\n', ...
        winnerIdx, candidates(winnerIdx).label, marginPct, MIN_MARGIN_PCT);
    if basinAmbiguous
        warning('run_layer12_total_channel:ambiguousBasin', ...
            ['Basin selection margin only %.1f%% (< %.0f%% threshold) -- AMBIGUOUS, review manually ' ...
             '(expected here more often than the per-channel driver -- no second channel to help resolve it).'], ...
            marginPct, MIN_MARGIN_PCT);
    end
    thetaE = fitRaw; thetaE.candidates = candidates;   % audit trail, matches fit_layer12_stage1_electrical.m's shape
    thetaFull = thetaFullCand{winnerIdx};
    thetaFullSwapped = thetaFullCand{otherIdx};
    if isempty(thetaFull)
        error('run_layer12_total_channel:noWinningFit', 'winning candidate''s Stage 2 fit failed: %s', stage2ErrCand{winnerIdx});
    end
    fprintf('\n');

    % ---- Step 7: diagnostics, cross-checked against bulk_mixed_models' own
    % "Total firing rate (LVN+RVN)" metric (NOT "Nerve firing rate (TOTAL)",
    % which does not exist in R.cells) ---------------------------------------
    figs = plot_layer12_fit_diagnostics(totalRows, 'TOTAL', thetaFull, R, saveDir, 'Total firing rate (LVN+RVN)');

    % ---- Step 8: summary + pattern check ------------------------------------
    paramNames = {'W_c','u_half','f_max_c','f_lo','f_hi','u_M_half','f_max_c2','w_1','w_2'};
    Summary = struct();
    for p = 1:numel(paramNames); Summary.(paramNames{p}) = thetaFull.(paramNames{p}); end
    Summary.check1_M50E100_gt_M100E100 = thetaFull.check1_M50E100_gt_M100E100;
    Summary.check2_M100E10_gt_M50E10   = thetaFull.check2_M100E10_gt_M50E10;
    SummaryTable = struct2table(Summary, 'AsArray', true);
    fprintf('==== Fitted parameter summary (TOTAL) ====\n');
    disp(SummaryTable);

    if ~isempty(saveDir)
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        csvPath = fullfile(saveDir, 'layer12_total_channel_summary.csv');
        writetable(SummaryTable, csvPath);
        matPath = fullfile(saveDir, 'layer12_total_channel_results.mat');
        save(matPath, 'SummaryTable', 'thetaE', 'thetaFull', 'thetaFullSwapped', 'R', 'totalRows');
        fprintf('[saved] %s\n[saved] %s\n', csvPath, matPath);
    end

    fprintf('\n================ TOTAL pattern check ================\n');
    fprintf('  Check1 (M50xE100>M100xE100) = %s   Check2 (M100xE10>M50xE10) = %s\n', ...
        bool2str_local(thetaFull.check1_M50E100_gt_M100E100), bool2str_local(thetaFull.check2_M100E10_gt_M50E10));
    fprintf('=======================================================\n');

    % ---- L17: (u_M_half, w_2) asymptotic ridge flag (post-fit diagnostic --
    % see fit_layer12_stage2_interaction.m / model_fitting_plan.md L17) -----
    fprintf('\n================ L17: (u_M_half, w_2) asymptotic ridge check ================\n');
    if isfield(thetaFull,'asymptoticRidgeFlag') && thetaFull.asymptoticRidgeFlag
        fprintf(2, ['  TOTAL : FLAGGED -- u_M_half=%.1f, w_2=%.4g individually unreliable ' ...
            '(data-coverage limit, L4/L17). kappa2=w_2/u_M_half=%.4g is the more robust quantity.\n'], ...
            thetaFull.u_M_half, thetaFull.w_2, thetaFull.kappa2);
    else
        fprintf('  TOTAL : not flagged this run.\n');
    end
    fprintf('===============================================================================\n');

    results = struct('thetaE', thetaE, 'thetaFull', thetaFull, 'thetaFullSwapped', thetaFullSwapped, ...
        'summary', SummaryTable, 'figs', figs, 'bulkMixedModels', R, 'totalRows', totalRows);
end

% ----------------------------------------------------------------------
function totalRows = build_total_rows_local(rawRows)
% One row per (animal, condition, phase) with label='TOTAL' and
% mean.rate = sum of that trial's RVN + LVN raw mean.rate. Mirrors
% bulk_mixed_models.m's build_total_rate_table grouping convention EXACTLY
% (match by stem+phase, fallback to animal+condition+phase) -- reused, not
% reinvented, per matlab_implementation_instructions.md Section 0.7.
    n = numel(rawRows);
    if isfield(rawRows,'stem') && ~all(cellfun(@isempty,{rawRows.stem}))
        rk = arrayfun(@(r) sprintf('%s|%s', char(rawRows(r).stem), rawRows(r).phase), 1:n, 'uni', 0);
    else
        rk = arrayfun(@(r) sprintf('%s|%s|%s', rawRows(r).animal, rawRows(r).condition, rawRows(r).phase), 1:n, 'uni', 0);
    end
    uk = unique(rk, 'stable');
    totalRows = repmat(empty_total_row_local(), 0, 1);
    for k = 1:numel(uk)
        idx = find(strcmp(rk, uk{k}));
        tot = 0; cnt = 0;
        for ii = idx(:)'
            if isfield(rawRows(ii),'mean') && isfield(rawRows(ii).mean,'rate') && isfinite(rawRows(ii).mean.rate)
                tot = tot + rawRows(ii).mean.rate; cnt = cnt + 1;
            end
        end
        if cnt == 0; continue; end
        r = empty_total_row_local();
        r.animal = rawRows(idx(1)).animal;
        r.condition = rawRows(idx(1)).condition;
        r.phase = rawRows(idx(1)).phase;
        r.label = 'TOTAL';
        if isfield(rawRows(idx(1)),'stem'); r.stem = rawRows(idx(1)).stem; end
        r.mean.rate = tot;
        totalRows(end+1) = r; %#ok<AGROW>
    end
end

function r = empty_total_row_local()
    r = struct('animal','','condition','','phase','','label','', 'mean',struct('rate',NaN), 'stem','');
end

% ----------------------------------------------------------------------
function thetaE = fit_stage1_single_local(totalRows)
% Single-group 3-param [W_c_TOTAL, u_half, f_max_c] fit against TOTAL
% E-alone trials. NOTE: with only one group, there is nothing to pool
% against -- this is the ORIGINAL (pre-L12) single-channel procedure, and
% the exact (u_half,f_max_c) swap symmetry (model_fitting_plan.md L12) is
% NOT resolved here. Resolved (if it resolves at all) via the dual-basin
% Stage 2 check in the caller, same as L13 validated for the per-channel case.
    isRec = strcmpi({totalRows.phase}, 'recovery');
    n = numel(totalRows);
    uM = nan(n,1); uE = nan(n,1);
    for i = 1:n; [uM(i), uE(i)] = parse_me_local(totalRows(i).condition); end
    isEalone = (uM == 0) & (uE > 0);
    sel = find(isRec(:) & isEalone(:));
    rate = arrayfun(@(i) totalRows(i).mean.rate, sel);
    uEsel = uE(sel);
    finite = isfinite(rate) & isfinite(uEsel);
    sel = sel(finite); rate = rate(finite); uEsel = uEsel(finite); %#ok<NASGU>
    fprintf('[stage1-TOTAL] %d E-alone trials (u_E levels: %s)\n', numel(sel), mat2str(unique(uEsel)'));
    if numel(sel) < 3
        error('run_layer12_total_channel:tooFew', 'only %d TOTAL E-alone trials -- need >= 3 to fit [W_c,u_half,f_max_c].', numel(sel));
    end

    resFun = @(th) rate(:) - th(1) .* model_layer12_equations('Phi_c', uEsel(:), th(3)) ...
                                    .* model_layer12_equations('h',     uEsel(:), th(2));
    lb = [0 0 0]; ub = [Inf Inf Inf];
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12, 'MaxFunctionEvaluations',5000);
    nRestarts = 20;
    scaleGuess = max(rate(isfinite(rate) & rate>0), [], 'omitnan');
    if isempty(scaleGuess) || ~isfinite(scaleGuess) || scaleGuess <= 0; scaleGuess = 1; end

    results = repmat(struct('theta',[nan nan nan],'resnorm',NaN,'exitflag',NaN), nRestarts, 1);
    for r = 1:nRestarts
        x0 = [ scaleGuess * 10.^(rand()*4-2), 10.^(rand()*4), 10.^(rand()*4) ];
        try
            [th, resnorm, ~, exitflag] = lsqnonlin(resFun, x0, lb, ub, opts);
        catch ME
            warning('run_layer12_total_channel:restart', 'restart %d failed (%s)', r, ME.message);
            th = [NaN NaN NaN]; resnorm = NaN; exitflag = -99;
        end
        results(r) = struct('theta',th,'resnorm',resnorm,'exitflag',exitflag);
    end
    resnorms = arrayfun(@(s) s.resnorm, results);
    [bestResnorm, bestIdx] = min(resnorms);
    bt = results(bestIdx).theta;
    fprintf('[stage1-TOTAL] %d/%d restarts converged. best resnorm=%.6g, theta=[W_c_TOTAL=%.4g, u_half=%.4g, f_max_c=%.4g]\n', ...
        nnz(arrayfun(@(s) s.exitflag>0, results)), nRestarts, bestResnorm, bt(1), bt(2), bt(3));
    fprintf(['[stage1-TOTAL] raw u_half/f_max_c/W_c_TOTAL NOT individually meaningful here (exact swap symmetry, ' ...
        'model_fitting_plan.md L12 -- no second channel to pool against). Identifiable invariant ' ...
        'sqrt(u_half*f_max_c) = %.4g.\n'], sqrt(bt(2)*bt(3)));

    thetaE = struct('W_c_TOTAL', bt(1), 'u_half', bt(2), 'f_max_c', bt(3), ...
        'rho', sqrt(bt(2)*bt(3)), 'resnorm', bestResnorm);
end

function [M,E] = parse_me_local(cond)
% Replicates parse_me in bulk_mixed_models.m exactly: M(\d+) / E(\d+),
% default 0 if the respective token is absent.
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end

function s = bool2str_local(tf)
    if isnan(tf); s = 'N/A'; elseif tf; s = 'PASS'; else; s = 'FAIL'; end
end
