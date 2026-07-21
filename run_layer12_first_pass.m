function results = run_layer12_first_pass(rawRows, saveDir)
%RUN_LAYER12_FIRST_PASS Top-level driver for the Layer 1-2 first-pass model fit.
%   results = run_layer12_first_pass(rawRows)            % rawRows from existing pipeline
%   results = run_layer12_first_pass(rawRows, saveDir)    % also saves figures/CSV
%
%   MAJOR REVISION per L13/L15 (model_fitting_plan.md): Stage 1's (u_half,
%   f_max_c) is PROVEN non-unique (L12, exact swap symmetry) -- Stage 1
%   returns BOTH equally-valid candidates (fit_layer12_stage1_electrical.m's
%   thetaE.candidates(1:2)), and Stage 2 reliably picks a winner between them
%   (L13). This driver evaluates Stage 2 for BOTH candidates x BOTH channels
%   (4 fits total), SELECTS the winning candidate by total resnorm summed
%   across channels (a joint decision, since u_half/f_max_c are shared), and
%   propagates ONLY the winner's already-computed fits to diagnostics/summary
%   -- it does not silently use whatever Stage 1's optimizer converged to
%   first, and it does not re-fit a third time.
%
%   rawRows must have the shape produced by run_pipeline_bulk.m's harvest
%   (one row per file x phase x channel, with .animal/.condition/.label/
%   .phase/.mean.rate/... -- see matlab_implementation_instructions.md
%   Section 0.5). If you already have out.raw from run_pipeline_bulk, pass
%   that directly as rawRows.
%
%   NOTE: this pass fits RATE ONLY (model_fitting_plan.md Section 5.4 --
%   CV2/FWHM extension -- is explicitly deferred until this fit is reviewed).

    if nargin < 2; saveDir = ''; end
    channels = {'RVN','LVN'};
    MIN_MARGIN_PCT = 10;   % below this, flag basin selection AMBIGUOUS rather than trusting it (L15)

    % ---- Step 1: unit tests gate -- do not fit with failing tests --------
    fprintf('==== Step 1: unit tests (test_layer12_equations) ====\n');
    testResults = runtests('test_layer12_equations');
    nFailed = nnz(~[testResults.Passed]);
    if nFailed > 0
        disp(testResults);   % printed for visibility before erroring
        error('run_layer12_first_pass:testsFailed', ...
            '%d/%d unit tests FAILED -- aborting before any fitting (do not fit against real data with failing unit tests).', ...
            nFailed, numel(testResults));
    end
    fprintf('All %d unit tests passed.\n\n', numel(testResults));

    % ---- Step 8 (bulk_mixed_models cross-check): computed here, EARLY, so
    % it's ready before Step 7's per-channel diagnostics need it (Section 6).
    fprintf('==== Step 8 (computed early -- needed by Step 7 diagnostics): bulk_mixed_models cross-check ====\n');
    R = [];
    if exist('bulk_mixed_models','file')==2 && exist('bulk_compile','file')==2
        try
            normRows = bulk_compile(rawRows);
            out = struct('raw', rawRows, 'norm', normRows);
            % includeSystemic=false: our cross-check only needs the nerve/rate
            % metric (plot_layer12_fit_diagnostics.m filters R.cells to
            % "Nerve firing rate (<channel>)"); skipping systemic (HR/HRV)
            % avoids touching whatever gemsplots_queue.mat/buildFileQueue state
            % happens to exist in the current directory (that loader can pop
            % its own interactive UI, unrelated to this fit).
            R = bulk_mixed_models(out, saveDir, [], false);
        catch ME
            warning('run_layer12_first_pass:bmm', ...
                'bulk_mixed_models cross-check failed (%s) -- proceeding without it.', ME.message);
        end
    else
        warning('run_layer12_first_pass:bmmMissing', ...
            'bulk_mixed_models/bulk_compile not on path -- proceeding without the cross-check.');
    end
    fprintf('\n');

    % ---- Step 2: Stage 1, POOLED across both channels (L12) --------------
    % Called ONCE. Returns thetaE.candidates(1:2) -- BOTH Stage-1-equivalent
    % solutions, not a single flat u_half/f_max_c (L12/L15).
    fprintf('==== Step 2: Stage 1 (pooled E-alone fit, both channels, both candidates) ====\n');
    thetaE = [];
    stage1Error = '';
    try
        thetaE = fit_layer12_stage1_electrical(rawRows);
    catch ME
        warning('run_layer12_first_pass:stage1','Stage 1 (pooled) failed (%s)', ME.message);
        stage1Error = ME.message;
    end
    fprintf('\n');

    % ---- Step 3: Stage 2 for BOTH candidates x BOTH channels (4 fits) -----
    fprintf('==== Step 3: Stage 2, both candidates x both channels ====\n');
    resnorm2 = nan(2, numel(channels));       % rows = candidate index, cols = channel index
    thetaFullAll = cell(2, numel(channels));  % same indexing
    stage2Error = cell(2, numel(channels));
    if ~isempty(thetaE)
        for ci = 1:2
            for k = 1:numel(channels)
                try
                    [thetaFullAll{ci,k}, resnorm2(ci,k)] = fit_layer12_stage2_interaction( ...
                        rawRows, channels{k}, thetaE.candidates(ci));
                catch ME
                    warning('run_layer12_first_pass:stage2', '%s candidate %d (%s): fit failed (%s)', ...
                        channels{k}, ci, thetaE.candidates(ci).label, ME.message);
                    stage2Error{ci,k} = ME.message;
                end
            end
        end
    end
    fprintf('\n');

    % ---- Step 4: select the winning candidate (L15), by TOTAL resnorm
    % summed across channels (joint/shared-parameter decision) --------------
    winnerIdx = []; otherIdx = []; marginPct = NaN; basinAmbiguous = NaN;
    if ~isempty(thetaE)
        resnormForSelection = resnorm2;
        resnormForSelection(isnan(resnormForSelection)) = Inf;   % a failed fit must not win by NaN-propagation
        totalResnorm = sum(resnormForSelection, 2);
        [bestTotal, winnerIdx] = min(totalResnorm);
        otherIdx = 3 - winnerIdx;
        if isfinite(totalResnorm(otherIdx)) && isfinite(bestTotal) && bestTotal > 0
            marginPct = 100 * (totalResnorm(otherIdx) - bestTotal) / bestTotal;
        else
            marginPct = Inf;   % other candidate failed entirely -- no real ambiguity
        end
        basinAmbiguous = marginPct < MIN_MARGIN_PCT;
        fprintf('==== Step 4: basin selection (L15) ====\n');
        fprintf('  candidate 1 (%s) total resnorm = %.6g\n', thetaE.candidates(1).label, totalResnorm(1));
        fprintf('  candidate 2 (%s) total resnorm = %.6g\n', thetaE.candidates(2).label, totalResnorm(2));
        fprintf('  SELECTED: candidate %d (%s), margin = %.1f%% (threshold %.0f%%)\n', ...
            winnerIdx, thetaE.candidates(winnerIdx).label, marginPct, MIN_MARGIN_PCT);
        if basinAmbiguous
            warning('run_layer12_first_pass:ambiguousBasin', ...
                ['Stage 1 basin selection margin only %.1f%% (< %.0f%% threshold) -- ' ...
                 'AMBIGUOUS, review manually before trusting auto-selected results.'], ...
                marginPct, MIN_MARGIN_PCT);
        end
        fprintf('\n');
    end

    % ---- Step 5: per-channel agreement check, every run (L15) --------------
    channelDisagreement = NaN;
    if ~isempty(thetaE)
        resnormForAgreement = resnorm2;
        resnormForAgreement(isnan(resnormForAgreement)) = Inf;
        [~, winnerPerChannel] = min(resnormForAgreement, [], 1);   % 1 x nChannels
        channelDisagreement = ~all(winnerPerChannel == winnerPerChannel(1));
        fprintf('==== Step 5: per-channel basin agreement check ====\n');
        for k = 1:numel(channels)
            fprintf('  %-4s prefers candidate %d (%s)\n', channels{k}, winnerPerChannel(k), ...
                thetaE.candidates(winnerPerChannel(k)).label);
        end
        if channelDisagreement
            warning('run_layer12_first_pass:channelDisagreement', ...
                ['Channels DISAGREE on which Stage 1 candidate is correct -- this is evidence AGAINST the ' ...
                 'L12 pooling assumption (shared u_half/f_max_c across channels), not just a selection ' ...
                 'nuance. Reporting prominently; NOT silently defaulting to the total''s answer.']);
        else
            fprintf('  Channels AGREE -- both independently prefer the same candidate (supports L12 pooling).\n');
        end
        fprintf('\n');
    end

    % ---- Step 6: promote the winning candidate's ALREADY-COMPUTED fits ---
    % (no third fit) -----------------------------------------------------
    perChannel = struct();
    for k = 1:numel(channels)
        ch = channels{k};
        entry = struct('thetaFull',[], 'resnorm2',NaN, 'figs',[], 'error','');
        if isempty(thetaE)
            entry.error = ['Stage 1 (pooled) failed, cannot run Stage 2: ' stage1Error];
        elseif isempty(thetaFullAll{winnerIdx,k})
            entry.error = ['winning-candidate Stage 2 fit failed: ' stage2Error{winnerIdx,k}];
        else
            entry.thetaFull = thetaFullAll{winnerIdx,k};
            entry.resnorm2 = resnorm2(winnerIdx,k);
        end
        perChannel.(ch) = entry;
    end

    % ---- Step 7: diagnostics, using ONLY the winning candidate's fits ----
    fprintf('==== Step 7: diagnostics (winning candidate only) ====\n');
    for k = 1:numel(channels)
        ch = channels{k};
        if ~isempty(perChannel.(ch).thetaFull)
            try
                perChannel.(ch).figs = plot_layer12_fit_diagnostics(rawRows, ch, perChannel.(ch).thetaFull, R, saveDir);
            catch ME
                warning('run_layer12_first_pass:diagnostics', '%s: diagnostics failed (%s)', ch, ME.message);
            end
        else
            fprintf('  %s: skipped (no winning-candidate fit available -- %s)\n', ch, perChannel.(ch).error);
        end
    end
    fprintf('\n');

    % ---- Step 9: summary tables + audit trail (both candidates saved) ----
    % u_half/f_max_c reported ONCE (SHARED, from the SELECTED candidate) --
    % NOT duplicated per channel, per L12/L15. --------------------------
    thetaE_selected = []; thetaE_rejected = [];
    if isempty(thetaE)
        sel_u_half = NaN; sel_f_max_c = NaN; selectedCandidateLabel = '';
    else
        thetaE_selected = thetaE.candidates(winnerIdx);
        thetaE_rejected = thetaE.candidates(otherIdx);
        sel_u_half = thetaE_selected.u_half;
        sel_f_max_c = thetaE_selected.f_max_c;
        selectedCandidateLabel = thetaE_selected.label;
    end
    % AsArray=true: without it, struct2table on a SCALAR struct tries to infer
    % a row count per field and chokes on char fields (0 rows when '', vs 1
    % row for the numeric scalars) -- forces "one row" instead.
    Shared = struct2table(struct( ...
        'u_half', sel_u_half, 'f_max_c', sel_f_max_c, ...
        'selectedCandidate', selectedCandidateLabel, ...
        'selectionMarginPct', marginPct, 'basinAmbiguous', basinAmbiguous, 'channelDisagreement', channelDisagreement, ...
        'resnormPooled', pick(thetaE,'resnormPooled'), ...
        'resnormIndependentTotal', pick(thetaE,'resnormIndependentTotal'), ...
        'pooledVsIndependentPctChange', pick(thetaE,'pooledVsIndependentPctChange'), ...
        'stage1Error', stage1Error), 'AsArray', true);

    perChanParamNames = {'W_c','f_lo','f_hi','u_M_half','f_max_c2','w_1','w_2'};
    rows = {};
    for k = 1:numel(channels)
        ch = channels{k}; entry = perChannel.(ch);
        r = struct('channel',ch);
        if isempty(entry.thetaFull)
            for p = 1:numel(perChanParamNames); r.(perChanParamNames{p}) = NaN; end
            r.check1_M50E100_gt_M100E100 = NaN; r.check2_M100E10_gt_M50E10 = NaN;
            r.fitError = entry.error;
        else
            tf = entry.thetaFull;
            for p = 1:numel(perChanParamNames); r.(perChanParamNames{p}) = tf.(perChanParamNames{p}); end
            r.check1_M50E100_gt_M100E100 = tf.check1_M50E100_gt_M100E100;
            r.check2_M100E10_gt_M50E10   = tf.check2_M100E10_gt_M50E10;
            r.fitError = '';
        end
        rows{end+1} = r; %#ok<AGROW>
    end
    PerChannel = struct2table(vertcat(rows{:}));

    % audit trail: BOTH candidates' resnorm2 across BOTH channels -- do not
    % discard the rejected candidate once selection is made (L15/Section 9).
    AuditRows = {};
    if ~isempty(thetaE)
        for ci = 1:2
            for k = 1:numel(channels)
                AuditRows{end+1} = struct('candidateIdx',ci, 'candidateLabel',thetaE.candidates(ci).label, ...
                    'channel',channels{k}, 'resnorm2',resnorm2(ci,k), 'selected',(ci==winnerIdx)); %#ok<AGROW>
            end
        end
    end
    if ~isempty(AuditRows); Audit = struct2table(vertcat(AuditRows{:})); else; Audit = table(); end

    fprintf('==== Step 9: fitted parameter summary ====\n');
    fprintf('-- shared (u_half, f_max_c; pooled across RVN+LVN, SELECTED candidate) --\n');
    disp(Shared);
    fprintf('-- per-channel (W_c + Stage 2 interaction params, SELECTED candidate) --\n');
    disp(PerChannel);
    fprintf('-- audit trail (both candidates x both channels, resnorm2) --\n');
    disp(Audit);

    if ~isempty(saveDir)
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        sharedCsvPath = fullfile(saveDir, 'layer12_first_pass_summary_shared.csv');
        perChanCsvPath = fullfile(saveDir, 'layer12_first_pass_summary_perchannel.csv');
        auditCsvPath = fullfile(saveDir, 'layer12_first_pass_basin_audit.csv');
        writetable(Shared, sharedCsvPath);
        writetable(PerChannel, perChanCsvPath);
        if ~isempty(Audit); writetable(Audit, auditCsvPath); end
        matPath = fullfile(saveDir, 'layer12_first_pass_results.mat');
        save(matPath, 'Shared', 'PerChannel', 'Audit', 'thetaE', 'thetaE_selected', 'thetaE_rejected', ...
            'perChannel', 'winnerIdx', 'R');
        fprintf('[saved] %s\n[saved] %s\n[saved] %s\n[saved] %s\n', sharedCsvPath, perChanCsvPath, auditCsvPath, matPath);
    end

    % ---- Step 10: prominent pass/fail report (from SELECTED candidate) ---
    fprintf('\n================ Section 5 pattern check (BOTH channels, SELECTED candidate) ================\n');
    for k = 1:numel(channels)
        ch = channels{k}; entry = perChannel.(ch);
        if isempty(entry.thetaFull)
            fprintf('  %-4s : FIT FAILED (%s)\n', ch, entry.error);
            continue;
        end
        tf = entry.thetaFull;
        fprintf('  %-4s : Check1 (M50xE100>M100xE100) = %s   Check2 (M100xE10>M50xE10) = %s\n', ...
            ch, bool2str(tf.check1_M50E100_gt_M100E100), bool2str(tf.check2_M100E10_gt_M50E10));
    end
    fprintf('===============================================================================================\n');

    % ---- Step 11: pooled-vs-independent Stage 1 residual comparison (L12) --
    fprintf('\n================ L12 validation: pooled-vs-independent Stage 1 ================\n');
    if isempty(thetaE)
        fprintf('  Stage 1 (pooled) FAILED -- no validation available (%s)\n', stage1Error);
    else
        fprintf('  pooled total resnorm      = %.6g\n', thetaE.resnormPooled);
        fprintf('  independent total resnorm = %.6g\n', thetaE.resnormIndependentTotal);
        fprintf('  pooled vs independent     = %+.1f%%\n', thetaE.pooledVsIndependentPctChange);
        if thetaE.pooledVsIndependentPctChange > 50
            fprintf(2, '  *** POOLING INCREASED RESIDUAL SHARPLY *** -- evidence AGAINST the shared-physiology assumption.\n');
        else
            fprintf('  pooling did not increase residual sharply -- supports the shared-physiology assumption.\n');
        end
    end
    fprintf('=================================================================================\n');

    % ---- Step 12: basin selection margin + ambiguity/disagreement warnings,
    % prominently (new, required reporting, not optional diagnostic output) --
    fprintf('\n================ L13/L15: basin selection outcome ================\n');
    if isempty(thetaE)
        fprintf('  Stage 1 (pooled) FAILED -- no basin selection available.\n');
    else
        fprintf('  Selected candidate %d (%s); rejected candidate %d (%s)\n', ...
            winnerIdx, thetaE.candidates(winnerIdx).label, otherIdx, thetaE.candidates(otherIdx).label);
        fprintf('  Selection margin: %.1f%% (threshold %.0f%%) -- %s\n', marginPct, MIN_MARGIN_PCT, ...
            ternary_str_local(basinAmbiguous, 'AMBIGUOUS, review manually', 'decisive'));
        fprintf('  Per-channel agreement: %s\n', ...
            ternary_str_local(channelDisagreement, 'CHANNELS DISAGREE -- pooling assumption in question', 'channels agree'));
        fprintf('  Identifiable invariants (same for both candidates, L12): sqrt(u_half*f_max_c)=%.3f, W_c_LVN/W_c_RVN=%.3f\n', ...
            thetaE.rho, thetaE.W_c_ratio);
    end
    fprintf('====================================================================\n');

    % ---- Step 12b: L17 (u_M_half, w_2) asymptotic ridge flag, BOTH channels --
    % Post-fit diagnostic only (see fit_layer12_stage2_interaction.m) -- report
    % prominently here too if flagged, per model_fitting_plan.md L17 / Section 9.
    fprintf('\n================ L17: (u_M_half, w_2) asymptotic ridge check ================\n');
    anyRidgeFlagged = false;
    for k = 1:numel(channels)
        ch = channels{k}; entry = perChannel.(ch);
        if isempty(entry.thetaFull); continue; end
        tf = entry.thetaFull;
        if isfield(tf,'asymptoticRidgeFlag') && tf.asymptoticRidgeFlag
            anyRidgeFlagged = true;
            fprintf(2, ['  %-4s : FLAGGED -- u_M_half=%.1f, w_2=%.4g individually unreliable ' ...
                '(data-coverage limit, L4/L17). kappa2=w_2/u_M_half=%.4g is the more robust quantity.\n'], ...
                ch, tf.u_M_half, tf.w_2, tf.kappa2);
        else
            fprintf('  %-4s : not flagged (u_M_half=%.1f within tested-range multiple)\n', ch, tf.u_M_half);
        end
    end
    if ~anyRidgeFlagged
        fprintf('  No channel flagged this run.\n');
    end
    fprintf('===============================================================================\n');

    results = struct('shared', Shared, 'perChannel', PerChannel, 'audit', Audit, 'thetaE', thetaE, ...
        'thetaE_selected', thetaE_selected, 'thetaE_rejected', thetaE_rejected, ...
        'selectionMarginPct', marginPct, 'basinAmbiguous', basinAmbiguous, 'channelDisagreement', channelDisagreement, ...
        'perChannelDetail', perChannel, 'bulkMixedModels', R);
end

% ----------------------------------------------------------------------
function v = pick(s, field)
    if isempty(s) || ~isfield(s, field); v = NaN; else; v = s.(field); end
end

function s = bool2str(tf)
    if isnan(tf); s = 'N/A'; elseif tf; s = 'PASS'; else; s = 'FAIL'; end
end

function s = ternary_str_local(cond, ifTrue, ifFalse)
    if isnan(cond); s = 'N/A'; elseif cond; s = ifTrue; else; s = ifFalse; end
end
