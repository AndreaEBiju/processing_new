function thetaE = fit_layer12_stage1_electrical(rawRows)
%FIT_LAYER12_STAGE1_ELECTRICAL Fit E-alone sub-model POOLED across RVN and LVN.
%
%   thetaE = fit_layer12_stage1_electrical(rawRows)
%
%   NOTE signature change (L12): no longer takes a channelLabel argument --
%   this function now ALWAYS fits both channels together. u_half and f_max_c
%   are SHARED across channels; W_c_RVN and W_c_LVN are fit separately.
%   Empirically forced: independent per-channel fits of (u_half, f_max_c)
%   landed in different parameter-space basins with near-identical residuals
%   -- direct evidence of the ridge anticipated in model_fitting_plan.md
%   Section 5.5 item 3. u_half/f_max_c describe fibre-level entrainment/
%   following physiology with no expected left/right difference; W_c absorbs
%   per-electrode coupling, which plausibly DOES differ by channel (L12).
%
%   Uses RAW (non-normalized) per-trial mean rate from rawRows(i).mean.rate,
%   phase == 'recovery', filtered to E-alone conditions, BOTH channels pooled
%   into ONE fit.
%
%   Objective: nonlinear least squares via lsqnonlin, jointly over both
%   channels:
%     minimize sum_over_BOTH_channels_and_trials[
%       (observed_rate - W_c_{channel}*Phi_c(u_E,f_max_c)*h(u_E,u_half)).^2 ]
%   -- 4 free parameters (u_half, f_max_c, W_c_RVN, W_c_LVN) against
%   3 E-alone conditions x 2 channels = 6 data-generating groups.
%
%   REQUIRED VALIDATION (L12, not optional): also runs the original
%   independent-per-channel version (2 separate 3-parameter fits) and reports
%   its total residual alongside the pooled fit's, so the pooling assumption
%   itself is checked rather than assumed.
%
%   RETURNS thetaE as a struct with a `candidates` field: a 1x2 struct array,
%   each element holding one of the two Stage-1-EQUALLY-VALID solutions --
%   candidates(1) ("raw", whatever lsqnonlin converged to) and candidates(2)
%   ("swapped", the analytically-computed alternative -- u_half<->f_max_c,
%   W_c_<channel> rescaled by f_max_c/u_half). This is NOT a flat
%   single-solution struct (L12/L15): raw u_half/f_max_c/W_c are not
%   individually identifiable, so returning only one candidate would silently
%   commit to an arbitrary, optimizer-initialization-dependent choice between
%   two exactly-equally-good fits. The caller (run_layer12_first_pass.m) MUST
%   evaluate Stage 2 for BOTH candidates and select the winner (L15) -- do not
%   pick one here. thetaE.rho and thetaE.W_c_ratio (the two identifiable
%   invariants, L12) are the SAME for both candidates by construction and are
%   stored once at the top level, not duplicated per candidate.

    channels = {'RVN','LVN'};

    % ---- gather E-alone trials for BOTH channels -------------------------
    isRec = strcmpi({rawRows.phase}, 'recovery');
    n = numel(rawRows);
    uM = nan(n,1); uE = nan(n,1);
    for i = 1:n; [uM(i), uE(i)] = parse_me_local(rawRows(i).condition); end
    isEalone = (uM == 0) & (uE > 0);

    perChanData = struct();
    for c = 1:numel(channels)
        ch = channels{c};
        isCh = strcmpi({rawRows.label}, ch);
        sel = find(isRec(:) & isCh(:) & isEalone(:));
        rate = arrayfun(@(i) rawRows(i).mean.rate, sel);
        uEsel = uE(sel);
        finite = isfinite(rate) & isfinite(uEsel);
        sel = sel(finite); rate = rate(finite); uEsel = uEsel(finite);
        fprintf('[stage1] %s: %d E-alone trials (u_E levels: %s)\n', ch, numel(sel), mat2str(unique(uEsel)'));
        perChanData.(ch) = struct('rate', rate, 'uE', uEsel);
    end

    totalN = numel(perChanData.RVN.rate) + numel(perChanData.LVN.rate);
    nParamsPooled = 4;
    if totalN < nParamsPooled
        error('fit_layer12_stage1_electrical:tooFew', ...
            'only %d E-alone trials (both channels combined) -- need >= %d to fit [u_half,f_max_c,W_c_RVN,W_c_LVN].', ...
            totalN, nParamsPooled);
    end
    for c = 1:numel(channels)
        ch = channels{c};
        if numel(unique(perChanData.(ch).uE)) < 3
            warning('fit_layer12_stage1_electrical:fewLevels', ...
                ['%s: only %d distinct u_E levels present among E-alone trials -- ' ...
                 'expect (u_half,f_max_c) to be weakly identified (see model_fitting_plan.md Sec 5.5 item 3).'], ...
                ch, numel(unique(perChanData.(ch).uE)));
        end
    end

    lbPooled = [0 0 0 0]; ubPooled = [Inf Inf Inf Inf];
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12, 'MaxFunctionEvaluations',5000);
    nRestarts = 20;

    % ---- POOLED fit: theta = [u_half, f_max_c, W_c_RVN, W_c_LVN] ----------
    resFunPooled = @(th) pooled_residuals(th, perChanData);
    scaleGuess = max([perChanData.RVN.rate(:); perChanData.LVN.rate(:)], [], 'omitnan');
    if isempty(scaleGuess) || ~isfinite(scaleGuess) || scaleGuess <= 0; scaleGuess = 1; end

    poolResults = repmat(struct('theta',nan(1,4),'resnorm',NaN,'exitflag',NaN), nRestarts, 1);
    for r = 1:nRestarts
        x0 = [ 10.^(rand()*4), ...                    % u_half   ~ [1, 1e4] Hz
               10.^(rand()*4), ...                    % f_max_c  ~ [1, 1e4] Hz
               scaleGuess * 10.^(rand()*4-2), ...      % W_c_RVN
               scaleGuess * 10.^(rand()*4-2) ];         % W_c_LVN
        try
            [th, resnorm, ~, exitflag] = lsqnonlin(resFunPooled, x0, lbPooled, ubPooled, opts);
        catch ME
            warning('fit_layer12_stage1_electrical:restart','pooled restart %d failed (%s)', r, ME.message);
            th = nan(1,4); resnorm = NaN; exitflag = -99;
        end
        poolResults(r) = struct('theta',th,'resnorm',resnorm,'exitflag',exitflag);
    end
    resnorms = arrayfun(@(s) s.resnorm, poolResults);
    [bestResnorm, bestIdx] = min(resnorms);
    bt = poolResults(bestIdx).theta;

    fprintf('[stage1-pooled] %d/%d restarts converged (exitflag>0).\n', ...
        nnz(arrayfun(@(s) s.exitflag>0, poolResults)), nRestarts);
    fprintf('[stage1-pooled] best resnorm=%.6g, exitflag=%d\n', bestResnorm, poolResults(bestIdx).exitflag);
    fprintf('[stage1-pooled] theta=[u_half=%.4g, f_max_c=%.4g, W_c_RVN=%.4g, W_c_LVN=%.4g]\n', ...
        bt(1), bt(2), bt(3), bt(4));

    near = resnorms <= bestResnorm * 1.01;
    thetaNear = vertcat(poolResults(near).theta);
    names = {'u_half','f_max_c','W_c_RVN','W_c_LVN'};
    fprintf('[stage1-pooled] spread among %d near-best restarts:\n', nnz(near));
    for k = 1:4
        fprintf('    %-10s [%.4g, %.4g]\n', names{k}, min(thetaNear(:,k)), max(thetaNear(:,k)));
    end
    fitRaw = struct('u_half', bt(1), 'f_max_c', bt(2), 'W_c_RVN', bt(3), 'W_c_LVN', bt(4), ...
        'resnormPooled', bestResnorm, 'nTrials', totalN, 'nRestartsConverged', nnz(near));

    % ---- L12: raw u_half/f_max_c/W_c are NOT individually identifiable ----
    % (PROVEN exact one-parameter symmetry, model_fitting_plan.md L12 -- swapping
    % u_half<->f_max_c and rescaling W_c -> W_c*(f_max_c/u_half) leaves the fit
    % EXACTLY unchanged for every u_E). Only these two combinations are
    % identifiable even with infinite noiseless data -- print both the raw
    % values (for debugging/reproducing this specific fit) AND the invariants,
    % explicitly labeled, so raw values are never mistaken for meaningful ones.
    % Both are the SAME for both candidates below -- that is the point of L12.
    rho   = sqrt(fitRaw.u_half * fitRaw.f_max_c);
    ratio = fitRaw.W_c_LVN / fitRaw.W_c_RVN;
    fprintf(['Stage 1 raw values (NOT individually meaningful -- exact swap ' ...
             'symmetry, see model_fitting_plan.md L12):\n']);
    fprintf('  u_half = %.3f, f_max_c = %.3f, W_c_RVN = %.3f, W_c_LVN = %.3f\n', ...
        fitRaw.u_half, fitRaw.f_max_c, fitRaw.W_c_RVN, fitRaw.W_c_LVN);
    fprintf(['Stage 1 IDENTIFIABLE invariants (these ARE meaningful, and are the\n' ...
             'SAME for both candidates below -- that is the point of L12):\n' ...
             '  sqrt(u_half*f_max_c) = %.3f\n  W_c_LVN/W_c_RVN = %.3f\n'], rho, ratio);

    % ---- L13/L15: package BOTH Stage-1-equivalent candidates -----------
    % Closed-form, no re-fit: swapping u_half<->f_max_c and rescaling W_c by
    % the compensating factor f_max_c/u_half is EXACTLY as good a Stage 1 fit
    % (L12). Stage 2 (via g1(u_M), which has no compensating freedom under
    % this swap) reliably picks a winner between the two (L13) -- but that
    % selection happens in the driver, once both channels' Stage 2 fits for
    % both candidates exist, not here.
    rescaleFactor = fitRaw.f_max_c / fitRaw.u_half;
    thetaE = struct();
    thetaE.candidates(1) = struct('u_half', fitRaw.u_half, 'f_max_c', fitRaw.f_max_c, ...
        'W_c_RVN', fitRaw.W_c_RVN, 'W_c_LVN', fitRaw.W_c_LVN, 'label', 'raw');
    thetaE.candidates(2) = struct('u_half', fitRaw.f_max_c, 'f_max_c', fitRaw.u_half, ...
        'W_c_RVN', fitRaw.W_c_RVN*rescaleFactor, 'W_c_LVN', fitRaw.W_c_LVN*rescaleFactor, 'label', 'swapped');
    thetaE.rho = rho;
    thetaE.W_c_ratio = ratio;
    thetaE.resnormPooled = fitRaw.resnormPooled;
    thetaE.nTrials = fitRaw.nTrials;
    thetaE.nRestartsConverged = fitRaw.nRestartsConverged;

    % ---- REQUIRED VALIDATION (L12): independent per-channel fit ----------
    % Reruns the ORIGINAL (now-superseded) per-channel 3-param procedure and
    % compares total residual against the pooled fit -- this is the actual
    % check for whether the shared-physiology pooling assumption is
    % defensible for this dataset, not something to assume.
    indepTotalResnorm = 0;
    indepTheta = struct();
    for c = 1:numel(channels)
        ch = channels{c};
        [thc, resnormc] = fit_one_channel_only(perChanData.(ch), scaleGuess, nRestarts, opts);
        indepTheta.(ch) = struct('W_c',thc(1), 'u_half',thc(2), 'f_max_c',thc(3), 'resnorm',resnormc);
        indepTotalResnorm = indepTotalResnorm + resnormc;
        fprintf('[stage1-independent] %s: theta=[W_c=%.4g, u_half=%.4g, f_max_c=%.4g], resnorm=%.6g\n', ...
            ch, thc(1), thc(2), thc(3), resnormc);
    end

    pctChange = 100 * (bestResnorm - indepTotalResnorm) / max(indepTotalResnorm, eps);
    fprintf(['\n[L12 VALIDATION] Pooled-vs-independent Stage 1 residual comparison:\n' ...
        '    pooled total resnorm     = %.6g  (u_half=%.4g, f_max_c=%.4g shared)\n' ...
        '    independent total resnorm = %.6g  (RVN: u_half=%.4g,f_max_c=%.4g | LVN: u_half=%.4g,f_max_c=%.4g)\n' ...
        '    pooled vs independent     = %+.1f%%\n'], ...
        bestResnorm, bt(1), bt(2), indepTotalResnorm, ...
        indepTheta.RVN.u_half, indepTheta.RVN.f_max_c, indepTheta.LVN.u_half, indepTheta.LVN.f_max_c, pctChange);
    if pctChange > 50
        fprintf(2, ['[L12 VALIDATION] *** POOLING INCREASED RESIDUAL SHARPLY (%+.1f%%) *** -- this is ' ...
            'evidence AGAINST the shared-physiology assumption for this dataset. Reporting prominently, ' ...
            'not silently preferring the pooled result.\n'], pctChange);
    else
        fprintf('[L12 VALIDATION] pooling did not increase residual sharply -- supports the shared-physiology assumption.\n');
    end
    thetaE.independentValidation = indepTheta;
    thetaE.resnormIndependentTotal = indepTotalResnorm;
    thetaE.pooledVsIndependentPctChange = pctChange;
end

% ----------------------------------------------------------------------
function res = pooled_residuals(th, perChanData)
    u_half = th(1); f_max_c = th(2); W_c_RVN = th(3); W_c_LVN = th(4);
    resRVN = perChanData.RVN.rate(:) - W_c_RVN .* ...
        model_layer12_equations('Phi_c', perChanData.RVN.uE(:), f_max_c) .* ...
        model_layer12_equations('h',     perChanData.RVN.uE(:), u_half);
    resLVN = perChanData.LVN.rate(:) - W_c_LVN .* ...
        model_layer12_equations('Phi_c', perChanData.LVN.uE(:), f_max_c) .* ...
        model_layer12_equations('h',     perChanData.LVN.uE(:), u_half);
    res = [resRVN; resLVN];
end

function [bestTheta, bestResnorm] = fit_one_channel_only(chanData, scaleGuess, nRestarts, opts)
% Original single-channel 3-param [W_c, u_half, f_max_c] fit -- retained ONLY
% as the L12 validation comparison, not the live Stage 1 procedure.
    rate = chanData.rate; uEsel = chanData.uE;
    resFun = @(th) rate(:) - th(1) .* model_layer12_equations('Phi_c', uEsel(:), th(3)) ...
                                    .* model_layer12_equations('h',     uEsel(:), th(2));
    lb = [0 0 0]; ub = [Inf Inf Inf];
    results = repmat(struct('theta',[nan nan nan],'resnorm',NaN), nRestarts, 1);
    for r = 1:nRestarts
        x0 = [ scaleGuess * 10.^(rand()*4-2), 10.^(rand()*4), 10.^(rand()*4) ];
        try
            [th, resnorm] = lsqnonlin(resFun, x0, lb, ub, opts);
        catch
            th = [NaN NaN NaN]; resnorm = NaN;
        end
        results(r) = struct('theta',th,'resnorm',resnorm);
    end
    resnorms = arrayfun(@(s) s.resnorm, results);
    [bestResnorm, bestIdx] = min(resnorms);
    bestTheta = results(bestIdx).theta;
end

function [M,E] = parse_me_local(cond)
% Replicates parse_me in bulk_mixed_models.m exactly: M(\d+) / E(\d+),
% default 0 if the respective token is absent (so 'E100' -> M=0,E=100).
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end
