function [thetaFull, resnorm2] = fit_layer12_stage2_interaction(rawRows, channelLabel, candidate)
%FIT_LAYER12_STAGE2_INTERACTION Fit [f_lo, f_hi, u_M_half, f_max_c2, w_1, w_2] to
%   the 9 combined M x E conditions for ONE channel, using ONE Stage 1
%   CANDIDATE (either the raw or swapped member of thetaE.candidates from
%   fit_layer12_stage1_electrical.m -- this function does not know or care
%   which; it is called once per candidate by run_layer12_first_pass.m,
%   which does the actual basin selection, L15).
%
%   candidate must have fields: u_half, f_max_c, W_c_RVN, W_c_LVN (i.e. one
%   element of thetaE.candidates(:)) -- or, for a single-group caller like
%   run_layer12_total_channel.m, u_half/f_max_c plus W_c_<channelLabel>.
%
%   [thetaFull, resnorm2] = fit_layer12_stage2_interaction(rawRows, 'RVN', thetaE.candidates(1))
%
%   resnorm2 (2nd output): the fitted Stage 2 residual (best resnorm across
%   restarts) -- required by run_layer12_first_pass.m's basin selection
%   logic (L15), which compares resnorm2 across candidates and channels. Do
%   not silently drop this second output when calling from new code.
%
%   L16 REPARAMETRIZATION (fix by construction, not by search): (f_lo, f_hi, w_1)
%   has the SAME exact swap symmetry as Stage 1's (u_half, f_max_c, W_c) (L12) --
%   g1(u_M; f_hi, f_lo) = (f_lo/f_hi) * g1(u_M; f_lo, f_hi) for every u_M (proof:
%   model_fitting_plan.md L16), so swapping f_lo<->f_hi and rescaling
%   w_1 -> w_1*(f_hi/f_lo) leaves term_1 exactly unchanged. UNLIKE (u_half,f_max_c),
%   there is no term_2-style analog here to break the tie (L13's mechanism doesn't
%   apply), and Stage 2 isn't pooled across channels the way Stage 1 is -- so
%   instead of fitting f_lo/f_hi directly and relying on a tie-break that doesn't
%   exist, the optimizer fits rho1=sqrt(f_lo*f_hi) (peak location, exactly
%   invariant) and kappa1=0.5*ln(f_hi/f_lo) (SIGN-CONSTRAINED >=0 -- a fixed
%   REPORTING CONVENTION, not a data-driven finding) directly. f_lo/f_hi are
%   recovered afterward as DERIVED display quantities in thetaFull, for
%   backward-compatible use by model_layer12_equations.m and every existing
%   caller -- but rho1/kappa1 (also stored in thetaFull) are the actual
%   identifiable fit outputs; do not report f_lo/f_hi individually as if they
%   were data-driven. Stage 1's (u_half,f_max_c) deliberately keeps its existing
%   L13/L15 tie-break approach rather than this same reparametrization -- that
%   tie-break is empirically validated (decisive on real data) and NOT wrong,
%   just a different (and, for that specific pair, viable) resolution mechanism;
%   see model_fitting_plan.md L16 for this documented design choice.
%
%   Objective (model_fitting_plan.md Section 5.3, step 2): nonlinear least
%   squares over ALL matching trials (one residual per trial, i.e. summed
%   over animals AND conditions -- repeat trials of the same animal/condition,
%   if present, each contribute their own residual rather than being
%   pre-averaged, consistent with Stage 1's convention):
%       sum( (observed_rate - rate_hat(u_M,u_E; theta)).^2 )
%   over the 9 combined (u_M>0 & u_E>0) conditions.

    if nargin < 3 || isempty(candidate)
        error('fit_layer12_stage2_interaction:args','candidate (one element of thetaE.candidates) is required.');
    end
    req = {'u_half','f_max_c'};
    for k = 1:numel(req)
        if ~isfield(candidate, req{k})
            error('fit_layer12_stage2_interaction:candidate','candidate missing field %s.', req{k});
        end
    end
    % Dynamic field lookup (not a hardcoded RVN/LVN switch) so this works for
    % any channel/group label Stage 1 produced a W_c_<LABEL> for -- e.g. the
    % 'TOTAL' pseudo-channel (run_layer12_total_channel.m), not just RVN/LVN.
    wcField = sprintf('W_c_%s', upper(channelLabel));
    if ~isfield(candidate, wcField)
        error('fit_layer12_stage2_interaction:badChannel', ...
            'candidate has no field %s for channelLabel ''%s''.', wcField, channelLabel);
    end
    W_c_this = candidate.(wcField);

    % ---- filter rows: recovery phase, requested channel, combined M x E ----
    isRec = strcmpi({rawRows.phase}, 'recovery');
    isCh  = strcmpi({rawRows.label}, channelLabel);
    uM = nan(numel(rawRows),1); uE = nan(numel(rawRows),1);
    for i = 1:numel(rawRows)
        [uM(i), uE(i)] = parse_me_local(rawRows(i).condition);
    end
    isCombined = (uM > 0) & (uE > 0);
    sel = find(isRec(:) & isCh(:) & isCombined(:));

    rate = arrayfun(@(i) rawRows(i).mean.rate, sel);
    uMsel = uM(sel); uEsel = uE(sel);
    finite = isfinite(rate) & isfinite(uMsel) & isfinite(uEsel);
    sel = sel(finite); rate = rate(finite); uMsel = uMsel(finite); uEsel = uEsel(finite);

    uCells = unique([uMsel(:) uEsel(:)], 'rows');
    fprintf('[stage2-%s-%s] %d combined-condition trials across %d distinct M x E cells:\n', ...
        channelLabel, candidate_label_local(candidate), numel(sel), size(uCells,1));
    for c = 1:size(uCells,1)
        n = nnz(uMsel==uCells(c,1) & uEsel==uCells(c,2));
        fprintf('    M%d x E%d : %d trial(s)\n', uCells(c,1), uCells(c,2), n);
    end
    if size(uCells,1) < 9
        warning('fit_layer12_stage2_interaction:missingCells', ...
            '%s: only %d/9 M x E cells have any trials -- reporting per-cell counts above, not silently excluding.', ...
            channelLabel, size(uCells,1));
    end

    nParams = 6;
    if numel(sel) < nParams
        error('fit_layer12_stage2_interaction:tooFew', ...
            '%s: only %d combined-condition trials with finite rate -- need >= %d to fit 6 params.', ...
            channelLabel, numel(sel), nParams);
    end

    % ---- objective: theta = [rho1, kappa1, u_M_half, f_max_c2, w_1, w_2] ---
    % (L16: rho1=sqrt(f_lo*f_hi), kappa1=0.5*ln(f_hi/f_lo) >= 0 -- see header.)
    resFun = @(th) rate(:) - predict_local(th, uMsel(:), uEsel(:), candidate, W_c_this);

    % bounds: all params > 0, EXCEPT kappa1 which is >= 0 by the L16 sign
    % convention (fixes f_hi >= f_lo by construction, not by data -- the
    % mirror-image family with f_hi < f_lo is the SAME shape with w_1 rescaled,
    % per the exact symmetry, so nothing is lost by excluding it). f_max_c2
    % constrained STRICTLY below the fixed f_max_c (enforced via bounds, not
    % just checked post-fit) so f1 and f2 cannot converge to the same corner
    % and collapse into proportional copies of each other (model_fitting_plan.md
    % Section 3 / Section 5, "MUST be constrained distinct from f_max_c"). This
    % f_max_c2 check is UNRELATED to the L16 kappa1 constraint -- do not conflate.
    fMaxC2Cap = 0.999 * candidate.f_max_c;
    KAPPA1_MAX = 10;   % generous cap (f_hi/f_lo up to exp(20) ~ 4.85e8); a hard
                        % safety bound only, not expected to bind in practice
    lb = [0    0           0    0            0    0   ];
    ub = [Inf  KAPPA1_MAX  Inf  fMaxC2Cap    Inf  Inf ];
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'FunctionTolerance',1e-12, 'StepTolerance',1e-12, 'MaxFunctionEvaluations',5000);

    nRestarts = 20;
    results = repmat(struct('theta',nan(1,6),'resnorm',NaN,'exitflag',NaN), nRestarts, 1);
    for r = 1:nRestarts
        x0 = [ 10.^(rand()*3), ...                  % rho1      ~ [1, 1e3]   (peak location)
               rand()*4, ...                         % kappa1    ~ [0, 4]     (ratio up to e^8 ~ 2981)
               10.^(rand()*3), ...                  % u_M_half  ~ [1, 1e3]
               rand() * fMaxC2Cap, ...               % f_max_c2  ~ [0, cap]
               10.^(rand()*4-1), ...                 % w_1       ~ [0.1, 1e3] -- unknown scale, kept wide
               10.^(rand()*4-1) ];                    % w_2       ~ [0.1, 1e3] -- unknown scale, kept wide
        try
            [th, resnorm, ~, exitflag] = lsqnonlin(resFun, x0, lb, ub, opts);
        catch ME
            warning('fit_layer12_stage2_interaction:restart','restart %d failed (%s)', r, ME.message);
            th = nan(1,6); resnorm = NaN; exitflag = -99;
        end
        results(r) = struct('theta',th,'resnorm',resnorm,'exitflag',exitflag);
    end

    resnorms = arrayfun(@(s) s.resnorm, results);
    [bestResnorm, bestIdx] = min(resnorms);
    bt = results(bestIdx).theta;
    resnorm2 = bestResnorm;
    rho1 = bt(1); kappa1 = bt(2);
    f_lo = rho1 * exp(-kappa1); f_hi = rho1 * exp(kappa1);   % DERIVED display quantities (L16)

    fprintf('[stage2-%s-%s] %d/%d restarts converged (exitflag>0).\n', channelLabel, candidate_label_local(candidate), ...
        nnz(arrayfun(@(s) s.exitflag>0, results)), nRestarts);
    fprintf('[stage2-%s-%s] best resnorm=%.6g, exitflag=%d\n', channelLabel, candidate_label_local(candidate), ...
        bestResnorm, results(bestIdx).exitflag);
    fprintf('[stage2-%s-%s] IDENTIFIABLE (L16): rho1=sqrt(f_lo*f_hi)=%.4g, kappa1=0.5*ln(f_hi/f_lo)=%.4g (sign-constrained convention)\n', ...
        channelLabel, candidate_label_local(candidate), rho1, kappa1);
    fprintf('[stage2-%s-%s] DERIVED display values: f_lo=%.4g, f_hi=%.4g -- theta=[u_M_half=%.4g, f_max_c2=%.4g, w_1=%.4g, w_2=%.4g]\n', ...
        channelLabel, candidate_label_local(candidate), f_lo, f_hi, bt(3), bt(4), bt(5), bt(6));

    near = resnorms <= bestResnorm * 1.01;
    thetaNear = vertcat(results(near).theta);
    names = {'rho1','kappa1','u_M_half','f_max_c2','w_1','w_2'};
    fprintf('[stage2-%s-%s] spread among %d near-best restarts:\n', channelLabel, candidate_label_local(candidate), nnz(near));
    for k = 1:6
        fprintf('    %-10s [%.4g, %.4g]\n', names{k}, min(thetaNear(:,k)), max(thetaNear(:,k)));
    end

    % thetaFull must have a scalar W_c field for model_layer12_equations to use
    % -- built here by copying this channel's W_c_this in as the single 'W_c'
    % field name model_layer12_equations expects. f_lo/f_hi included for
    % backward-compatible use by model_layer12_equations.m and every existing
    % caller; rho1/kappa1 (L16) are the actual identifiable fit outputs.
    thetaFull = struct('W_c',W_c_this, 'u_half',candidate.u_half, 'f_max_c',candidate.f_max_c, ...
        'f_lo',f_lo, 'f_hi',f_hi, 'rho1',rho1, 'kappa1',kappa1, ...
        'u_M_half',bt(3), 'f_max_c2',bt(4), 'w_1',bt(5), 'w_2',bt(6), ...
        'resnorm',bestResnorm, 'nTrials',numel(sel), 'nRestartsConverged',nnz(near));
    thetaFull.allRestarts = results;   % diagnostic: full per-restart theta/resnorm/exitflag table

    % ---- f_max_c2 distinctness check (explicit report, not just a bound) --
    sep = candidate.f_max_c - thetaFull.f_max_c2;
    fprintf('[stage2-%s-%s] f_max_c2 distinctness: f_max_c=%.4g, f_max_c2=%.4g, separation=%.4g (%s)\n', ...
        channelLabel, candidate_label_local(candidate), candidate.f_max_c, thetaFull.f_max_c2, sep, ...
        string(sep > 1e-6));
    if sep <= 1e-6
        warning('fit_layer12_stage2_interaction:notDistinct', ...
            ['%s: f_max_c2 converged to (near-)the same value as f_max_c -- f1 and f2 are effectively ' ...
             'proportional, and the model has lost the moving-optimum property motivating two mechanical ' ...
             'coupling terms (model_fitting_plan.md Section 3). Reporting, not silently accepting.'], channelLabel);
    end

    % ---- L17: (u_M_half, w_2) asymptotic ridge flag (post-fit diagnostic
    % ONLY -- no optimizer/bound change; a genuinely large true u_M_half is a
    % real finding, not an error to constrain away, see model_fitting_plan.md
    % L17). When u_M_half lands far beyond the tested u_M range, g2 becomes
    % locally linear and only kappa2=w_2/u_M_half is identified, not the two
    % individually -- this is a DATA-COVERAGE limitation (L4: coarse 3-point
    % u_M sweep), NOT an exact algebraic symmetry like L12/L16, so there is no
    % reparametrization that removes it; monitor and report, don't suppress.
    maxUM_tested = max(uMsel);   % from THIS fit's actual data, not hardcoded
    if thetaFull.u_M_half > 5 * maxUM_tested
        kappa2 = thetaFull.w_2 / thetaFull.u_M_half;
        warning('fit_layer12_stage2_interaction:asymptoticRidge', ...
            ['[%s-%s] u_M_half = %.1f is > 5x max tested u_M (%.0f) -- L17 asymptotic ' ...
             'ridge risk. u_M_half and w_2 individually unreliable in this regime; ' ...
             'report kappa2 = w_2/u_M_half = %.4g as the more robust quantity ' ...
             'alongside (not instead of) the raw values.'], channelLabel, candidate_label_local(candidate), ...
            thetaFull.u_M_half, maxUM_tested, kappa2);
        thetaFull.asymptoticRidgeFlag = true;
        thetaFull.kappa2 = kappa2;
    else
        thetaFull.asymptoticRidgeFlag = false;
        thetaFull.kappa2 = thetaFull.w_2 / thetaFull.u_M_half;   % still stored -- harmless/robust either way
    end

    % ---- REQUIRED pass/fail pattern check (model_fitting_plan.md Sec 5.5 /
    % matlab_implementation_instructions.md Section 5) -- print explicitly,
    % do not bury, and return the fit regardless of pass/fail. -------------
    pred_M50E100  = model_layer12_equations('r_vagus', 50, 100, thetaFull);
    pred_M100E100 = model_layer12_equations('r_vagus', 100, 100, thetaFull);
    pred_M100E10  = model_layer12_equations('r_vagus', 100, 10, thetaFull);
    pred_M50E10   = model_layer12_equations('r_vagus', 50, 10, thetaFull);
    check1 = pred_M50E100 > pred_M100E100;
    check2 = pred_M100E10 > pred_M50E10;
    fprintf('[%s-%s] Check 1 (M50xE100 > M100xE100): %s (%.3f vs %.3f)\n', channelLabel, candidate_label_local(candidate), ...
        string(check1), pred_M50E100, pred_M100E100);
    fprintf('[%s-%s] Check 2 (M100xE10 > M50xE10):   %s (%.3f vs %.3f)\n', channelLabel, candidate_label_local(candidate), ...
        string(check2), pred_M100E10, pred_M50E10);
    if ~check1 || ~check2
        fprintf(2, ['[stage2-%s-%s] *** PATTERN CHECK FAILED *** -- the two-term structure''s ' ...
            'one required job (moving optimum) was NOT reproduced by this fit. Reporting prominently, ' ...
            'NOT retrying with different starting points to force a pass.\n'], channelLabel, candidate_label_local(candidate));
    end
    thetaFull.check1_M50E100_gt_M100E100 = check1;
    thetaFull.check2_M100E10_gt_M50E10   = check2;
end

% ----------------------------------------------------------------------
function y = predict_local(th, uMv, uEv, candidate, W_c_this)
% th = [rho1, kappa1, u_M_half, f_max_c2, w_1, w_2] (L16) -- convert rho1/kappa1
% to f_lo/f_hi here, at the one place the optimizer's search vector meets
% model_layer12_equations.m (which is unchanged and still expects f_lo/f_hi).
    rho1 = th(1); kappa1 = th(2);
    f_lo = rho1 .* exp(-kappa1); f_hi = rho1 .* exp(kappa1);
    p = struct('W_c',W_c_this, 'u_half',candidate.u_half, 'f_max_c',candidate.f_max_c, ...
        'f_lo',f_lo, 'f_hi',f_hi, 'u_M_half',th(3), 'f_max_c2',th(4), 'w_1',th(5), 'w_2',th(6));
    y = model_layer12_equations('r_vagus', uMv, uEv, p);
end

function s = candidate_label_local(candidate)
% Cosmetic only (console labeling) -- falls back to 'n/a' for callers (e.g.
% run_layer12_total_channel.m's single-group candidate) that don't set .label.
    if isfield(candidate,'label') && ~isempty(candidate.label); s = candidate.label; else; s = 'n/a'; end
end

function [M,E] = parse_me_local(cond)
% Replicates parse_me in bulk_mixed_models.m exactly: M(\d+) / E(\d+),
% default 0 if the respective token is absent.
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end
