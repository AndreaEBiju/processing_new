%TEST_LAYER12_SYNTHETIC_SMOKE Synthetic-data smoke test for the Layer 1-2
%   first-pass pipeline (model_layer12_equations.m, fit_layer12_stage1/2,
%   plot_layer12_fit_diagnostics.m, run_layer12_first_pass.m).
%
%   NOT one of the 6 files in matlab_implementation_instructions.md Section 2
%   -- this is throwaway test infrastructure to exercise the whole pipeline
%   end-to-end BEFORE real data is available, since no real rawRows exists
%   yet. Delete/ignore once real data arrives; run_layer12_first_pass.m is
%   the actual deliverable this validates.
%
%   Builds a synthetic rawRows struct array with the EXACT shape produced by
%   run_pipeline_bulk.m's harvest_channel/empty_row (see
%   matlab_implementation_instructions.md Section 0.5), generated from a
%   KNOWN ground-truth theta with small noise, then runs
%   run_layer12_first_pass on it and checks the fitted parameters land close
%   to the ground truth and the Section 5 pattern checks come out as
%   expected for the chosen ground truth.
%
%   L12 UPDATE: ground truth now gives RVN and LVN a DELIBERATELY DIFFERENT
%   true W_c (2 vs 3.5) while SHARING u_half/f_max_c and every Stage 2 shape
%   param -- this is exactly the scenario Stage 1's pooled fit
%   (fit_layer12_stage1_electrical.m) is supposed to handle correctly:
%   recover the one shared (u_half,f_max_c) pair AND each channel's distinct
%   W_c, rather than landing in two different (u_half,f_max_c) basins the way
%   the old independent per-channel fits empirically did on this same kind of
%   data (see model_fitting_plan.md L12).
%
%   Ground truth's shape params (f_lo,f_hi,u_M_half,f_max_c2,w_1,w_2) were
%   chosen by a random search over the model's own parameter space for the
%   LARGEST achievable relative margin on both hotspot checks simultaneously
%   -- the best found across 500k random draws was only ~11% relative
%   margin. That margin is a property of the model structure itself (worth
%   flagging: even at its best-case parameter choice, this two-term
%   structure only barely threads the M50xE100 / M100xE10 needle --
%   consistent with L9's "near-zero residual df" warning).
%   Noise here is deliberately kept small (2%) so this smoke test validates
%   PIPELINE MECHANICS, not robustness to noise -- real data may legitimately
%   fail the pattern check even with a correctly-implemented pipeline.

clear; clc;
rng(42);   % reproducible synthetic data

thetaTrueShared = struct('u_half',84, 'f_max_c',1850);           % SHARED, pooled across channels (L12)
thetaTrueOther  = struct('f_lo',1.15, 'f_hi',32.3, 'u_M_half',250, ...
                          'f_max_c2',14.8, 'w_1',7.3, 'w_2',27.4); % per-channel Stage 2 shape (same true value both channels here)
thetaTrueWc     = struct('RVN',2, 'LVN',3.5);                     % DELIBERATELY DIFFERENT per channel (L12)

animals  = {'a1','a2','a3','a4','a5','a6'};
channels = {'RVN','LVN'};
eAloneLevels = [10 100 1000];
combinedM = [10 50 100]; combinedE = [10 100 1000];
noiseFracRecovery = 0.02;    % 2% relative noise on recovery-phase mean rate
baselineRateMean  = 0.6;     % Hz, arbitrary (baseline is not part of the r_vagus fit target)
noiseFracBaseline = 0.10;
nSampPerTrial = 20;

rows = repmat(empty_row_local(), 0, 1);
for a = 1:numel(animals)
    animal = animals{a};
    for c = 1:numel(channels)
        label = channels{c};
        thetaTrueCh = full_theta_local(label, thetaTrueShared, thetaTrueOther, thetaTrueWc);
        % M-alone (u_E=0): NOT part of the r_vagus fit target (L1 -- excluded
        % by construction, r_vagus(u_M,0)==0 exactly), but bulk_mixed_models'
        % saturated cross-check model NEEDS these rows: without them, m_k
        % (M-presence indicator) is a perfect linear combination of that row's
        % M x E interaction columns for every remaining row (E-alone rows have
        % m_k=0 identically, combined rows have m_k == sum_l(m_k:e_l) exactly),
        % making the fixed-effects design matrix rank-deficient. Recovery rate
        % here is generated as "flat vs baseline" (matching the model's own
        % "observed flat mechanical-alone data" assumption), NOT via r_vagus.
        for m = combinedM
            cond = sprintf('M%d', m);
            rows = add_trial_local(rows, animal, cond, label, baselineRateMean, noiseFracBaseline, ...
                baselineRateMean, noiseFracBaseline, nSampPerTrial);
        end
        for e = eAloneLevels
            cond = sprintf('E%d', e);
            trueRate = model_layer12_equations('r_vagus', 0, e, thetaTrueCh);
            rows = add_trial_local(rows, animal, cond, label, trueRate, noiseFracRecovery, ...
                baselineRateMean, noiseFracBaseline, nSampPerTrial);
        end
        for m = combinedM
            for e = combinedE
                cond = sprintf('M%dE%d', m, e);
                trueRate = model_layer12_equations('r_vagus', m, e, thetaTrueCh);
                rows = add_trial_local(rows, animal, cond, label, trueRate, noiseFracRecovery, ...
                    baselineRateMean, noiseFracBaseline, nSampPerTrial);
            end
        end
    end
end

nCondPerChannel = numel(combinedM) + numel(eAloneLevels) + numel(combinedM)*numel(combinedE);
fprintf('Synthetic rawRows built: %d rows (%d animals x %d channels x %d conditions x 2 phases)\n', ...
    numel(rows), numel(animals), numel(channels), nCondPerChannel);
fprintf('True W_c: RVN=%.4g, LVN=%.4g (deliberately different); shared true u_half=%.4g, f_max_c=%.4g\n\n', ...
    thetaTrueWc.RVN, thetaTrueWc.LVN, thetaTrueShared.u_half, thetaTrueShared.f_max_c);

results = run_layer12_first_pass(rows, '');

% ---- validation against ground truth ----------------------------------
fprintf('\n==== Ground-truth recovery check (synthetic data only) ====\n');
fprintf(['NOTE: post-L16, f_lo/f_hi/w_1 recovery should be much closer to ground truth than earlier\n' ...
    'runs showed (that was the (f_lo,f_hi,w_1) exact swap symmetry, model_fitting_plan.md L16 --\n' ...
    'now fixed by reparametrizing to rho1=sqrt(f_lo*f_hi)/kappa1 instead of fitting f_lo/f_hi\n' ...
    'directly). Residual %%err here should reflect genuine noise + parameter correlation, not the\n' ...
    'old bimodal-restart symptom -- compare against rho1/kappa1 below, and check that near-best\n' ...
    'restarts converge tightly (printed by fit_layer12_stage2_interaction.m) rather than spanning\n' ...
    'orders of magnitude. Judge by resnorm + pattern checks first regardless.\n\n']);

fprintf('-- L15 basin selection outcome --\n');
if isempty(results.thetaE)
    fprintf('   Stage 1 (pooled) FAILED.\n');
else
    fprintf('   selected=%s  rejected=%s  margin=%.1f%%  ambiguous=%d  channelDisagreement=%d\n', ...
        results.thetaE_selected.label, results.thetaE_rejected.label, ...
        results.selectionMarginPct, results.basinAmbiguous, results.channelDisagreement);
end

fprintf('-- Stage 1 SHARED params (u_half, f_max_c; pooled across RVN+LVN, SELECTED candidate) --\n');
if isempty(results.thetaE)
    fprintf('   Stage 1 (pooled) FAILED.\n');
else
    % NOTE: the selected candidate may be either the "raw" or "swapped" basin
    % (L12) -- comparing its raw u_half/f_max_c against ground truth directly
    % can show a large %%err even on a CORRECT fit if the swapped basin won.
    % rho=sqrt(u_half*f_max_c) and the W_c ratio are the values that are
    % actually meaningful to compare (both identical for either basin).
    sharedNames = {'u_half','f_max_c'};
    for p = 1:numel(sharedNames)
        truthVal = thetaTrueShared.(sharedNames{p});
        fitVal = results.thetaE_selected.(sharedNames{p});
        pctErr = 100*abs(fitVal-truthVal)/abs(truthVal);
        fprintf('   %-10s true=%10.4g  fit=%10.4g  %%err=%6.1f%%  (raw value -- see rho/ratio below for the meaningful check)\n', ...
            sharedNames{p}, truthVal, fitVal, pctErr);
    end
    fprintf('   W_c_RVN    true=%10.4g  fit=%10.4g\n', thetaTrueWc.RVN, results.thetaE_selected.W_c_RVN);
    fprintf('   W_c_LVN    true=%10.4g  fit=%10.4g\n', thetaTrueWc.LVN, results.thetaE_selected.W_c_LVN);
    rhoTrue = sqrt(thetaTrueShared.u_half * thetaTrueShared.f_max_c);
    ratioTrue = thetaTrueWc.LVN / thetaTrueWc.RVN;
    fprintf('   rho        true=%10.4g  fit=%10.4g  (the actually-identifiable invariant, L12)\n', rhoTrue, results.thetaE.rho);
    fprintf('   W_c ratio  true=%10.4g  fit=%10.4g  (the actually-identifiable invariant, L12)\n', ratioTrue, results.thetaE.W_c_ratio);
    fprintf('   L12 validation: pooled resnorm=%.6g, independent resnorm=%.6g, %+.1f%%\n', ...
        results.thetaE.resnormPooled, results.thetaE.resnormIndependentTotal, results.thetaE.pooledVsIndependentPctChange);
end

perChanParamNames = {'f_lo','f_hi','u_M_half','f_max_c2','w_1','w_2'};
for c = 1:numel(channels)
    ch = channels{c};
    entry = results.perChannelDetail.(ch);
    if isempty(entry.thetaFull)
        fprintf('\n%s: FIT FAILED -- %s\n', ch, entry.error);
        continue;
    end
    fprintf('\n-- %s Stage 2 (per-channel) --\n', ch);
    fprintf('   %-10s true=%10.4g  fit=%10.4g\n', 'W_c', thetaTrueWc.(ch), entry.thetaFull.W_c);
    rho1True = sqrt(thetaTrueOther.f_lo * thetaTrueOther.f_hi);
    kappa1True = abs(0.5*log(thetaTrueOther.f_hi/thetaTrueOther.f_lo));   % thetaTrue's own sign is arbitrary; fit enforces kappa1>=0 (L16)
    fprintf('   %-10s true=%10.4g  fit=%10.4g  (L16 identifiable quantity)\n', 'rho1', rho1True, entry.thetaFull.rho1);
    fprintf('   %-10s true=%10.4g  fit=%10.4g  (L16 identifiable quantity, sign-constrained)\n', 'kappa1', kappa1True, entry.thetaFull.kappa1);
    for p = 1:numel(perChanParamNames)
        truthVal = thetaTrueOther.(perChanParamNames{p});
        fitVal = entry.thetaFull.(perChanParamNames{p});
        pctErr = 100*abs(fitVal-truthVal)/abs(truthVal);
        flag = ''; if pctErr > 30; flag = '  <-- check resnorm/restart-spread before treating as a fitting problem'; end
        fprintf('   %-10s true=%10.4g  fit=%10.4g  %%err=%6.1f%%%s\n', ...
            perChanParamNames{p}, truthVal, fitVal, pctErr, flag);
    end
    fprintf('   Check1 (M50xE100>M100xE100): %s\n', mat2str(entry.thetaFull.check1_M50E100_gt_M100E100));
    fprintf('   Check2 (M100xE10>M50xE10):   %s\n', mat2str(entry.thetaFull.check2_M100E10_gt_M50E10));
end

% ---- ground-truth-known expectation (from the search above) -----------
% Shape params are shared across channels in this synthetic truth, so the
% pattern check's ground-truth answer is the same for both -- compute once
% using RVN's full true theta.
thetaTrueRVN = full_theta_local('RVN', thetaTrueShared, thetaTrueOther, thetaTrueWc);
trueCheck1 = model_layer12_equations('r_vagus',50,100,thetaTrueRVN) > model_layer12_equations('r_vagus',100,100,thetaTrueRVN);
trueCheck2 = model_layer12_equations('r_vagus',100,10,thetaTrueRVN) > model_layer12_equations('r_vagus',50,10,thetaTrueRVN);
fprintf('\nGround truth itself satisfies: Check1=%s, Check2=%s (both should be true by construction)\n', ...
    mat2str(trueCheck1), mat2str(trueCheck2));

% ==== run_layer12_total_channel.m validation (same synthetic rows) ========
% Ground truth for this exact synthetic dataset (RVN/LVN share every shape
% param except W_c) has a closed form for the TOTAL series: since r_vagus is
% linear in (W_c, w_1, w_2) for fixed u_half/f_max_c/shape, sum(RVN,LVN) is
% EXACTLY equal to one r_vagus with W_c_total=W_c_RVN+W_c_LVN, w_1*2, w_2*2,
% same u_half/f_max_c/f_lo/f_hi/u_M_half/f_max_c2 -- verified numerically in
% Python (see conversation) to floating-point precision before trusting this.
fprintf('\n\n==== run_layer12_total_channel.m validation (same synthetic rows) ====\n');
thetaTrueTotal = thetaTrueOther;
thetaTrueTotal.u_half = thetaTrueShared.u_half;
thetaTrueTotal.f_max_c = thetaTrueShared.f_max_c;
thetaTrueTotal.W_c = thetaTrueWc.RVN + thetaTrueWc.LVN;
thetaTrueTotal.w_1 = 2 * thetaTrueOther.w_1;
thetaTrueTotal.w_2 = 2 * thetaTrueOther.w_2;

resultsTotal = run_layer12_total_channel(rows, '');

fprintf('\n-- TOTAL ground-truth recovery check --\n');
if isempty(resultsTotal.thetaFull)
    fprintf('   TOTAL fit FAILED.\n');
else
    fprintf('   %-10s true=%10.4g  fit=%10.4g\n', 'rho(shared)', ...
        sqrt(thetaTrueShared.u_half*thetaTrueShared.f_max_c), resultsTotal.thetaE.rho);
    rho1TrueTotal = sqrt(thetaTrueOther.f_lo * thetaTrueOther.f_hi);
    kappa1TrueTotal = abs(0.5*log(thetaTrueOther.f_hi/thetaTrueOther.f_lo));
    fprintf('   %-10s true=%10.4g  fit=%10.4g  (L16 identifiable quantity)\n', 'rho1', rho1TrueTotal, resultsTotal.thetaFull.rho1);
    fprintf('   %-10s true=%10.4g  fit=%10.4g  (L16 identifiable quantity, sign-constrained)\n', 'kappa1', kappa1TrueTotal, resultsTotal.thetaFull.kappa1);
    totalParamNames = {'W_c','f_lo','f_hi','u_M_half','f_max_c2','w_1','w_2'};
    for p = 1:numel(totalParamNames)
        truthVal = thetaTrueTotal.(totalParamNames{p});
        fitVal = resultsTotal.thetaFull.(totalParamNames{p});
        pctErr = 100*abs(fitVal-truthVal)/abs(truthVal);
        flag = ''; if pctErr > 30; flag = '  <-- check resnorm/restart-spread before treating as a fitting problem'; end
        fprintf('   %-10s true=%10.4g  fit=%10.4g  %%err=%6.1f%%%s\n', ...
            totalParamNames{p}, truthVal, fitVal, pctErr, flag);
    end
    fprintf('   Check1 (M50xE100>M100xE100): %s\n', mat2str(resultsTotal.thetaFull.check1_M50E100_gt_M100E100));
    fprintf('   Check2 (M100xE10>M50xE10):   %s\n', mat2str(resultsTotal.thetaFull.check2_M100E10_gt_M50E10));
    if ~isempty(resultsTotal.thetaFullSwapped)
        fprintf('   Dual-basin resnorm: normal=%.6g  swapped=%.6g\n', ...
            resultsTotal.thetaFull.resnorm, resultsTotal.thetaFullSwapped.resnorm);
    end
end
trueCheck1Total = model_layer12_equations('r_vagus',50,100,thetaTrueTotal) > model_layer12_equations('r_vagus',100,100,thetaTrueTotal);
trueCheck2Total = model_layer12_equations('r_vagus',100,10,thetaTrueTotal) > model_layer12_equations('r_vagus',50,10,thetaTrueTotal);
fprintf('   Ground truth (TOTAL) itself satisfies: Check1=%s, Check2=%s\n', ...
    mat2str(trueCheck1Total), mat2str(trueCheck2Total));

% ----------------------------------------------------------------------
function p = full_theta_local(channelLabel, shared, other, wc)
% Composes a full 9-field r_vagus param struct for one channel from the
% SHARED (u_half,f_max_c), the per-channel-but-same-here shape params, and
% that channel's own (deliberately distinct) W_c.
    p = other;
    p.u_half = shared.u_half;
    p.f_max_c = shared.f_max_c;
    p.W_c = wc.(channelLabel);
end

function r = empty_row_local()
% Matches run_pipeline_bulk.m's empty_row() EXACTLY -- this is the real
% shape fit_layer12_stage1/2 and plot_layer12_fit_diagnostics.m expect.
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'rate_t',[],'vpp',[],'fwhm',[],'spk_t',[],'cv2',[],'cv2_t',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN, 'stem','');
end

function rows = add_trial_local(rows, animal, cond, label, trueRate, noiseFracRec, ...
        baselineMean, noiseFracBase, nSamp)
% Appends ONE baseline row + ONE recovery row (sharing a stem, per
% bulk_compile.m's trial-pairing convention) for (animal, cond, label).
    stem = sprintf('%s_%s', animal, cond);

    recRate = max(trueRate * (1 + noiseFracRec*randn()), 1e-3);
    rRec = empty_row_local();
    rRec.animal = animal; rRec.condition = cond; rRec.phase = 'recovery';
    rRec.label = label; rRec.stem = stem;
    rRec.dist.rate = max(recRate*(1+noiseFracRec*randn(nSamp,1)), 0);
    rRec.dist.rate_t = linspace(0,60,nSamp)';
    rRec.dist.vpp  = 50 + 5*randn(nSamp,1);
    rRec.dist.fwhm = 0.5 + 0.05*randn(nSamp,1);
    rRec.dist.cv2  = abs(0.8 + 0.1*randn(nSamp,1));
    rRec.mean.rate = recRate;
    rRec.mean.excess = 10 + randn();
    rRec.mean.vpp  = median(rRec.dist.vpp);
    rRec.mean.fwhm = median(rRec.dist.fwhm);
    rRec.mean.cv2  = mean(rRec.dist.cv2);
    rRec.fanoSlope = -0.10 + 0.05*randn();
    rows(end+1) = rRec; %#ok<AGROW>

    baseRate = max(baselineMean*(1+noiseFracBase*randn()), 1e-3);
    rBase = empty_row_local();
    rBase.animal = animal; rBase.condition = cond; rBase.phase = 'baseline';
    rBase.label = label; rBase.stem = stem;
    rBase.dist.rate = max(baseRate*(1+noiseFracBase*randn(nSamp,1)), 0);
    rBase.dist.rate_t = linspace(0,60,nSamp)';
    rBase.dist.vpp  = 48 + 5*randn(nSamp,1);
    rBase.dist.fwhm = 0.48 + 0.05*randn(nSamp,1);
    rBase.dist.cv2  = abs(0.85 + 0.1*randn(nSamp,1));
    rBase.mean.rate = baseRate;
    rBase.mean.excess = 9 + randn();
    rBase.mean.vpp  = median(rBase.dist.vpp);
    rBase.mean.fwhm = median(rBase.dist.fwhm);
    rBase.mean.cv2  = mean(rBase.dist.cv2);
    rBase.fanoSlope = -0.12 + 0.05*randn();
    rows(end+1) = rBase; %#ok<AGROW>
end
