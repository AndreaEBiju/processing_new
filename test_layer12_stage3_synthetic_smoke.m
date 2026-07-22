rng('shuffle');   % fresh matlab -batch sessions otherwise reuse the SAME default seed
% Synthetic ground-truth check for fit_layer12_stage3_shape.m's FULLY
% decoupled architecture: v_1/v_2 (not w_1/w_2) AND fresh rho_1/kappa_1 (not
% frozen from Stage 2, which was shown unstable given w_1~0). This
% specifically exercises the untested case: term_1 (v_1, rho_1, kappa_1)
% genuinely active for CV2/FWHM, with its OWN shape distinct from whatever
% arbitrary f_lo/f_hi Stage 2's dead (w_1~0) fit happened to freeze.

thetaTrueRate = struct('W_c',2.0,'u_half',84,'f_max_c',1850, ...
    'f_lo',7.9e-10,'f_hi',4.78e5,'u_M_half',61,'f_max_c2',15.75, ...    % deliberately the kind of ARBITRARY, unconstrained value Stage 2's flat objective could freeze given w_1~0 -- Stage 3 must NOT use this
    'w_1',2.22e-14,'w_2',2.22e-14);
thetaTrueShape = struct('v_1',6.0,'v_2',4.0, ...                        % BOTH genuinely nonzero now (previously only v_2 was tested)
    'rho_1',45,'kappa_1',1.1, ...                                       % term_1's OWN shape for CV2/FWHM, distinct from thetaTrueRate.f_lo/f_hi above
    'CV2_cap',0.25,'CV2_1',0.85,'CV2_2',0.45, ...
    'FWHM_cap',0.7,'FWHM_1',1.5,'FWHM_2',1.0);
f_lo_true = thetaTrueShape.rho_1 * exp(-thetaTrueShape.kappa_1);
f_hi_true = thetaTrueShape.rho_1 * exp(thetaTrueShape.kappa_1);

uMlevels = [10 50 100]; uElevels = [10 100 1000];
rawRows = repmat(empty_row_local(), 0, 1);
noiseSD_rate = 0.02; noiseSD_shape = 0.02;

pShapeBase = thetaTrueRate;
pShapeBase.v_1 = thetaTrueShape.v_1; pShapeBase.v_2 = thetaTrueShape.v_2;
pShapeBase.f_lo = f_lo_true; pShapeBase.f_hi = f_hi_true;   % term_1's OWN shape, not thetaTrueRate.f_lo/f_hi
pShapeBase.CV2_cap = thetaTrueShape.CV2_cap; pShapeBase.CV2_1 = thetaTrueShape.CV2_1; pShapeBase.CV2_2 = thetaTrueShape.CV2_2;
pShapeBase.FWHM_cap = thetaTrueShape.FWHM_cap; pShapeBase.FWHM_1 = thetaTrueShape.FWHM_1; pShapeBase.FWHM_2 = thetaTrueShape.FWHM_2;

for a = 1:4
    % E-alone rows (u_M=0, u_E>0) -- NEW: needed now that Stage 3 measures
    % CV2_cap/FWHM_cap from these rather than fitting them. At u_M=0,
    % cv2_hat/fwhm_hat equal CV2_cap/FWHM_cap EXACTLY (g1(0,.)=g2(0,.)=0), so
    % this synthetic generation is just adding noise around the same ground truth.
    for ei = 1:numel(uElevels)
        for trial = 1:3
            cv2Hat  = model_layer12_equations('cv2_hat',  0, uElevels(ei), pShapeBase);
            fwhmHat = model_layer12_equations('fwhm_hat', 0, uElevels(ei), pShapeBase);
            r = empty_row_local();
            r.animal = sprintf('a%d', a);
            r.condition = sprintf('E%d', uElevels(ei));
            r.phase = 'recovery'; r.label = 'RVN';
            r.mean.rate = model_layer12_equations('r_vagus', 0, uElevels(ei), thetaTrueRate) * (1 + noiseSD_rate*randn());
            r.mean.cv2  = cv2Hat  * (1 + noiseSD_shape*randn());
            r.mean.fwhm = fwhmHat * (1 + noiseSD_shape*randn());
            rawRows(end+1) = r; %#ok<AGROW>
        end
    end
    for mi = 1:numel(uMlevels)
        for ei = 1:numel(uElevels)
            for trial = 1:3
                rateHat = model_layer12_equations('r_vagus', uMlevels(mi), uElevels(ei), thetaTrueRate);
                cv2Hat  = model_layer12_equations('cv2_hat',  uMlevels(mi), uElevels(ei), pShapeBase);
                fwhmHat = model_layer12_equations('fwhm_hat', uMlevels(mi), uElevels(ei), pShapeBase);

                r = empty_row_local();
                r.animal = sprintf('a%d', a);
                r.condition = sprintf('M%dE%d', uMlevels(mi), uElevels(ei));
                r.phase = 'recovery'; r.label = 'RVN';
                r.mean.rate = rateHat * (1 + noiseSD_rate*randn());
                r.mean.cv2  = cv2Hat  * (1 + noiseSD_shape*randn());
                r.mean.fwhm = fwhmHat * (1 + noiseSD_shape*randn());
                rawRows(end+1) = r; %#ok<AGROW>
            end
        end
    end
end

fprintf('Synthetic rawRows built: %d rows. rate w_1=%.4g, w_2=%.4g (near-zero); rate''s OWN f_lo/f_hi are the ARBITRARY %.4g/%.4g (must NOT be used by Stage 3).\n', ...
    numel(rawRows), thetaTrueRate.w_1, thetaTrueRate.w_2, thetaTrueRate.f_lo, thetaTrueRate.f_hi);
fprintf('term_1''s TRUE, distinct shape for CV2/FWHM: f_lo=%.4g, f_hi=%.4g (rho_1=%.4g, kappa_1=%.4g)\n', ...
    f_lo_true, f_hi_true, thetaTrueShape.rho_1, thetaTrueShape.kappa_1);

thetaShape = fit_layer12_stage3_shape(rawRows, 'RVN', thetaTrueRate);

fprintf('\n================ GROUND TRUTH RECOVERY ================\n');
fprintf('  %-10s true=%.4f  measured=%.4f (from E-alone means, not fit)\n', 'CV2_cap', thetaTrueShape.CV2_cap, thetaShape.CV2_cap);
fprintf('  %-10s true=%.4f  measured=%.4f (from E-alone means, not fit)\n', 'FWHM_cap', thetaTrueShape.FWHM_cap, thetaShape.FWHM_cap);
checkFields = {'v_1','v_2','rho_1','kappa_1','CV2_1','CV2_2','FWHM_1','FWHM_2'};
allClose = true;
for f = 1:numel(checkFields)
    fn = checkFields{f};
    truth = thetaTrueShape.(fn); fit = thetaShape.(fn);
    pctErr = 100*abs(fit-truth)/truth;
    ok = pctErr < 20;
    allClose = allClose && ok;
    fprintf('  %-10s true=%.4f  fit=%.4f  pctErr=%.2f%%  %s\n', fn, truth, fit, pctErr, string(ok));
end
fprintf('  %-10s true=%.4f  fit=%.4f (derived from rho_1/kappa_1)\n', 'f_lo', f_lo_true, thetaShape.f_lo);
fprintf('  %-10s true=%.4f  fit=%.4f (derived from rho_1/kappa_1)\n', 'f_hi', f_hi_true, thetaShape.f_hi);
fprintf('  NOT using rate''s dead-fit f_lo/f_hi (%.4g/%.4g) -- confirms decoupling from Stage 2.\n', thetaTrueRate.f_lo, thetaTrueRate.f_hi);

fprintf('\n================ REQUIRED DIAGNOSTIC (a): restart-spread flatness on rho_1/kappa_1 ================\n');
ar = thetaShape.stage3_allRestarts;
resn = arrayfun(@(s) s.resnorm, ar);
best = min(resn); near = resn <= best*1.01;
thetas = vertcat(ar.theta);
rho1Near = thetas(near,3); kappa1Near = thetas(near,4);
v1Near = thetas(near,1); v2Near = thetas(near,2);
fprintf('  v_1 range among near-best restarts:     %.4g (true=%.2f)\n', max(v1Near)-min(v1Near), thetaTrueShape.v_1);
fprintf('  v_2 range among near-best restarts:     %.4g (true=%.2f)\n', max(v2Near)-min(v2Near), thetaTrueShape.v_2);
fprintf('  rho_1 range among near-best restarts:   %.4g (true=%.2f)\n', max(rho1Near)-min(rho1Near), thetaTrueShape.rho_1);
fprintf('  kappa_1 range among near-best restarts: %.4g (true=%.2f)\n', max(kappa1Near)-min(kappa1Near), thetaTrueShape.kappa_1);
flatnessOK = (max(v1Near)-min(v1Near) < 0.5*thetaTrueShape.v_1) && (max(v2Near)-min(v2Near) < 0.5*thetaTrueShape.v_2) ...
    && (max(rho1Near)-min(rho1Near) < 0.5*thetaTrueShape.rho_1) && (max(kappa1Near)-min(kappa1Near) < 0.5*thetaTrueShape.kappa_1);
fprintf('  Non-flat (well-constrained), as expected when term_1 is genuinely active: %s\n', string(flatnessOK));

fprintf('\n================ REQUIRED DIAGNOSTIC (b): resnorm with vs without term_1/term_2 ================\n');
fprintf('  CV2  improvement: %.2f%%\n', thetaShape.stage3_cv2ImprovementPct);
fprintf('  FWHM improvement: %.2f%%\n', thetaShape.stage3_fwhmImprovementPct);
improvementOK = thetaShape.stage3_cv2ImprovementPct > 20 && thetaShape.stage3_fwhmImprovementPct > 20;
fprintf('  Real resnorm improvement (not 0.00%%): %s\n', string(improvementOK));

fprintf('\n================ OVERALL ================\n');
if allClose && flatnessOK && improvementOK
    fprintf('PASS: ground truth recovered (including term_1''s OWN rho_1/kappa_1, distinct from rate''s dead-fit values), restarts well-constrained, real resnorm improvement.\n');
    fprintf('Safe to proceed to real data.\n');
else
    fprintf(2, 'FAIL: at least one check did not pass -- DO NOT proceed to real data.\n');
end

function r = empty_row_local()
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'rate_t',[],'vpp',[],'fwhm',[],'spk_t',[],'cv2',[],'cv2_t',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN, 'stem','');
end
