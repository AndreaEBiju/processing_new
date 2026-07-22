rng_note = 'reuse the exact same synthetic generation as test_layer12_stage3_synthetic_smoke.m';

thetaTrueRate = struct('W_c',2.0,'u_half',84,'f_max_c',1850, ...
    'f_lo',7.9e-10,'f_hi',4.78e5,'u_M_half',61,'f_max_c2',15.75, ...
    'w_1',2.22e-14,'w_2',2.22e-14);
thetaTrueShape = struct('v_1',6.0,'v_2',4.0, ...
    'rho_1',45,'kappa_1',1.1, ...
    'CV2_cap',0.25,'CV2_1',0.85,'CV2_2',0.45, ...
    'FWHM_cap',0.7,'FWHM_1',1.5,'FWHM_2',1.0);
f_lo_true = thetaTrueShape.rho_1 * exp(-thetaTrueShape.kappa_1);
f_hi_true = thetaTrueShape.rho_1 * exp(thetaTrueShape.kappa_1);

uMlevels = [10 50 100]; uElevels = [10 100 1000];
rawRows = repmat(empty_row_local(), 0, 1);
noiseSD_rate = 0.02; noiseSD_shape = 0.02;
rngSeeded = false; %#ok<NASGU>
for a = 1:4
    for mi = 1:numel(uMlevels)
        for ei = 1:numel(uElevels)
            for trial = 1:3
                rateHat = model_layer12_equations('r_vagus', uMlevels(mi), uElevels(ei), thetaTrueRate);
                pShape = thetaTrueRate;
                pShape.v_1 = thetaTrueShape.v_1; pShape.v_2 = thetaTrueShape.v_2;
                pShape.f_lo = f_lo_true; pShape.f_hi = f_hi_true;
                pShape.CV2_cap = thetaTrueShape.CV2_cap; pShape.CV2_1 = thetaTrueShape.CV2_1; pShape.CV2_2 = thetaTrueShape.CV2_2;
                pShape.FWHM_cap = thetaTrueShape.FWHM_cap; pShape.FWHM_1 = thetaTrueShape.FWHM_1; pShape.FWHM_2 = thetaTrueShape.FWHM_2;
                cv2Hat  = model_layer12_equations('cv2_hat',  uMlevels(mi), uElevels(ei), pShape);
                fwhmHat = model_layer12_equations('fwhm_hat', uMlevels(mi), uElevels(ei), pShape);
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

% Recompute wCV2/wFWHM exactly as fit_layer12_stage3_shape.m does
isCombined = true(size(rawRows));
cv2v = arrayfun(@(r) r.mean.cv2, rawRows);
fwhmv = arrayfun(@(r) r.mean.fwhm, rawRows);
uMv = arrayfun(@(r) sscanf(regexp(r.condition,'M(\d+)','match','once'),'M%d'), rawRows);
uEv = arrayfun(@(r) sscanf(regexp(r.condition,'E(\d+)','match','once'),'E%d'), rawRows);
wCV2 = 1/var(cv2v); wFWHM = 1/var(fwhmv);

% resnorm at TRUE params
pTrue = thetaTrueRate;
pTrue.v_1 = thetaTrueShape.v_1; pTrue.v_2 = thetaTrueShape.v_2;
pTrue.f_lo = f_lo_true; pTrue.f_hi = f_hi_true;
pTrue.CV2_cap = thetaTrueShape.CV2_cap; pTrue.CV2_1 = thetaTrueShape.CV2_1; pTrue.CV2_2 = thetaTrueShape.CV2_2;
pTrue.FWHM_cap = thetaTrueShape.FWHM_cap; pTrue.FWHM_1 = thetaTrueShape.FWHM_1; pTrue.FWHM_2 = thetaTrueShape.FWHM_2;
cv2HatTrue = model_layer12_equations('cv2_hat', uMv, uEv, pTrue);
fwhmHatTrue = model_layer12_equations('fwhm_hat', uMv, uEv, pTrue);
resnormAtTruth = wCV2*sum((cv2v(:)-cv2HatTrue(:)).^2) + wFWHM*sum((fwhmv(:)-fwhmHatTrue(:)).^2);
fprintf('resnorm (CV2+FWHM only, weighted) AT TRUE params = %.6g\n', resnormAtTruth);

thetaShape = fit_layer12_stage3_shape(rawRows, 'RVN', thetaTrueRate);
fprintf('resnorm (CV2+FWHM only, weighted) AT FITTED best = %.6g\n', thetaShape.stage3_resnorm - sum((wCV2*0)));
% stage3_resnorm includes the fixed rate offset; strip it for fair comparison
rateHatAll = model_layer12_equations('r_vagus', uMv, uEv, thetaTrueRate);
rateVals = arrayfun(@(r) r.mean.rate, rawRows);
wRate = 1/var(rateVals);
rateOffset = wRate * sum((rateVals(:) - rateHatAll(:)).^2);
fprintf('resnorm (CV2+FWHM only) at fitted best (stripped of rate offset) = %.6g\n', thetaShape.stage3_resnorm - rateOffset);
fprintf('Difference: fitted - truth = %.6g\n', (thetaShape.stage3_resnorm - rateOffset) - resnormAtTruth);

function r = empty_row_local()
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'rate_t',[],'vpp',[],'fwhm',[],'spk_t',[],'cv2',[],'cv2_t',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN, 'stem','');
end
