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
noiseSD_rate = 0.0001; noiseSD_shape = 0.0001;
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

thetaShape = fit_layer12_stage3_shape(rawRows, 'RVN', thetaTrueRate);

fprintf('\n================ Step 1: candidate invariant combinations ================\n');
v1f = thetaShape.v_1; v2f = thetaShape.v_2; kap1f = thetaShape.kappa_1;
CV2capf = thetaShape.CV2_cap; CV21f = thetaShape.CV2_1; CV22f = thetaShape.CV2_2;
FWHMcapf = thetaShape.FWHM_cap; FWHM1f = thetaShape.FWHM_1; FWHM2f = thetaShape.FWHM_2;

v1t = thetaTrueShape.v_1; v2t = thetaTrueShape.v_2; kap1t = thetaTrueShape.kappa_1;
CV2capt = thetaTrueShape.CV2_cap; CV21t = thetaTrueShape.CV2_1; CV22t = thetaTrueShape.CV2_2;
FWHMcapt = thetaTrueShape.FWHM_cap; FWHM1t = thetaTrueShape.FWHM_1; FWHM2t = thetaTrueShape.FWHM_2;

cand1 = v1f*(CV21f-CV2capf); cand1_true = v1t*(CV21t-CV2capt);
cand2 = v1f*(FWHM1f-FWHMcapf); cand2_true = v1t*(FWHM1t-FWHMcapt);
cand3 = v2f*(CV22f-CV2capf); cand3_true = v2t*(CV22t-CV2capt);
cand4 = v2f*(FWHM2f-FWHMcapf); cand4_true = v2t*(FWHM2t-FWHMcapt);

names = {'v1*(CV2_1-CV2_cap)','v1*(FWHM_1-FWHM_cap)','v2*(CV2_2-CV2_cap)','v2*(FWHM_2-FWHM_cap)'};
cands = [cand1 cand2 cand3 cand4]; candsTrue = [cand1_true cand2_true cand3_true cand4_true];
for i = 1:4
    pctErr = 100*abs(cands(i)-candsTrue(i))/abs(candsTrue(i));
    fprintf('  %-24s fit=%.6g  true=%.6g  pctErr=%.2f%%\n', names{i}, cands(i), candsTrue(i), pctErr);
end

fprintf('\n================ Step 2: numerical Jacobian + SVD on [v_1,v_2,kappa_1,CV2_cap,FWHM_cap] ================\n');
% fixed (well-recovered) params held at their FITTED values
rho_1 = thetaShape.rho_1; CV2_1 = thetaShape.CV2_1; CV2_2 = thetaShape.CV2_2;
FWHM_1 = thetaShape.FWHM_1; FWHM_2 = thetaShape.FWHM_2;
baseP = thetaTrueRate;   % W_c, u_half, f_max_c, u_M_half, f_max_c2 (fixed, unchanged)

[uMg, uEg] = meshgrid(uMlevels, uElevels);
uMv = uMg(:); uEv = uEg(:);

params0 = [thetaShape.v_1, thetaShape.v_2, thetaShape.kappa_1, thetaShape.CV2_cap, thetaShape.FWHM_cap];
J = zeros(18, 5);
epsStep = 1e-4;
for k = 1:5
    pPlus = params0; pPlus(k) = pPlus(k) * (1+epsStep);
    pMinus = params0; pMinus(k) = pMinus(k) * (1-epsStep);
    predsPlus  = compute_all_preds(pPlus,  baseP, rho_1, CV2_1, CV2_2, FWHM_1, FWHM_2, uMv, uEv);
    predsMinus = compute_all_preds(pMinus, baseP, rho_1, CV2_1, CV2_2, FWHM_1, FWHM_2, uMv, uEv);
    J(:,k) = (predsPlus - predsMinus) ./ (2*epsStep);   % d(pred)/d(log-param) since step is relative
end

[~,Sv,V] = svd(J);
singvals = diag(Sv);
fprintf('Singular values: %s\n', mat2str(singvals', 4));
fprintf('Relative to largest: %s\n', mat2str((singvals/max(singvals))', 4));
paramNames = {'v_1','v_2','kappa_1','CV2_cap','FWHM_cap'};
for k = 1:5
    fprintf('  V(:,%d) [%s]: %s  (singval=%.4g, rel=%.4g)\n', k, strjoin(paramNames,','), mat2str(V(:,k)',4), singvals(min(k,numel(singvals))), singvals(min(k,numel(singvals)))/max(singvals));
end

function y = compute_all_preds(vars, baseP, rho_1, CV2_1, CV2_2, FWHM_1, FWHM_2, uMv, uEv)
    p = baseP;
    p.v_1 = vars(1); p.v_2 = vars(2);
    kappa_1 = vars(3);
    p.f_lo = rho_1*exp(-kappa_1); p.f_hi = rho_1*exp(kappa_1);
    p.CV2_cap = vars(4); p.CV2_1 = CV2_1; p.CV2_2 = CV2_2;
    p.FWHM_cap = vars(5); p.FWHM_1 = FWHM_1; p.FWHM_2 = FWHM_2;
    cv2v = model_layer12_equations('cv2_hat', uMv, uEv, p);
    fwhmv = model_layer12_equations('fwhm_hat', uMv, uEv, p);
    y = [cv2v(:); fwhmv(:)];
end

function r = empty_row_local()
    r = struct('animal','','condition','','phase','','label','', ...
        'dist',struct('rate',[],'rate_t',[],'vpp',[],'fwhm',[],'spk_t',[],'cv2',[],'cv2_t',[]), ...
        'mean',struct('rate',NaN,'excess',NaN,'vpp',NaN,'fwhm',NaN,'cv2',NaN), ...
        'fanoSlope',NaN, 'stem','');
end
