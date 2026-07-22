S = load('real_layer12_results.mat');

fprintf('================ Step 1: f_lo/f_hi (rho1/kappa1) restart stability given w_1~0 ================\n');

% RVN, LVN from the per-channel driver
pc = S.resultsPerChannel.perChannelDetail;
channels = {'RVN','LVN'};
for k = 1:numel(channels)
    ch = channels{k};
    tf = pc.(ch).thetaFull;
    report_stability(ch, tf);
end

% TOTAL from the total-channel driver
tfTotal = S.resultsTotal.thetaFull;
report_stability('TOTAL', tfTotal);

function report_stability(ch, tf)
    if ~isfield(tf,'allRestarts')
        fprintf('\n-- %s: no allRestarts field captured, cannot run this diagnostic --\n', ch);
        return;
    end
    ar = tf.allRestarts;
    thetas = vertcat(ar.theta);   % columns: [rho1, kappa1, u_M_half, f_max_c2, w_1, w_2]
    resn = arrayfun(@(s) s.resnorm, ar);
    rho1 = thetas(:,1); kappa1 = thetas(:,2);
    w1 = thetas(:,5); w2 = thetas(:,6);
    f_lo = rho1 .* exp(-kappa1); f_hi = rho1 .* exp(kappa1);

    best = min(resn);
    near = resn <= best * 1.01;

    fprintf('\n-- %s (fitted w_1=%.4g, w_2=%.4g -- near-zero, per the established rate null) --\n', ch, tf.w_1, tf.w_2);
    fprintf('  best resnorm = %.6g, resnorm range across ALL %d restarts: [%.6g, %.6g]\n', ...
        best, numel(resn), min(resn), max(resn));
    fprintf('  %d of %d restarts are near-best (<=1.01x best)\n', nnz(near), numel(resn));
    fprintf('  Among near-best restarts:\n');
    fprintf('    rho1   range: [%.4g, %.4g]\n', min(rho1(near)), max(rho1(near)));
    fprintf('    kappa1 range: [%.4g, %.4g]\n', min(kappa1(near)), max(kappa1(near)));
    fprintf('    f_lo   range: [%.4g, %.4g]  (derived)\n', min(f_lo(near)), max(f_lo(near)));
    fprintf('    f_hi   range: [%.4g, %.4g]  (derived)\n', min(f_hi(near)), max(f_hi(near)));
    fprintf('    w_1    range: [%.4g, %.4g]\n', min(w1(near)), max(w1(near)));
    fprintf('    w_2    range: [%.4g, %.4g]\n', min(w2(near)), max(w2(near)));

    rho1RelSpread = (max(rho1(near))-min(rho1(near))) / max(median(rho1(near)), eps);
    resnormRelSpread = (max(resn(near))-min(resn(near))) / max(best, eps);
    fprintf('  rho1 relative spread = %.2f  |  resnorm relative spread among near-best = %.2e\n', rho1RelSpread, resnormRelSpread);
    if rho1RelSpread > 0.1 && resnormRelSpread < 1e-3
        fprintf(2, '  VERDICT: f_lo/f_hi (rho1/kappa1) swing widely while resnorm stays flat -- SAME L17 signature as u_M_half/w_2. UNSTABLE, not data-informed. Frozen Stage 3 values are potentially arbitrary.\n');
    else
        fprintf('  VERDICT: f_lo/f_hi (rho1/kappa1) were stable despite w_1~0 -- data DID constrain them via term_c/f1 sharing u_half/f_max_c even with w_1 tiny. No fix needed; term_1-contributes-nothing conclusion stands.\n');
    end
end
