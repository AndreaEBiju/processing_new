S = load('real_layer12_results.mat');
R = load('real_rawRows.mat');
rawRows = R.out.raw;

fprintf('================ Q1: full restart distribution for w_1, w_2 ================\n');
pc = S.resultsPerChannel.perChannelDetail;
channels = fieldnames(pc);
for k = 1:numel(channels)
    ch = channels{k};
    entry = pc.(ch);
    tf = entry.thetaFull;
    if ~isfield(tf,'allRestarts'); continue; end
    ar = tf.allRestarts;
    resn = arrayfun(@(s) s.resnorm, ar);
    best = min(resn);
    thetas = vertcat(ar.theta);
    w1 = thetas(:,5); w2 = thetas(:,6);
    nearMask = resn <= best*1.01;
    fprintf('\n-- %s (winning candidate) --\n', ch);
    fprintf('  best resnorm = %.6g\n', best);
    fprintf('  all 20 restarts: resnorm range [%.6g, %.6g]\n', min(resn), max(resn));
    fprintf('  %-3d of 20 restarts are near-best (<=1.01x best)\n', nnz(nearMask));
    fprintf('  Among near-best restarts:\n');
    fprintf('    w_1: %s\n', mat2str(sort(w1(nearMask))', 4));
    fprintf('    w_2: %s\n', mat2str(sort(w2(nearMask))', 4));
    fprintf('  Among ALL 20 restarts (resnorm, w_1, w_2), sorted by resnorm:\n');
    [~, ord] = sort(resn);
    for i = 1:numel(ord)
        j = ord(i);
        fprintf('    resnorm=%.6g  w_1=%.4g  w_2=%.4g  exitflag=%d\n', resn(j), w1(j), w2(j), ar(j).exitflag);
    end
end

fprintf('\n================ Q2: does term_1/term_2 improve resnorm vs term_c alone? ================\n');
for k = 1:numel(channels)
    ch = channels{k};
    entry = pc.(ch);
    tf = entry.thetaFull;
    if isempty(tf); continue; end

    % pull this channel's combined-condition (M>0 & E>0) trials directly from rawRows
    sel = arrayfun(@(r) strcmpi(r.label, ch), rawRows);
    rr = rawRows(sel);
    uM = nan(numel(rr),1); uE = nan(numel(rr),1); rate = nan(numel(rr),1);
    for i = 1:numel(rr)
        c = upper(strtrim(rr(i).condition));
        m = regexp(c,'M(\d+)','tokens','once'); if isempty(m); mm=0; else; mm=str2double(m{1}); end
        e = regexp(c,'E(\d+)','tokens','once'); if isempty(e); ee=0; else; ee=str2double(e{1}); end
        uM(i) = mm; uE(i) = ee; rate(i) = rr(i).mean.rate;
    end
    comb = uM>0 & uE>0 & ~isnan(rate);
    uMc = uM(comb); uEc = uE(comb); ratec = rate(comb);

    % full model prediction (fitted w_1, w_2)
    predFull = model_layer12_equations('r_vagus', uMc, uEc, tf);
    resnormFull = sum((ratec - predFull).^2);

    % term_c-only prediction: same theta but w_1=0, w_2=0
    tf0 = tf; tf0.w_1 = 0; tf0.w_2 = 0;
    predC = model_layer12_equations('r_vagus', uMc, uEc, tf0);
    resnormC = sum((ratec - predC).^2);

    pctImprovement = 100*(resnormC - resnormFull)/resnormC;
    fprintf('\n-- %s --\n', ch);
    fprintf('  n combined-condition trials = %d\n', numel(ratec));
    fprintf('  resnorm, term_c ONLY        = %.6g\n', resnormC);
    fprintf('  resnorm, full model (fitted)= %.6g\n', resnormFull);
    fprintf('  improvement from adding term_1+term_2 = %.4f%%\n', pctImprovement);
end
