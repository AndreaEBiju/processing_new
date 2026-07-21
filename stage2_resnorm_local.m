function rn = stage2_resnorm_local(rows, channelLabel, thetaFull)
%STAGE2_RESNORM_LOCAL Sum of squared residuals for one channel's Stage 2 fit,
%   computed independently of fit_layer12_stage2_interaction's own return
%   values (self-contained, no assumption about its output signature).
    rn = 0;
    for i = 1:numel(rows)
        if ~strcmpi(rows(i).label, channelLabel); continue; end
        if ~strcmpi(rows(i).phase, 'recovery'); continue; end
        tok = regexp(rows(i).condition, '^M(\d+)E(\d+)$', 'tokens', 'once');
        if isempty(tok); continue; end   % skips M-alone / E-alone rows
        uM = str2double(tok{1}); uE = str2double(tok{2});
        predicted = model_layer12_equations('r_vagus', uM, uE, thetaFull);
        rn = rn + (rows(i).mean.rate - predicted)^2;
    end
end