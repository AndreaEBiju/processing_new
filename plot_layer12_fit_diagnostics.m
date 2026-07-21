function figs = plot_layer12_fit_diagnostics(rawRows, channelLabel, thetaFull, R, saveDir, metricLabelOverride)
%PLOT_LAYER12_FIT_DIAGNOSTICS Diagnostics for the Layer 1-2 fit, one channel.
%
%   figs = plot_layer12_fit_diagnostics(rawRows, 'RVN', thetaFull)
%   figs = plot_layer12_fit_diagnostics(rawRows, 'RVN', thetaFull, R)          % R = bulk_mixed_models(...) output, adds synergy cross-check
%   figs = plot_layer12_fit_diagnostics(rawRows, 'RVN', thetaFull, R, saveDir) % also saves png/svg/fig
%   figs = plot_layer12_fit_diagnostics(rawRows, 'TOTAL', thetaFull, R, saveDir, 'Total firing rate (LVN+RVN)')
%       % metricLabelOverride: use when R.cells' metric name doesn't follow
%       % the default 'Nerve firing rate (<channelLabel>)' pattern -- e.g.
%       % run_layer12_total_channel.m's TOTAL pseudo-channel, whose matching
%       % bulk_mixed_models row is actually named 'Total firing rate (LVN+RVN)'.
%
% Produces, per matlab_implementation_instructions.md Section 6:
%   1. Observed vs predicted scatter -- Stage 1 (E-alone) / Stage 2 (interaction) panels
%   2. Predicted M x E heatmap (this model) side by side with bulk_mixed_models'
%      per-cell synergy heatmap for the same metric/channel (if R given)
%   3. Residuals vs u_M and vs u_E, separately
%   4. Section 5 pass/fail pattern check, printed as text overlay + console
%
% Reuses pubfig_setup / me_heatmap_render / save_one_figure conventions
% (matlab_implementation_instructions.md Section 0.8) -- do not restyle.

    if nargin < 4; R = []; end
    if nargin < 5; saveDir = ''; end
    if nargin < 6; metricLabelOverride = ''; end
    formats = {'png','svg','fig'};

    if exist('pubfig_setup','file')==2
        try; pubfig_setup('Theme','light','BaseFontSize',14,'LineWidth',2.0,'MarkerSize',10,'EnableLaTeX',false);
        catch ME; warning('plot_layer12:pubfig','pubfig_setup failed (%s)', ME.message); end
    end

    % ---- gather trials for this channel, recovery phase ---------------
    isRec = strcmpi({rawRows.phase}, 'recovery');
    isCh  = strcmpi({rawRows.label}, channelLabel);
    n = numel(rawRows);
    uM = nan(n,1); uE = nan(n,1);
    for i = 1:n; [uM(i), uE(i)] = parse_me_local(rawRows(i).condition); end
    rate = nan(n,1);
    for i = 1:n
        if isfield(rawRows(i),'mean') && isfield(rawRows(i).mean,'rate'); rate(i) = rawRows(i).mean.rate; end
    end
    sel = find(isRec(:) & isCh(:) & isfinite(rate) & isfinite(uM) & isfinite(uE) & (uE > 0));
    uMs = uM(sel); uEs = uE(sel); obs = rate(sel);
    pred = model_layer12_equations('r_vagus', uMs, uEs, thetaFull);
    resid = obs - pred;
    isEalone   = uMs == 0;
    isCombined = uMs > 0;

    check1 = NaN; check2 = NaN;
    if isfield(thetaFull,'check1_M50E100_gt_M100E100'); check1 = thetaFull.check1_M50E100_gt_M100E100; end
    if isfield(thetaFull,'check2_M100E10_gt_M50E10');   check2 = thetaFull.check2_M100E10_gt_M50E10;   end
    checkStr = sprintf('Check1 (M50xE100>M100xE100): %s   Check2 (M100xE10>M50xE10): %s', ...
        bool2str(check1), bool2str(check2));
    fprintf('[diagnostics-%s] %s\n', channelLabel, checkStr);

    figs = struct();

    % ==== Figure 1: observed vs predicted, Stage1 / Stage2 panels ======
    f1 = figure('Color','w','Name',sprintf('%s: observed vs predicted', channelLabel), ...
        'Position',[100 100 1000 480]);
    tl = tiledlayout(f1, 1, 2, 'Padding','compact');
    title(tl, sprintf('%s -- observed vs predicted rate', channelLabel), 'Interpreter','none');
    subtitle(tl, checkStr, 'Interpreter','none');

    nexttile(tl); hold on;
    scatter(pred(isEalone), obs(isEalone), 60, 'filled');
    add_identity_line();
    xlabel('predicted rate\_hat (Hz)'); ylabel('observed rate (Hz)');
    title('Stage 1: E-alone', 'Interpreter','none');

    nexttile(tl); hold on;
    scatter(pred(isCombined), obs(isCombined), 60, 'filled');
    add_identity_line();
    xlabel('predicted rate\_hat (Hz)'); ylabel('observed rate (Hz)');
    title('Stage 2: combined M x E', 'Interpreter','none');
    figs.obsPred = f1;

    % ==== Figure 2: predicted M x E heatmap (+ bulk_mixed_models cross-check) ====
    mLevAll = [10 50 100]; eLevAll = [10 100 1000];
    Zpred = nan(numel(mLevAll)+1, numel(eLevAll)+1);
    for j = 1:numel(eLevAll)
        Zpred(1,1+j) = model_layer12_equations('r_vagus', 0, eLevAll(j), thetaFull);
    end
    Zpred(2:end,1) = 0;   % M-alone: exactly 0 by construction (L1) -- see model_layer12_equations.m
    for i = 1:numel(mLevAll)
        for j = 1:numel(eLevAll)
            Zpred(1+i,1+j) = model_layer12_equations('r_vagus', mLevAll(i), eLevAll(j), thetaFull);
        end
    end
    fPredHm = me_heatmap_render(Zpred, mLevAll, eLevAll, ...
        sprintf('%s: predicted rate\\_hat (Hz)', channelLabel), ...
        'predicted rate\_hat (Hz)  (M-alone = 0 by construction, L1)', false(size(Zpred)), '%.2f');
    figs.predictedHeatmap = fPredHm;

    figs.synergyHeatmap = [];
    if ~isempty(R) && isfield(R,'cells')
        if ~isempty(metricLabelOverride)
            metricLabel = metricLabelOverride;
        else
            metricLabel = sprintf('Nerve firing rate (%s)', channelLabel);
        end
        rows = R.cells(strcmpi(R.cells.metric, metricLabel), :);
        if isempty(rows)
            warning('plot_layer12:noSynergyRows', ...
                '%s: no bulk_mixed_models cells found for metric "%s" -- skipping synergy cross-check.', ...
                channelLabel, metricLabel);
        else
            Zsyn = nan(numel(mLevAll)+1, numel(eLevAll)+1);
            Ssig = false(size(Zsyn));
            for r = 1:height(rows)
                i = find(mLevAll==rows.M(r)); j = find(eLevAll==rows.E(r));
                if isempty(i) || isempty(j); continue; end
                Zsyn(1+i,1+j) = 100*rows.synergy(r);
                Ssig(1+i,1+j) = rows.pFDR(r) < 0.05;
            end
            fSynHm = me_heatmap_render(Zsyn, mLevAll, eLevAll, ...
                sprintf('%s: bulk\\_mixed\\_models interaction (%% synergy)', channelLabel), ...
                'M x E interaction coefficient (%, * FDR<0.05) -- marginals not available from R.cells', ...
                Ssig);
            figs.synergyHeatmap = fSynHm;
        end
    else
        fprintf('[diagnostics-%s] no bulk_mixed_models output (R) supplied -- synergy cross-check skipped.\n', channelLabel);
    end

    % ==== Figure 3: residuals vs u_M, vs u_E ============================
    f3 = figure('Color','w','Name',sprintf('%s: residuals', channelLabel), 'Position',[100 100 1000 480]);
    tl3 = tiledlayout(f3, 1, 2, 'Padding','compact');
    title(tl3, sprintf('%s -- residuals (observed - predicted)', channelLabel), 'Interpreter','none');

    nexttile(tl3); hold on;
    scatter(uMs, resid, 60, 'filled');
    yline(0,'k--');
    xlabel('u\_M (Hz)'); ylabel('residual (Hz)');
    title('Residual vs u\_M', 'Interpreter','none');

    nexttile(tl3); hold on;
    scatter(uEs, resid, 60, 'filled');
    yline(0,'k--');
    xlabel('u\_E (Hz)'); ylabel('residual (Hz)');
    title('Residual vs u\_E', 'Interpreter','none');
    figs.residuals = f3;

    if ~isempty(saveDir) && exist('save_one_figure','file')==2
        save_one_figure(f1, saveDir, sprintf('layer12_%s_obsPred', channelLabel), formats);
        save_one_figure(fPredHm, saveDir, sprintf('layer12_%s_predictedHeatmap', channelLabel), formats);
        if ~isempty(figs.synergyHeatmap)
            save_one_figure(figs.synergyHeatmap, saveDir, sprintf('layer12_%s_synergyHeatmap', channelLabel), formats);
        end
        save_one_figure(f3, saveDir, sprintf('layer12_%s_residuals', channelLabel), formats);
    end
end

% ----------------------------------------------------------------------
function add_identity_line()
    ax = gca; xl = xlim(ax); yl = ylim(ax);
    lo = min([xl yl]); hi = max([xl yl]);
    plot(ax, [lo hi], [lo hi], 'k--');
    xlim(ax, xl); ylim(ax, yl);
end

function s = bool2str(tf)
    if isnan(tf); s = 'N/A'; elseif tf; s = 'PASS'; else; s = 'FAIL'; end
end

function [M,E] = parse_me_local(cond)
% Replicates parse_me in bulk_mixed_models.m exactly: M(\d+) / E(\d+),
% default 0 if the respective token is absent.
    c = upper(strtrim(char(cond))); M=0; E=0;
    t = regexp(c,'M(\d+)','tokens','once'); if ~isempty(t); M=str2double(t{1}); end
    t = regexp(c,'E(\d+)','tokens','once'); if ~isempty(t); E=str2double(t{1}); end
end
