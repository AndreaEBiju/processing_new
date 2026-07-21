function D = step7_dfa_report(D, P, plotMode)
% STEP7_DFA_REPORT  Gap-aware DFA (fractal scaling) on the Step 6
% blanking-corrected firing-rate trace, per neural channel.
%
%   D = step7_dfa_report(D, P, plotMode)
%
% Requires D.metrics(k).fr_t / .fr_hz / .fr_validFrac (Step 6). Heartbeat
% blanking is already dead-time-corrected inside firing_rate()'s
% denominator and is invisible here; motion-artifact blanking is the only
% DFA-relevant gap source for this channel, and shows up as bins whose
% fr_validFrac falls below the floor already used for the Step 3b
% envelope (vfrac >= 0.5, see step3b_envelope.m).

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'metrics') || ~isfield(D.metrics, 'fr_hz')
        error('step7_dfa_report:missing', 'Run step6_spike_report first.');
    end

    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    validFracFloor = 0.5; % reuse step3b_envelope's existing floor (vfrac >= 0.5), not a new threshold

    D.dfa = repmat(emptyDfa(P), 1, nCh);

    for k = 1:nCh
        M = D.metrics(k);
        binValidMask = M.fr_validFrac(:) >= validFracFloor;

        try
            D.dfa(k) = dfaGapAware(M.fr_hz(:), binValidMask, P);
            fprintf(['[step7] ch %d (%s): alpha1=%.3f (R2=%.2f) | alpha2=%.3f (R2=%.2f) | ' ...
                     'alphaFull=%.3f (R2=%.2f) | %d/%d scales excluded\n'], ...
                ch(k), lab{k}, D.dfa(k).alpha1, D.dfa(k).R2_1, D.dfa(k).alpha2, D.dfa(k).R2_2, ...
                D.dfa(k).alphaFull, D.dfa(k).R2_full, numel(D.dfa(k).excludedScales), numel(D.dfa(k).nVals) + numel(D.dfa(k).excludedScales));
        catch ME
            if strcmp(ME.identifier, 'dfaGapAware:noValidRuns')
                warning('step7:noValidRuns', 'ch %d (%s): %s', ch(k), lab{k}, ME.message);
                D.dfa(k) = emptyDfa(P);
            else
                rethrow(ME);
            end
        end

        if plotMode
            plot_dfa(M.fr_t(:), M.fr_hz(:), binValidMask, D.dfa(k), lab{k});
        end
    end
end

% ========================================================================
function d = emptyDfa(P)
    d = struct('nVals', [], 'F', [], 'nWindows', [], 'alpha1', NaN, 'alpha2', NaN, ...
        'alphaFull', NaN, 'R2_1', NaN, 'R2_2', NaN, 'R2_full', NaN, 'nCross', NaN, ...
        'localAlpha', struct('n', [], 'slope', []), 'runs', [], 'excludedScales', [], ...
        'mu', NaN, 'P', P);
end

% ========================================================================
function plot_dfa(t, fr_hz, validMask, dfaOut, label)
    fig = figure('Color', 'w', 'Name', sprintf('Step 7 — %s', label), ...
        'Position', [160 100 1000 780]);
    tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 7 gap-aware DFA  |  %s', label), 'Interpreter', 'none');

    % ---- Panel 1: firing-rate trace with gap shading + run boundaries ----
    ax1 = nexttile(tl);
    plot(ax1, t, fr_hz, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold(ax1, 'on');
    shade_runs(ax1, t, ~validMask);
    if ~isempty(dfaOut.runs)
        for r = 1:size(dfaOut.runs, 1)
            xline(ax1, t(min(dfaOut.runs(r,1), numel(t))), 'g:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            xline(ax1, t(min(dfaOut.runs(r,2), numel(t))), 'g:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    title(ax1, sprintf('Firing-rate bins (%d bins, %.0f%% valid, %d run(s))', ...
        numel(t), 100*mean(validMask), size(dfaOut.runs,1)), 'Interpreter', 'none');
    ylabel(ax1, 'fr (Hz)'); xlabel(ax1, 'Time (s)');

    % ---- Panel 2: log-log F(n) with excluded scales + crossover ----
    ax2 = nexttile(tl);
    if ~isempty(dfaOut.nVals)
        loglog(ax2, dfaOut.nVals, dfaOut.F, 'ko-', 'MarkerFaceColor', 'k'); hold(ax2, 'on');
        if isfinite(dfaOut.alpha1)
            idxS = dfaOut.nVals >= dfaOut.P.dfaShortRange(1) & dfaOut.nVals <= dfaOut.P.dfaShortRange(2);
            plot(ax2, dfaOut.nVals(idxS), dfaOut.F(idxS), 'r-', 'LineWidth', 2);
        end
        if isfinite(dfaOut.alpha2)
            idxL = dfaOut.nVals > dfaOut.P.dfaLongRangeMinScale;
            plot(ax2, dfaOut.nVals(idxL), dfaOut.F(idxL), 'b-', 'LineWidth', 2);
        end
        if isfinite(dfaOut.nCross)
            xline(ax2, dfaOut.nCross, '--', sprintf('crossover n=%.1f', dfaOut.nCross));
        end
    end
    if ~isempty(dfaOut.excludedScales)
        yl2 = ylim(ax2);
        for es = dfaOut.excludedScales(:)'
            plot(ax2, es, yl2(1), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    grid(ax2, 'on'); xlabel(ax2, 'n (log scale)'); ylabel(ax2, 'F(n) (log scale)');
    title(ax2, sprintf('alpha1=%.3f (R2=%.2f)  alpha2=%.3f (R2=%.2f)  alphaFull=%.3f (R2=%.2f)', ...
        dfaOut.alpha1, dfaOut.R2_1, dfaOut.alpha2, dfaOut.R2_2, dfaOut.alphaFull, dfaOut.R2_full), ...
        'Interpreter', 'none', 'FontSize', 9);

    % ---- Panel 3: pooled window count per scale (QC floor) + local slope ----
    ax3 = nexttile(tl);
    if ~isempty(dfaOut.nVals)
        yyaxis(ax3, 'left');
        bar(ax3, dfaOut.nVals, dfaOut.nWindows, 'FaceColor', [0.7 0.7 0.9]); hold(ax3, 'on');
        yline(ax3, dfaOut.P.dfaMinWindowsPerScale, 'r--');
        ylabel(ax3, 'pooled windows (QC)');
        yyaxis(ax3, 'right');
        plot(ax3, dfaOut.localAlpha.n, dfaOut.localAlpha.slope, 'g.-');
        ylabel(ax3, 'local slope alpha(n)');
    end
    xlabel(ax3, 'n');
    title(ax3, 'Window-count floor (bars, red line = P.dfaMinWindowsPerScale) and local exponent (line)', ...
        'Interpreter', 'none', 'FontSize', 8);
end

% ========================================================================
function shade_runs(ax, t, mask)
% Shade contiguous true-runs of mask as translucent patches (matches
% step2_noise_sigma.m's shade_runs idiom).
    if ~any(mask); return; end
    yl = ylim(ax);
    d = diff([0; mask(:); 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    maxRuns = 300;
    for i = 1:min(numel(starts), maxRuns)
        x0 = t(starts(i)); x1 = t(min(ends(i), numel(t)));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.4 0.4], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end
