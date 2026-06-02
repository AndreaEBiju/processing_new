function plot_spike_raster(D, tWin)
% PLOT_SPIKE_RASTER  Raster of firing events for both neural channels.
%
%   plot_spike_raster(D)            % full recording
%   plot_spike_raster(D, [t0 t1])  % only the time window [t0 t1] seconds
%
% One row per neural channel; each detected (accepted) spike is drawn as a
% vertical line at its time. Invalid/blanked regions are lightly shaded.
%
% Requires D.spikes(k).alignedTimes (step 4).

    if ~isfield(D, 'spikes') || ~isfield(D.spikes, 'alignedTimes')
        error('plot_spike_raster:missing', 'Run step4_waveforms first.');
    end
    fs = D.fs; N = size(D.filtered, 1); Tend = (N - 1) / fs;
    if nargin < 2 || isempty(tWin); tWin = [0 Tend]; end
    ch = D.neuralChannels; nCh = numel(ch); lab = D.channelLabels;
    cond = ''; if isfield(D, 'condition'); cond = [' | ' D.condition]; end

    figure('Color', 'w', 'Name', 'Spike raster', 'Position', [120 200 1250 360]);
    ax = axes; hold(ax, 'on');

    % light shading of invalid regions (use channel 1's mask as reference)
    if isfield(D, 'validMask')
        t = (0:N-1) / fs;
        shade_runs(ax, t, ~D.validMask(:, 1), [0.5 nCh + 0.5]);
    end

    for k = 1:nCh
        st = D.spikes(k).alignedTimes(:);
        st = st(st >= tWin(1) & st <= tWin(2));
        y0 = k - 0.4; y1 = k + 0.4;
        X = [st.'; st.'; nan(1, numel(st))];
        Y = [repmat(y0, 1, numel(st)); repmat(y1, 1, numel(st)); nan(1, numel(st))];
        plot(ax, X(:), Y(:), 'k', 'LineWidth', 0.2);
    end

    set(ax, 'YTick', 1:nCh, 'YTickLabel', lab, 'YDir', 'reverse');
    ylim(ax, [0.5 nCh + 0.5]); xlim(ax, tWin);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Channel');
    title(ax, sprintf('Spike raster — firing events%s', cond), 'Interpreter', 'none');
end

% ========================================================================
function shade_runs(ax, t, mask, yl)
    if ~any(mask); return; end
    d = diff([0; mask(:); 0]);
    s = find(d == 1); e = find(d == -1) - 1;
    for i = 1:min(numel(s), 500)
        x0 = t(s(i)); x1 = t(min(e(i), numel(t)));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.5 0.5], ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
