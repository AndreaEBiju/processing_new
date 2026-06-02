function D = step1_bandpass(D, P, plotMode)
% STEP1_BANDPASS  Zero-phase bandpass filter of the neural channels.
%
%   D = step1_bandpass(D, P, plotMode)
%
% Applies a 4th-order Butterworth bandpass (P.bandpassLow..P.bandpassHigh,
% high corner clamped to fs/2-1) with filtfilt (zero phase) to each neural
% channel in D.y(:, D.neuralChannels). NaN-blanked samples are filled for
% filtering then restored to NaN so blanks do not leak spectral energy.
%
% Adds to D:
%   D.filtered      samples x nNeural   bandpassed signal (NaN at blanks)
%   D.bandInfo      struct with the corners actually used
%
% When plotMode is true (single-file mode), for EVERY neural channel it shows
% the raw FFT magnitude before/after and the full-length time trace
% before/after. Set plotMode=false in bulk mode to suppress all plotting.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end

    fs   = D.fs;
    nyq  = fs / 2;
    lo   = P.bandpassLow;
    hi   = min(P.bandpassHigh, nyq - 1);
    if lo <= 0 || lo >= hi
        error('step1_bandpass:band', 'Bad band: low=%g hi=%g (fs=%g).', lo, hi, fs);
    end

    [b, a] = butter(P.filterOrder, [lo hi] / nyq, 'bandpass');

    ch     = D.neuralChannels;
    labels = D.channelLabels;
    nCh    = numel(ch);
    N      = size(D.y, 1);

    D.filtered = nan(N, nCh);
    D.bandInfo = struct('low', lo, 'high', hi, 'order', P.filterOrder, ...
                        'nyquist', nyq, 'requestedHigh', P.bandpassHigh, ...
                        'clamped', P.bandpassHigh > nyq - 1);

    fprintf('[step1] Bandpass %.1f-%.1f Hz (order %d, zero-phase). fs=%g, Nyq=%g.\n', ...
        lo, hi, P.filterOrder, fs, nyq);
    if D.bandInfo.clamped
        fprintf('        NOTE: requested high corner %.1f Hz clamped to %.1f Hz.\n', ...
            P.bandpassHigh, hi);
    end

    for k = 1:nCh
        x = D.y(:, ch(k));
        invalid = isnan(x);
        if any(invalid)
            xfill = fillmissing(x, 'linear', 'EndValues', 'nearest');
        else
            xfill = x;
        end
        xf = filtfilt(b, a, xfill);
        xf(invalid) = NaN;
        D.filtered(:, k) = xf;
        fprintf('        ch %d (%s): filtered (%d blanked samples)\n', ...
            ch(k), labels{k}, nnz(invalid));
    end

    if plotMode
        for k = 1:nCh
            plot_before_after(D.t, D.y(:, ch(k)), D.filtered(:, k), ...
                fs, labels{k}, lo, hi, P.plotMaxPoints);
        end
    end
end

% ========================================================================
function plot_before_after(t, xraw, xfilt, fs, label, lo, hi, maxPts)
% One figure per channel: row 1 = FFT magnitude before/after (overlaid),
% row 2 = full-length time trace before/after (overlaid, display-decimated).

    fig = figure('Color', 'w', 'Name', sprintf('Step 1 — %s', label), ...
        'Position', [120 120 1100 720]);
    tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 1 bandpass %.0f-%.0f Hz  |  %s', lo, hi, label), ...
        'Interpreter', 'none');

    % ---- FFT before / after ----
    [fR, mR] = single_sided_fft(xraw,  fs);
    [fF, mF] = single_sided_fft(xfilt, fs);

    nexttile;
    loglog(fR(2:end), mR(2:end), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.6); hold on;
    loglog(fF(2:end), mF(2:end), 'b', 'LineWidth', 0.8);
    xline(lo, 'r--'); xline(hi, 'r--');
    grid on; xlabel('Frequency (Hz)'); ylabel('|X(f)| (a.u.)');
    legend({'before', 'after', 'band'}, 'Location', 'southwest');
    title('FFT magnitude — before vs after');
    xlim([max(1, fR(2)) fs/2]);

    % ---- time trace before / after ----
    nexttile;
    n = numel(t);
    stride = max(1, floor(n / maxPts));
    idx = 1:stride:n;
    rawU = xraw(idx)  * 1e6;     % to uV (assumes volts; harmless scaling otherwise)
    filU = xfilt(idx) * 1e6;
    plot(t(idx), rawU, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4); hold on;
    plot(t(idx), filU, 'b', 'LineWidth', 0.4);
    grid on; xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    legend({'before', 'after'}, 'Location', 'northeast');
    if stride > 1
        title(sprintf('Full time trace — before vs after (display every %d^{th} sample)', stride));
    else
        title('Full time trace — before vs after');
    end
    xlim([t(1) t(end)]);
end

% ========================================================================
function [f, mag] = single_sided_fft(x, fs)
% Single-sided amplitude spectrum of the full signal. NaNs (blanks) are
% zero-filled so they do not propagate; this is for display only.
    x = x(:);
    x(~isfinite(x)) = 0;
    n = numel(x);
    X = fft(x);
    half = floor(n/2) + 1;
    mag = abs(X(1:half)) / n;
    if half > 2
        mag(2:end-1) = 2 * mag(2:end-1);
    end
    f = (0:half-1)' * (fs / n);
end
