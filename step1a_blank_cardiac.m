function D = step1a_blank_cardiac(D, P, plotMode)
% STEP1A_BLANK_CARDIAC  NaN-blank a window around every R-peak in the RAW signal.
%
%   D = step1a_blank_cardiac(D, P, plotMode)
%
% Run BEFORE step1_bandpass. Removing the sharp QRS from the raw signal
% *before* filtering prevents filter ringing (a sharp R-peak rings through a
% bandpass), and the blanked windows become invalid everywhere downstream:
%   * step1_bandpass interpolates across the NaNs to filter, then restores the
%     NaNs -- so the filter never sees the QRS and does not ring;
%   * step2's validity mask excludes the NaN samples;
%   * detection, rate, ISI and envelope all respect the validity mask.
%
% This replaces the old cardiac handling (template subtraction in step 1b and
% post-hoc spike removal in step 9), which are no longer needed. The window is
% P.cardiacRemoveWinMs (uniform, every channel).
%
% Requires D.y (raw, step 0) and D.rpeakSamples.
% Adds D.cardiacBlank (logical, blanked samples) and D.cardiacBlankWinMs.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'rpeakSamples') || isempty(D.rpeakSamples)
        warning('step1a_blank_cardiac:noR', 'No R-peaks available; nothing blanked.');
        D.cardiacBlank = false(size(D.y,1),1); D.cardiacBlankWinMs = 0; return;
    end

    fs = D.fs; N = size(D.y, 1);
    w  = round(P.cardiacRemoveWinMs * 1e-3 * fs);
    rs = round(D.rpeakSamples(:)); rs = rs(rs >= 1 & rs <= N);

    blank = false(N, 1);
    for j = 1:numel(rs)
        a = max(1, rs(j)-w); b = min(N, rs(j)+w); blank(a:b) = true;
    end

    % capture a representative window BEFORE blanking, for the diagnostic
    pd = [];
    if plotMode; pd = capture_window(D, blank, fs); end

    D.y(blank, :) = NaN;
    D.cardiacBlank = blank;
    D.cardiacBlankWinMs = P.cardiacRemoveWinMs;

    fprintf('[step1a] Cardiac blank: +/-%g ms around %d R-peaks -> %.1f%% of samples NaN.\n', ...
        P.cardiacRemoveWinMs, numel(rs), 100*nnz(blank)/N);

    if plotMode && ~isempty(pd); plot_blank(pd); end
end

% ======================================================================
function pd = capture_window(D, blank, fs)
    ch = D.neuralChannels(1); N = size(D.y,1);
    if isfield(D,'rpeakTimes') && ~isempty(D.rpeakTimes)
        c = median(D.rpeakTimes);
    else
        c = (N-1)/fs/2;
    end
    lo = max(0, c-2.5); hi = min((N-1)/fs, lo+5);
    i0 = max(1, floor(lo*fs)+1); i1 = min(N, floor(hi*fs)+1);
    pd.t = (i0-1:i1-1)/fs;
    pd.raw = D.y(i0:i1, ch) * 1e6;     % original (pre-blank), uV
    pd.blank = blank(i0:i1);
    pd.label = D.channelLabels{1};
end

function plot_blank(pd)
    figure('Color','w','Name','Step 1a — cardiac blank','Position',[200 200 1050 360]);
    plot(pd.t, pd.raw, 'Color', [0.3 0.3 0.3]); hold on;
    yl = ylim;
    d = diff([0; pd.blank(:); 0]); s = find(d==1); e = find(d==-1)-1;
    for i = 1:numel(s)
        x0 = pd.t(s(i)); x1 = pd.t(min(e(i), numel(pd.t)));
        patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [1 0.5 0.5], ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    ylim(yl); grid on; xlabel('Time (s)'); ylabel('Raw (\muV)');
    title(sprintf('%s: raw signal with QRS windows blanked (shaded) before filtering', pd.label), ...
        'Interpreter', 'none');
end
