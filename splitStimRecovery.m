function [x_stim, x_recovery, stimMask, recoveryMask, dig_aligned, detectInfo] = ...
    splitStimRecovery(vib, fs_vib, x, fs_sig, baseFileName, doPlot, threshold)
% splitStimRecovery
% Split signal matrix into stimulation and recovery parts using a noisy
% oscillatory control signal recorded at a different sampling rate.
%
% This version uses envelope-based detection:
%   1) smooth vib lightly
%   2) compute moving-RMS envelope
%   3) threshold envelope
%   4) remove short flickers
%   5) keep only the largest ON segment
%   6) align mask to the signal time base
%
% INPUTS
%   vib          : [Nv x 1] or [1 x Nv] noisy control/on-off signal
%   fs_vib       : sampling rate of vib
%   x            : [N x nChan] signal matrix to split
%   fs_sig       : sampling rate of x
%   baseFileName : base filename for saving, without suffix
%   doPlot       : true/false, whether to plot before/after (optional)
%   threshold    : envelope threshold (optional)
%
% OUTPUTS
%   x_stim       : samples where aligned digital signal is ON
%   x_recovery   : samples where aligned digital signal is OFF
%   stimMask     : logical mask on signal time base
%   recoveryMask : logical mask on signal time base
%   dig_aligned  : aligned logical stimulation vector on signal time base
%   detectInfo   : struct with debug information
%
% SAVED FILES
%   [baseFileName '_stim.mat']
%   [baseFileName '_recovery.mat']

    if nargin < 6 || isempty(doPlot)
        doPlot = true;
    end
    if nargin < 7
        threshold = [];
    end

    validateattributes(vib, {'numeric','logical'}, {'nonempty','vector'});
    validateattributes(fs_vib, {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(x, {'numeric'}, {'2d','nonempty','real','finite'});
    validateattributes(fs_sig, {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(baseFileName, {'char','string'}, {'nonempty'});

    vib = vib(:);
    baseFileName = char(baseFileName);

    [N, ~] = size(x);

    [dig_aligned, detectInfo] = alignDigitalToSignalEnvelope( ...
        vib, fs_vib, N, fs_sig, threshold);

    stimMask = dig_aligned(:) ~= 0;
    recoveryMask = ~stimMask;

    x_stim = x(stimMask, :);
    x_recovery = x(recoveryMask, :);

    stimFile = [baseFileName '_stim.mat'];
    recoveryFile = [baseFileName '_recovery.mat'];

    save(stimFile, 'x_stim', 'stimMask', 'dig_aligned', 'fs_sig', 'detectInfo');
    save(recoveryFile, 'x_recovery', 'recoveryMask', 'dig_aligned', 'fs_sig', 'detectInfo');

    if doPlot
        figure('Name', 'Stim/Recovery Split Check', ...
               'NumberTitle', 'off', ...
               'Color', 'w');

        t_vib = (0:numel(vib)-1) / fs_vib;
        t_sig = (0:N-1) / fs_sig;

        subplot(5,1,1)
        plot(t_vib, vib)
        ylabel('vib')
        title('Original noisy vib')

        subplot(5,1,2)
        plot(t_vib, detectInfo.vib_sm)
        ylabel('smoothed')
        title('Smoothed vib')

        subplot(5,1,3)
        plot(t_vib, detectInfo.vib_env)
        hold on
        yline(detectInfo.threshold, 'r--')
        ylabel('envelope')
        title('Envelope-based ON/OFF detection')

        subplot(5,1,4)
        plot(t_sig, double(dig_aligned))
        ylabel('dig')
        title('Aligned digital mask on signal time base')

        subplot(5,1,5)
        hold on
        if ~isempty(x_stim)
            plot(find(stimMask), x(stimMask,1), '.')
        end
        if ~isempty(x_recovery)
            plot(find(recoveryMask), x(recoveryMask,1), '.')
        end
        ylabel('Ch 1')
        xlabel('Sample index')
        title('Stim and recovery samples from channel 1')
        legend({'Stim','Recovery'})
    end
end

function [dig_aligned, info] = alignDigitalToSignalEnvelope(vib, fs_vib, Nsig, fs_sig, threshold)

    vib = vib(:);
    Nv = numel(vib);

    % Light smoothing of raw vib
    smoothWin = max(5, round(0.01 * fs_vib));
    vib_sm = movmean(vib, smoothWin);

    % Envelope / energy estimate
    envWin = max(5, round(0.1 * fs_vib));   % 100 ms window
    vib_env = sqrt(movmean(vib_sm.^2, envWin));

    % Automatic threshold if not provided
    if nargin < 5 || isempty(threshold)
        lo = prctile(vib_env, 20);
        hi = prctile(vib_env, 80);
        threshold = (lo + hi) / 2;
    end

    % Initial binary mask from envelope
    dig_vib = vib_env > threshold;

    % Remove short flickers first
    minDurSec = 10;
    minDur = max(1, round(minDurSec * fs_vib));
    dig_vib = removeShortSegments(dig_vib, minDur);

    % Keep only the largest ON segment
    d = diff([false; dig_vib; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    if ~isempty(starts)
        lengths = stops - starts + 1;
        [~, idx] = max(lengths);

        roughStart = starts(idx);
        roughStop = stops(idx);

        % Refine ON segment edges using the actual envelope crossing
        searchPadSec = 2.0;  % look +/- 2 s around rough edges
        searchPad = round(searchPadSec * fs_vib);

        s1 = max(1, roughStart - searchPad);
        s2 = min(Nv, roughStart + searchPad);
        e1 = max(1, roughStop  - searchPad);
        e2 = min(Nv, roughStop  + searchPad);

        % Refine start = first threshold crossing in start search window
        startRegion = vib_env(s1:s2) > threshold;
        startCandidates = find(startRegion, 1, 'first');
        if ~isempty(startCandidates)
            refinedStart = s1 + startCandidates - 1;
        else
            refinedStart = roughStart;
        end

        % Refine stop = last threshold crossing in end search window
        endRegion = vib_env(e1:e2) > threshold;
        endCandidates = find(endRegion, 1, 'last');
        if ~isempty(endCandidates)
            refinedStop = e1 + endCandidates - 1;
        else
            refinedStop = roughStop;
        end

        cleanMask = false(size(dig_vib));
        if refinedStop >= refinedStart
            cleanMask(refinedStart:refinedStop) = true;
        end
        dig_vib = cleanMask;
    else
        dig_vib = false(size(dig_vib));
        roughStart = [];
        roughStop = [];
        refinedStart = [];
        refinedStop = [];
    end

    % Align to target signal time base
    t_vib = (0:Nv-1)' / fs_vib;
    t_sig = (0:Nsig-1)' / fs_sig;
    dig_aligned = interp1(t_vib, double(dig_vib), t_sig, 'previous', 0) > 0;

    info = struct();
    info.vib_sm = vib_sm;
    info.vib_env = vib_env;
    info.threshold = threshold;
    info.dig_vib = dig_vib;
    info.smoothWin = smoothWin;
    info.envWin = envWin;
    info.minDurSec = minDurSec;
    info.roughStart = roughStart;
    info.roughStop = roughStop;
    info.refinedStart = refinedStart;
    info.refinedStop = refinedStop;
end

function mask_out = removeShortSegments(mask, minLen)
    mask = logical(mask(:));
    mask_out = mask;

    % Remove short ON segments
    d = diff([false; mask_out; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for i = 1:numel(starts)
        if (stops(i) - starts(i) + 1) < minLen
            mask_out(starts(i):stops(i)) = false;
        end
    end

    % Remove short OFF gaps too
    invMask = ~mask_out;
    d = diff([false; invMask; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for i = 1:numel(starts)
        if (stops(i) - starts(i) + 1) < minLen
            mask_out(starts(i):stops(i)) = true;
        end
    end
end