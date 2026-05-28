function [x_stim, x_recovery, stimMask, recoveryMask, stimEndTimeSec] = ...
    splitStimRecoveryManual(x, fs, baseFileName, doPlot)
% splitStimRecoveryManual
% Manually split signal into stim and recovery by asking user for stim end time.
%
% INPUTS
%   x            : [N x nChan] signal matrix
%   fs           : sampling rate
%   baseFileName : base filename for saving, without suffix
%   doPlot       : true/false
%
% OUTPUTS
%   x_stim           : signal from start up to stim end time
%   x_recovery       : signal after stim end time
%   stimMask         : logical mask for stim samples
%   recoveryMask     : logical mask for recovery samples
%   stimEndTimeSec   : user-entered stim end time in seconds

    if nargin < 4 || isempty(doPlot)
        doPlot = true;
    end

    validateattributes(x, {'numeric'}, {'2d','nonempty','real','finite'});
    validateattributes(fs, {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(baseFileName, {'char','string'}, {'nonempty'});

    [N, ~] = size(x);
    t = (0:N-1)' / fs;
    baseFileName = char(baseFileName);

    stimEndTimeSec = promptStimEndTime(t(end));
    if isempty(stimEndTimeSec)
        error('Stim/recovery split was canceled by user.');
    end

    stimEndIdx = max(1, min(N, floor(stimEndTimeSec * fs) + 1));

    stimMask = false(N,1);
    stimMask(1:stimEndIdx) = true;
    recoveryMask = ~stimMask;

    x_stim = x(stimMask, :);
    x_recovery = x(recoveryMask, :);

    save([baseFileName '_stim.mat'], 'x_stim', 'stimMask', 'fs', 'stimEndTimeSec');
    save([baseFileName '_recovery.mat'], 'x_recovery', 'recoveryMask', 'fs', 'stimEndTimeSec');

    if doPlot
        figure('Name', 'Manual Stim/Recovery Split', ...
               'NumberTitle', 'off', ...
               'Color', 'w');

        plot(t, x(:,1), 'k'); hold on
        xline(stimEndTimeSec, 'r--', 'LineWidth', 1.5)
        xlabel('Time (s)')
        ylabel('Ch 1')
        title('Manual stim/recovery split on channel 1')
        legend({'Signal', 'Stim End'})
    end
end

function stimEndTimeSec = promptStimEndTime(maxTimeSec)
    stimEndTimeSec = [];

    while true
        answer = inputdlg( ...
            {sprintf('Enter stim end time in seconds (0 to %.3f):', maxTimeSec)}, ...
            'Stim End Time', ...
            [1 50], ...
            {'120'});

        if isempty(answer)
            return;
        end

        stimEndTimeSec = str2double(strtrim(answer{1}));

        if isnan(stimEndTimeSec) || ~isfinite(stimEndTimeSec) || ...
                stimEndTimeSec <= 0 || stimEndTimeSec >= maxTimeSec
            uiwait(errordlg( ...
                sprintf('Stim end time must be a number between 0 and %.3f seconds.', maxTimeSec), ...
                'Invalid Stim End Time', 'modal'));
            continue;
        end

        return;
    end
end