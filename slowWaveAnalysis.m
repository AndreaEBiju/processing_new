function results = slowWaveAnalysis(data, lowPassOn, lowPassCutoff, lowPassOrder, fs, ...
    window, figtrue, folderpath, condition, mode, blankIdx, edgeBufferSec)
% slowWaveAnalysis
%
% Extract slow waves from stomach EMG data for either:
%   1) a single signal matrix
%   2) a stim/recovery struct
%
% INPUTS
%   data          :
%                   - single mode: [N x nChan] matrix
%                   - stim_recovery mode: struct with fields .stim and .recovery
%   lowPassOn     : logical, apply low-pass filter
%   lowPassCutoff : scalar, low-pass cutoff frequency (Hz)
%   lowPassOrder  : scalar, low-pass filter order
%   fs            : scalar, sampling frequency
%   window        : scalar, smoothing window in seconds
%   figtrue       : logical, generate figure
%   folderpath    : folder to save outputs
%   condition     : base condition name for saving
%   mode          : 'single' or 'stim_recovery'
%   blankIdx      :
%                   - single mode: [nSeg x 2] blanked sample index ranges
%                   - stim_recovery mode: struct with fields .stim and .recovery
%   edgeBufferSec : remove peaks within this many seconds of blank edges
%                   default = 0.25 s
%
% OUTPUT
%   results       : struct containing outputs

    if nargin < 10 || isempty(mode)
        mode = 'single';
    end
    if nargin < 11 || isempty(blankIdx)
        blankIdx = [];
    end
    if nargin < 12 || isempty(edgeBufferSec)
        edgeBufferSec = 0.25;
    end

    mode = lower(string(mode));
    results = struct();

    switch mode
        case "single"
            results = analyzeOneSignal(data, lowPassOn, lowPassCutoff, lowPassOrder, ...
                fs, window, figtrue, folderpath, condition, blankIdx, edgeBufferSec);

        case "stim_rec"
            if ~isstruct(data) || ~isfield(data, 'stim') || ~isfield(data, 'recovery')
                error('For stim_recovery mode, data must have fields .stim and .recovery');
            end
            if ~isstruct(blankIdx) || ~isfield(blankIdx, 'stim') || ~isfield(blankIdx, 'recovery')
                error('For stim_recovery mode, blankIdx must have fields .stim and .recovery');
            end

            results.stim = analyzeOneSignal(data.stim, lowPassOn, lowPassCutoff, lowPassOrder, ...
                fs, window, figtrue, folderpath, [char(string(condition)) '_stim'], ...
                blankIdx.stim, edgeBufferSec);

            results.recovery = analyzeOneSignal(data.recovery, lowPassOn, lowPassCutoff, lowPassOrder, ...
                fs, window, figtrue, folderpath, [char(string(condition)) '_recovery'], ...
                blankIdx.recovery, edgeBufferSec);

        otherwise
            error('Unknown mode: %s. Use ''single'' or ''stim_recovery''.', mode);
    end
end

function out = analyzeOneSignal(signal, lowPassOn, lowPassCutoff, lowPassOrder, ...
    fs, window, figtrue, folderpath, condition, blankIdx, edgeBufferSec)

    condition = string(condition);
    condition = condition(1);
    conditionLabel = char(condition);
    folderpath = char(string(folderpath));

    validateattributes(signal, {'numeric'}, {'2d','nonempty'});
    N = size(signal,1);
    nchan = size(signal,2);
    t = (0:N-1)'/fs;

    if isempty(blankIdx)
        blankIdx = zeros(0,2);
    end

    % ---------------- invalid / edge masks ----------------
    invalidMask = any(isnan(signal), 2);

    if ~isempty(blankIdx)
        blankIdx = round(blankIdx);
        blankIdx(:,1) = max(blankIdx(:,1), 1);
        blankIdx(:,2) = min(blankIdx(:,2), N);

        for k = 1:size(blankIdx,1)
            invalidMask(blankIdx(k,1):blankIdx(k,2)) = true;
        end
    end

    edgeBufferSamp = round(edgeBufferSec * fs);
    edgeMask = false(N,1);

    for k = 1:size(blankIdx,1)
        s1 = max(1, blankIdx(k,1) - edgeBufferSamp);
        s2 = min(N, blankIdx(k,2) + edgeBufferSamp);
        edgeMask(s1:s2) = true;
    end

    % ---------------- filtering ----------------
    xFill = signal;
    for ch = 1:nchan
        xFill(:,ch) = fillmissing(xFill(:,ch), 'linear', 'EndValues', 'nearest');
    end

    if lowPassOn
        Wn = lowPassCutoff/(fs/2);
        [z,p,k] = butter(lowPassOrder, Wn, 'low');
        [sos,g] = zp2sos(z, p, k);
        filtSignal = filtfilt(sos, g, xFill);
        filteredSignal = detrend(filtSignal);
    else
        filteredSignal = xFill;
    end

    % Put NaNs back into invalid regions for bookkeeping/display
    filteredSignal(invalidMask,:) = NaN;

    % ---------------- smoothing ----------------
    windowlen = max(1, floor(window * fs));
    slowWaveTimeSeries = smoothdata(filteredSignal, 1, 'gaussian', windowlen);
    slowWaveTimeSeries(invalidMask,:) = NaN;

    % ---------------- outputs allocation ----------------
    rateWinSec = 60;
    rateWinSamp = max(1, round(rateWinSec * fs));
    outLen = max(N - rateWinSamp + 1, 1);

    slowWaveRateSeries = NaN(outLen, nchan);
    avgSlowWave = NaN(1, nchan);
    slowWavePeakLocs = cell(1, nchan);
    validWindowSec = NaN(outLen, nchan);

    if figtrue
        fig = figure('Name', sprintf('%s slow waves', conditionLabel), ...
            'NumberTitle', 'off', ...
            'Color', 'w');
        ax = gobjects(nchan,1);
    end

    % ---------------- peak detection per channel ----------------
    for i = 1:nchan
        y = slowWaveTimeSeries(:,i);

        % For peak detection, fill invalid areas temporarily, then reject peaks later
        yFill = fillmissing(y, 'linear', 'EndValues', 'nearest');

        [~, locsRaw] = findpeaks(yFill, 'MinPeakProminence', 5e-7);

        validPeakMask = ~invalidMask(locsRaw) & ~edgeMask(locsRaw);
        locs = locsRaw(validPeakMask);
        slowWavePeakLocs{i} = locs;

        slowWaveTrain = zeros(N,1);
        slowWaveTrain(locs) = 1;
        slowWaveTrain(invalidMask) = 0;

        [slowWaveRateSeries(:,i), validWindowSec(:,i)] = ...
            movingRateIgnoringInvalid(slowWaveTrain, invalidMask, fs, rateWinSec);

        avgSlowWave(i) = mean(slowWaveRateSeries(:,i), 'omitnan');
        fprintf('%s average slow wave rate at channel %d: %.2f per min\n', ...
            conditionLabel, i, avgSlowWave(i));

        if figtrue
            ax(i) = subplot(nchan,1,i, 'Parent', fig); hold(ax(i), 'on');
            plot(ax(i), t, slowWaveTimeSeries(:,i));
            plot(ax(i), t(locs), slowWaveTimeSeries(locs,i), 'v');
            title(ax(i), sprintf('Slow wave timeseries - %d', i));
            ylabel(ax(i), 'y (mV)');
            grid(ax(i), 'on');
        end
    end

    if figtrue
        xlabel(ax(end), 'Time (s)');
        linkaxes(ax, 'x');
        savefig(fig, fullfile(folderpath, sprintf('%s_slowwave.fig', conditionLabel)));
        close(fig);
    end

    swFile = fullfile(folderpath, sprintf('%s_slowWaves.mat', conditionLabel));
    save(swFile, ...
        'slowWaveTimeSeries', 'slowWaveRateSeries', 'avgSlowWave', ...
        'slowWavePeakLocs', 'invalidMask', 'edgeMask', 'blankIdx', ...
        'validWindowSec', 't', 'fs');

    out = struct();
    out.slowWaveTimeSeries = slowWaveTimeSeries;
    out.slowWaveRateSeries = slowWaveRateSeries;
    out.avgSlowWave = avgSlowWave;
    out.slowWavePeakLocs = slowWavePeakLocs;
    out.invalidMask = invalidMask;
    out.edgeMask = edgeMask;
    out.blankIdx = blankIdx;
    out.validWindowSec = validWindowSec;
    out.t = t;
    out.fs = fs;
end

function [rateSeries, validWindowSec] = movingRateIgnoringInvalid(peakTrain, invalidMask, fs, winSec)
% Compute per-minute rate using only valid time in each window.

    N = numel(peakTrain);
    winSamp = max(1, round(winSec * fs));
    validMask = ~invalidMask;

    if N < winSamp
        counts = sum(peakTrain(validMask));
        validSamples = sum(validMask);
        validWindowSec = validSamples / fs;
        if validSamples == 0
            rateSeries = NaN;
        else
            rateSeries = counts / validWindowSec * 60;
        end
        return;
    end

    counts = movsum(peakTrain, [winSamp-1 0], 'Endpoints', 'discard');
    validSamples = movsum(double(validMask), [winSamp-1 0], 'Endpoints', 'discard');

    validWindowSec = validSamples / fs;
    rateSeries = counts ./ validWindowSec * 60;
    rateSeries(validSamples == 0) = NaN;
end