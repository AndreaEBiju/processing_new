function results = HR_BR_HRVAnalysis(data, fs, cutoff, order, folderpath, condition, chanidx, mode, blankIdx, edgeBufferSec)
% HR_BR_HRVAnalysis
%
% Analyze heart rate, breath rate, and HRV for either:
%   1) a single signal matrix
%   2) a stim/recovery struct

    if nargin < 8 || isempty(mode)
        mode = 'single';
    end
    if nargin < 9 || isempty(blankIdx)
        blankIdx = [];
    end
    if nargin < 10 || isempty(edgeBufferSec)
        edgeBufferSec = 0.25;
    end

    mode = lower(string(mode));
    results = struct();

    switch mode
        case "single"
            results = analyzeOneSignal(data, fs, cutoff, order, folderpath, ...
                condition, chanidx, blankIdx, edgeBufferSec);

        case "stim_rec"
            if ~isstruct(data) || ~isfield(data, 'stim') || ~isfield(data, 'recovery')
                error('For stim_recovery mode, data must have fields .stim and .recovery');
            end
            if ~isstruct(blankIdx) || ~isfield(blankIdx, 'stim') || ~isfield(blankIdx, 'recovery')
                error('For stim_recovery mode, blankIdx must have fields .stim and .recovery');
            end

            if isstruct(chanidx)
                stimChanIdx = chanidx.stim;
                recoveryChanIdx = chanidx.recovery;
            else
                stimChanIdx = chanidx;
                recoveryChanIdx = chanidx;
            end

%             results.stim = analyzeOneSignal(data.stim, fs, cutoff, order, folderpath, ...
%                 [char(string(condition)) '_stim'], stimChanIdx, blankIdx.stim, edgeBufferSec);

            results.recovery = analyzeOneSignal(data.recovery, fs, cutoff, order, folderpath, ...
                [char(string(condition)) '_recovery'], recoveryChanIdx, blankIdx.recovery, edgeBufferSec);

        otherwise
            error('Unknown mode: %s. Use ''single'' or ''stim_recovery''.', mode);
    end
end

function out = analyzeOneSignal(signal, fs, cutoff, order, folderpath, condition, chanidx, blankIdx, edgeBufferSec)

    condition = string(condition);
    condition = condition(1);
    folderpath = char(string(folderpath));
    conditionLabel = char(condition);

    validateattributes(signal, {'numeric'}, {'2d','nonempty'});
    if chanidx > size(signal,2)
        error('chanidx exceeds number of channels.');
    end
    if isempty(blankIdx)
        blankIdx = zeros(0,2);
    end

    x = signal(:, chanidx);
    N = length(x);
    t = (0:N-1)'/fs;

    % ---------------- Valid / invalid masks ----------------
    nanMask = isnan(x);
    invalidMask = nanMask;

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

    % ---------------- Filter heartbeat channel ----------------
    xFill = fillmissing(x, 'linear', 'EndValues', 'nearest');

    highFreqCutoff = cutoff/(fs/2);
    [z,p,k] = butter(order, highFreqCutoff, 'low');
    [sos,g] = zp2sos(z, p, k);
    yFilt = filtfilt(sos, g, xFill);
    heartBeatSeries = detrend(yFilt);

    heartBeatSeries(invalidMask) = NaN;

    % ---------------- Heart peaks ----------------
    [~, heartlocsRaw] = findpeaks(yFilt);

    validHeartPeakMask = ~invalidMask(heartlocsRaw) & ~edgeMask(heartlocsRaw);
    heartlocs = heartlocsRaw(validHeartPeakMask);

    heartPeakTrain = zeros(N,1);
    heartPeakTrain(heartlocs) = 1;
    heartPeakTrain(invalidMask) = 0;

    [heartRateSeries, validHeartWindowSec] = movingRateIgnoringInvalid(heartPeakTrain, invalidMask, fs, 60);
    avgHeartRate = mean(heartRateSeries, 'omitnan');

    fprintf('%s Average Heart Rate: %.2f bpm\n', conditionLabel, avgHeartRate);

    % ---------------- Breath peaks ----------------
    if isempty(heartlocs)
        br_locs_true = [];
        breathRateSeries = nan(max(N - round(fs*60) + 1, 1), 1);
        avgBreathRate = NaN;
    else
        candidateVals = abs(xFill(heartlocs));
        if all(candidateVals == 0) || isempty(candidateVals)
            br_locs_true = [];
        else
            [~, br_locs] = findpeaks(candidateVals, ...
                'MinPeakHeight', 0.1 * max(candidateVals));
            br_locs_true = heartlocs(br_locs);
        end

        br_locs_true = br_locs_true(~invalidMask(br_locs_true) & ~edgeMask(br_locs_true));

        breathPeakTrain = zeros(N,1);
        breathPeakTrain(br_locs_true) = 1;
        breathPeakTrain(invalidMask) = 0;

        [breathRateSeries, validBreathWindowSec] = movingRateIgnoringInvalid(breathPeakTrain, invalidMask, fs, 60);
        avgBreathRate = mean(breathRateSeries, 'omitnan');
    end

    fprintf('%s Average Breath Rate: %.2f bpm\n', conditionLabel, avgBreathRate);

    % ---------------- RR intervals and times ----------------
    [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs);

    % Moving HRV series (event-based, gap-aware)
    [hrv_t, hrv_series, nRR_used] = movingHRVIgnoringGaps( ...
        RR_intervals, RR_times, invalidMask, fs, 60, 1, 3);

    if isempty(RR_intervals) || numel(RR_intervals) < 2
        hrv = NaN;
        rmssd = NaN;
        pnn5 = NaN;
        sd1 = NaN;
        sd2 = NaN;
        sampEn = NaN;
        appxEn = NaN;
        alpha1 = NaN;
        alpha2 = NaN;
        nVals = [];
        F = [];
    else
        hrv = std(RR_intervals, 'omitnan');
        fprintf('%s Heart Rate Variability (HRV): %.2f s\n', conditionLabel, hrv);

        diffRR = diff(RR_intervals);
        rmssd = sqrt(mean(diffRR.^2, 'omitnan'));
        fprintf('%s RMSSD: %.2f s\n', conditionLabel, rmssd);

        absdiffRR = abs(diff(RR_intervals .* 1000));
        count_over_5 = sum(absdiffRR > 5, 'omitnan');
        pnn5 = (count_over_5 / numel(absdiffRR)) * 100;
        fprintf('%s pNN5: %.2f%%\n', conditionLabel, pnn5);

        sd1 = sqrt(0.5) * std(diffRR, 'omitnan');
        sd2 = sqrt(max(0, 2*hrv^2 - sd1^2));
        fprintf('%s Poincare sd1: %.2f s\n', conditionLabel, sd1);
        fprintf('%s Poincare sd2: %.2f s\n', conditionLabel, sd2);

        sampEn = sampleEntropy(RR_intervals);
        fprintf('%s Sample entropy: %.4f\n', conditionLabel, sampEn);

        appxEn = approximateEntropy(RR_intervals);
        fprintf('%s Approximate entropy: %.4f\n', conditionLabel, appxEn);

        [alpha1, alpha2, nVals, F] = dfaRR(RR_intervals);
        fprintf('%s DFA alpha1: %.4f\n', conditionLabel, alpha1);
        fprintf('%s DFA alpha2: %.4f\n', conditionLabel, alpha2);
    end

    % ---------------- Save outputs ----------------
    hrbrFile = fullfile(folderpath, sprintf('%s_HRBR.mat', conditionLabel));
    hrvFile  = fullfile(folderpath, sprintf('%s_HRVMeasures.mat', conditionLabel));

    save(hrbrFile, ...
        'chanidx', 'heartBeatSeries', 'heartlocs', ...
        'heartRateSeries', 'avgHeartRate', ...
        'breathRateSeries', 'avgBreathRate', ...
        'invalidMask', 'edgeMask', 'blankIdx', 'edgeBufferSec', ...
        't');

    save(hrvFile, ...
        'hrv', 'rmssd', 'pnn5', 'sd1', 'sd2', ...
        'sampEn', 'appxEn', 'alpha1', 'alpha2', 'nVals', 'F', ...
        'RR_intervals', 'RR_times', 'heartlocs', ...
        'invalidMask', 'edgeMask', 'blankIdx', 'edgeBufferSec', ...
        'hrv_t', 'hrv_series', 'nRR_used', 't');

    % ---------------- Output struct ----------------
    out = struct();
    out.chanidx = chanidx;
    out.heartBeatSeries = heartBeatSeries;
    out.heartlocs = heartlocs;
    out.heartRateSeries = heartRateSeries;
    out.avgHeartRate = avgHeartRate;
    out.breathRateSeries = breathRateSeries;
    out.avgBreathRate = avgBreathRate;
    out.RR_intervals = RR_intervals;
    out.RR_times = RR_times;
    out.hrv = hrv;
    out.rmssd = rmssd;
    out.pnn5 = pnn5;
    out.sd1 = sd1;
    out.sd2 = sd2;
    out.sampEn = sampEn;
    out.appxEn = appxEn;
    out.alpha1 = alpha1;
    out.alpha2 = alpha2;
    out.nVals = nVals;
    out.F = F;
    out.hrv_t = hrv_t;
    out.hrv_series = hrv_series;
    out.nRR_used = nRR_used;
    out.invalidMask = invalidMask;
    out.edgeMask = edgeMask;
    out.blankIdx = blankIdx;
    out.t = t;
end

function [rateSeries, validWindowSec] = movingRateIgnoringInvalid(peakTrain, invalidMask, fs, winSec)
% Compute rate per minute using only valid time within each window.

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

function [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs)
% Keep RR intervals only when the samples strictly between consecutive peaks
% contain no invalid values.
%
% RR_times are assigned to the first peak of each valid RR interval.

    if numel(heartlocs) < 2
        RR_intervals = [];
        RR_times = [];
        return;
    end

    csumInvalid = [0; cumsum(invalidMask(:))];
    RR_intervals = [];
    RR_times = [];

    for i = 1:numel(heartlocs)-1
        s1 = heartlocs(i);
        s2 = heartlocs(i+1);

        if s2 <= s1 + 1
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1) = s1 / fs; %#ok<AGROW>
            continue;
        end

        nInvalidBetween = csumInvalid(s2) - csumInvalid(s1+1);
        if nInvalidBetween == 0
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1) = s1 / fs; %#ok<AGROW>
        end
    end
end

function [hrv_t, hrv_series, nRR_used] = movingHRVIgnoringGaps(RR_intervals, RR_times, invalidMask, fs, winSec, stepSec, minRR)
% Compute moving HRV (std of RR intervals) over time while respecting blanked gaps.

    if nargin < 7 || isempty(minRR)
        minRR = 3;
    end

    if isempty(RR_intervals) || isempty(RR_times)
        hrv_t = [];
        hrv_series = [];
        nRR_used = [];
        return;
    end

    sigDurSec = (numel(invalidMask)-1) / fs;
    hrv_t = (0:stepSec:sigDurSec)';
    hrv_series = NaN(size(hrv_t));
    nRR_used = zeros(size(hrv_t));

    halfWin = winSec / 2;

    csumInvalid = [0; cumsum(invalidMask(:))];

    for i = 1:numel(hrv_t)
        t0 = max(0, hrv_t(i) - halfWin);
        t1 = min(sigDurSec, hrv_t(i) + halfWin);

        idx0 = max(1, floor(t0 * fs) + 1);
        idx1 = min(numel(invalidMask), ceil(t1 * fs) + 1);

        nInvalid = csumInvalid(idx1+1) - csumInvalid(idx0);
        if nInvalid > 0
            continue;
        end

        keep = RR_times >= t0 & RR_times <= t1;
        nRR_used(i) = sum(keep);

        if nRR_used(i) >= minRR
            hrv_series(i) = std(RR_intervals(keep), 'omitnan');
        end
    end
end