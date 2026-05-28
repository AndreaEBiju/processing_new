function results = detectSortNerveSpikesECAP(blankedMatFile, outDir, varargin)
% detectSortNerveSpikesECAP
% Detect nerve spikes from blanked data, sort them, compute firing rates,
% and optionally ECAPs. Plotting is adapted for sparse baseline data.

p = inputParser;
addRequired(p, 'blankedMatFile', @(x) ischar(x) || isstring(x));
addRequired(p, 'outDir', @(x) ischar(x) || isstring(x));

addParameter(p, 'NerveChannels', [1 2], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'ChannelLabels', {'Ch1','Ch2'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'BandpassLow', 300, @isscalar);
addParameter(p, 'BandpassHigh', 5000, @isscalar);
addParameter(p, 'FilterOrder', 3, @isscalar);
addParameter(p, 'DetectionPolarity', 'neg', @(x) ischar(x) || isstring(x));
addParameter(p, 'ThreshSigma', 6, @isscalar);
addParameter(p, 'MaxThreshSigma', 40, @isscalar);
addParameter(p, 'PreMs', 0.6, @isscalar);
addParameter(p, 'PostMs', 1.0, @isscalar);
addParameter(p, 'RefractoryMs', 1.5, @isscalar);
addParameter(p, 'EdgeBufferMs', 5, @isscalar);
addParameter(p, 'MinAmpUV', 8, @isscalar);
addParameter(p, 'MaxAmpUV', 150, @isscalar);
addParameter(p, 'MinWidthMs', 0.15, @isscalar);
addParameter(p, 'MaxWidthMs', 1.5, @isscalar);
addParameter(p, 'FRBinSec', 1, @isscalar);
addParameter(p, 'SmoothFRSec', 5, @isscalar);

addParameter(p, 'DoSorting', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'NumClusters', 3, @isscalar);
addParameter(p, 'NumPCs', 3, @isscalar);
addParameter(p, 'MinClusterSize', 10, @isscalar);

addParameter(p, 'BurstSmoothSec', 10, @isscalar);
addParameter(p, 'BurstPeakMinDistanceSec', 5, @isscalar);
addParameter(p, 'BurstPeakPromFrac', 0.25, @isscalar);
addParameter(p, 'MinSpikesForBurst', 500, @isscalar);
addParameter(p, 'MinMeanRateForBurst', 0.5, @isscalar); % spikes/s

addParameter(p, 'StimTimesSec', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'StimVector', [], @(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(p, 'StimThreshold', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'ECAPPreMs', 2, @isscalar);
addParameter(p, 'ECAPPostMs', 10, @isscalar);

% ECAP-specific artifact handling
addParameter(p, 'ECAPArtifactPreMs', 5, @isscalar);
addParameter(p, 'ECAPArtifactPostMs', 5, @isscalar);

% Optional: use blanked segment centers as stim times
addParameter(p, 'UseBlankSegmentsAsStimTimes', false, @(x) islogical(x) || isnumeric(x));

addParameter(p, 'MakePlots', true, @(x) islogical(x) || isnumeric(x));

parse(p, blankedMatFile, outDir, varargin{:});
opts = p.Results;

blankedMatFile = char(string(blankedMatFile));
outDir = char(string(outDir));
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

S = load(blankedMatFile);
if ~isfield(S, 'yOut') || ~isfield(S, 'fs')
    error('Input file must contain yOut and fs.');
end

y = double(S.yOut);
fs = S.fs;

if isfield(S, 't')
    t = S.t(:);
else
    t = (0:size(y,1)-1)'/fs;
end

if isfield(S, 'removedSegmentIdx')
    removedSegmentIdx = S.removedSegmentIdx;
else
    removedSegmentIdx = zeros(0,2);
end

stimTimesSec = opts.StimTimesSec(:);

if isempty(stimTimesSec) && ~isempty(opts.StimVector)
    stimTimesSec = detectStimTimes(opts.StimVector(:), fs, opts.StimThreshold);
end

% Use blanked segment centers as stim times if requested
% Only use this if removedSegmentIdx corresponds to stimulation artifacts.
if isfield(S, 'removedSegmentIdx_stim') && ~isempty(S.removedSegmentIdx_stim)
    stimSegmentIdx = S.removedSegmentIdx_stim;
else
    stimSegmentIdx = removedSegmentIdx;
end

if isempty(stimTimesSec) && logical(opts.UseBlankSegmentsAsStimTimes) && ~isempty(stimSegmentIdx)
    stimIdx = round(mean(stimSegmentIdx, 2));
    stimIdx = stimIdx(stimIdx >= 1 & stimIdx <= size(y,1));
    stimTimesSec = (stimIdx(:) - 1) / fs;
end

chanIdx = opts.NerveChannels(:)';
labels = cellstr(opts.ChannelLabels);
if numel(labels) ~= numel(chanIdx)
    error('ChannelLabels must match NerveChannels.');
end
if max(chanIdx) > size(y,2)
    error('NerveChannels exceeds columns in yOut.');
end

bpHigh = min(opts.BandpassHigh, fs/2 - 1);
[b,a] = butter(opts.FilterOrder, [opts.BandpassLow bpHigh]/(fs/2), 'bandpass');

results = struct([]);

for k = 1:numel(chanIdx)
    ch = chanIdx(k);
    label = labels{k};

    fprintf('\nProcessing channel %d (%s)\n', ch, label);

    x = y(:,ch);
    invalidMask = isnan(x);

    edgeMask = false(size(x));
    edgePad = round(opts.EdgeBufferMs * 1e-3 * fs);
    if ~isempty(removedSegmentIdx)
        segs = round(removedSegmentIdx);
        segs(:,1) = max(segs(:,1), 1);
        segs(:,2) = min(segs(:,2), numel(x));
        for i = 1:size(segs,1)
            s1 = max(1, segs(i,1) - edgePad);
            s2 = min(numel(x), segs(i,2) + edgePad);
            edgeMask(s1:s2) = true;
        end
    else
        edgeMask = movmax(double(invalidMask), [edgePad edgePad]) > 0;
    end

    if any(invalidMask)
        xFill = fillmissing(x, 'linear', 'EndValues', 'nearest');
    else
        xFill = x;
    end

    xf = filtfilt(b, a, xFill);
    xf(invalidMask) = NaN;

    quietMask = ~invalidMask & ~edgeMask;
    if nnz(quietMask) < max(1000, round(2*fs))
        quietMask = ~invalidMask;
    end

    sigma = median(abs(xf(quietMask)), 'omitnan') / 0.6745;
    if ~isfinite(sigma) || sigma <= 0
        sigma = std(xf(quietMask), 'omitnan');
    end
    if ~isfinite(sigma) || sigma <= 0
        error('Could not estimate noise sigma for channel %d.', ch);
    end

    thr = opts.ThreshSigma * sigma;
    thrMax = opts.MaxThreshSigma * sigma;

    npre = round(opts.PreMs * 1e-3 * fs);
    npost = round(opts.PostMs * 1e-3 * fs);
    refractory = round(opts.RefractoryMs * 1e-3 * fs);

    [locs, amps] = detectCandidates(xf, thr, refractory, opts.DetectionPolarity);

    fprintf('  candidates before masking: %d\n', numel(locs));

    good = locs >= 1 & locs <= numel(x);
    good = good & ~invalidMask(locs) & ~edgeMask(locs);
    locs = locs(good);
    amps = amps(good);

    good = abs(amps) < thrMax;
    locs = locs(good);
    amps = amps(good);

    good = locs > npre & locs < numel(x)-npost;
    locs = locs(good);
    amps = amps(good);

    fprintf('  candidates after masks:   %d\n', numel(locs));

    centers = [];
    waveforms = [];
    widths_ms = [];
    peakAmps_uv = [];

    for i = 1:numel(locs)
        win = xf(locs(i)-npre : locs(i)+npost);
        if all(~isfinite(win))
            continue
        end

        [~, peakIdx] = max(abs(win));
        center = locs(i) - npre + peakIdx - 1;

        if center <= npre || center >= numel(x)-npost
            continue
        end
        if invalidMask(center) || edgeMask(center)
            continue
        end

        idx = center-npre:center+npost;
        w = xf(idx);
        if any(~isfinite(w))
            continue
        end

        amp_uv = (max(w) - min(w)) * 1e6;
        if amp_uv < opts.MinAmpUV || amp_uv > opts.MaxAmpUV
            continue
        end

        wa = abs(w);
        pk = max(wa);
        if pk <= 0
            continue
        end

        above = find(wa >= 0.5 * pk);
        if isempty(above)
            continue
        end

        width_ms = (above(end) - above(1) + 1) / fs * 1e3;
        if width_ms < opts.MinWidthMs || width_ms > opts.MaxWidthMs
            continue
        end

        centers(end+1,1) = center; %#ok<AGROW>
        waveforms(end+1,:) = w; %#ok<AGROW>
        widths_ms(end+1,1) = width_ms; %#ok<AGROW>
        peakAmps_uv(end+1,1) = amp_uv; %#ok<AGROW>
    end

    if ~isempty(centers)
        [centers, keepIdx] = enforceRefractoryStrongest(centers, peakAmps_uv, refractory);
        waveforms = waveforms(keepIdx,:);
        widths_ms = widths_ms(keepIdx);
        peakAmps_uv = peakAmps_uv(keepIdx);
    end

    spikeTimes = (centers - 1) / fs;
    fprintf('  accepted spikes:         %d\n', numel(spikeTimes));

    if ~isempty(waveforms)
        meanWave = mean(waveforms,1);
        stdWave = std(waveforms,0,1);
        Vpp_uv = (max(meanWave)-min(meanWave))*1e6;
    else
        meanWave = [];
        stdWave = [];
        Vpp_uv = NaN;
    end

    quietMask2 = quietMask;
    pad = round(2e-3 * fs);
    for i = 1:numel(centers)
        s1 = max(1, centers(i)-npre-pad);
        s2 = min(numel(x), centers(i)+npost+pad);
        quietMask2(s1:s2) = false;
    end
    quietSig = xf(quietMask2);
    if numel(quietSig) < round(fs)
        quietSig = xf(quietMask);
    end
    noiseStd = median(abs(quietSig), 'omitnan') / 0.6745;
    if ~isfinite(noiseStd) || noiseStd <= 0
        noiseStd = std(quietSig, 'omitnan');
    end
    SNR = Vpp_uv / (2 * noiseStd * 1e6);

    isi_ms = diff(spikeTimes) * 1000;
    refracViolFrac = mean(isi_ms < opts.RefractoryMs);

    [fr_t, fr_raw, fr_smooth, validFrac] = computeMaskedFiringRate( ...
        centers, numel(x), fs, opts.FRBinSec, opts.SmoothFRSec, invalidMask | edgeMask);

    meanRate = numel(spikeTimes) / (sum(~(invalidMask | edgeMask)) / fs);
    doBurst = numel(spikeTimes) >= opts.MinSpikesForBurst && meanRate >= opts.MinMeanRateForBurst;

    if doBurst
        burst = computeBurstRate(fr_t, fr_smooth, ...
            opts.BurstSmoothSec, opts.BurstPeakMinDistanceSec, opts.BurstPeakPromFrac);
    else
        burst = emptyBurstStruct(fr_t);
    end

    res = struct();
    res.channel = ch;
    res.label = label;
    res.fs = fs;
    res.t = t;
    res.filtered = xf;
    res.invalidMask = invalidMask;
    res.edgeMask = edgeMask;
    res.spikeCenters = centers;
    res.spikeTimes = spikeTimes;
    res.waveforms = waveforms;
    res.meanWaveform = meanWave;
    res.stdWaveform = stdWave;
    res.widths_ms = widths_ms;
    res.peakAmps_uv = peakAmps_uv;
    res.nSpikes = numel(spikeTimes);
    res.Vpp_uv = Vpp_uv;
    res.SNR = SNR;
    res.isi_ms = isi_ms;
    res.refracViolFrac = refracViolFrac;
    res.fr_t = fr_t;
    res.fr_raw = fr_raw;
    res.fr_smooth = fr_smooth;
    res.validFrac = validFrac;
    res.meanRate_spk_per_s = meanRate;
    res.burst = burst;
    res.params = opts;

    if logical(opts.DoSorting) && size(waveforms,1) >= max(opts.MinClusterSize, opts.NumClusters)
        [clusterID, clusterInfo] = sortSpikesSimple(waveforms, spikeTimes, peakAmps_uv, ...
            opts.NumClusters, opts.NumPCs, opts.MinClusterSize);

        res.clusterID = clusterID;
        res.clusterInfo = clusterInfo;

        for c = 1:numel(clusterInfo)
            idxc = (clusterID == clusterInfo(c).id);
            [c_t, c_raw, c_smooth, c_valid] = computeMaskedFiringRate( ...
                centers(idxc), numel(x), fs, opts.FRBinSec, opts.SmoothFRSec, invalidMask | edgeMask);

            cMeanRate = sum(idxc) / (sum(~(invalidMask | edgeMask)) / fs);
            doClusterBurst = sum(idxc) >= opts.MinSpikesForBurst && cMeanRate >= opts.MinMeanRateForBurst;

            if doClusterBurst
                cBurst = computeBurstRate(c_t, c_smooth, ...
                    opts.BurstSmoothSec, opts.BurstPeakMinDistanceSec, opts.BurstPeakPromFrac);
            else
                cBurst = emptyBurstStruct(c_t);
            end

            res.clusterInfo(c).fr_t = c_t;
            res.clusterInfo(c).fr_raw = c_raw;
            res.clusterInfo(c).fr_smooth = c_smooth;
            res.clusterInfo(c).validFrac = c_valid;
            res.clusterInfo(c).meanRate_spk_per_s = cMeanRate;
            res.clusterInfo(c).burst = cBurst;
        end
    else
        res.clusterID = [];
        res.clusterInfo = struct([]);
    end

    if ~isempty(stimTimesSec)
        res.ecap = computeECAPs(xf, fs, stimTimesSec, invalidMask | edgeMask, ...
            opts.ECAPPreMs, opts.ECAPPostMs, ...
            opts.ECAPArtifactPreMs, opts.ECAPArtifactPostMs, label);
    else
        res.ecap = struct();
    end

    if isempty(results)
        results = res;
    else
        results(end+1) = res; %#ok<AGROW>
    end

    save(fullfile(outDir, sprintf('%s_results.mat', label)), 'res');

    if logical(opts.MakePlots)
        makeIntegratedPlot(res, outDir);
    end
end

save(fullfile(outDir, 'detectSortNerveSpikesECAP_allResults.mat'), 'results');
end

function [locs, amps] = detectCandidates(xf, thr, refractory, polarity)
x0 = xf;
x0(~isfinite(x0)) = 0;

switch lower(string(polarity))
    case "neg"
        [pks, locs] = findpeaks(-x0, 'MinPeakHeight', thr, 'MinPeakDistance', refractory);
        amps = -pks;
    case "pos"
        [amps, locs] = findpeaks(x0, 'MinPeakHeight', thr, 'MinPeakDistance', refractory);
    case "both"
        [pksPos, locsPos] = findpeaks(x0, 'MinPeakHeight', thr, 'MinPeakDistance', refractory);
        [pksNeg, locsNeg] = findpeaks(-x0, 'MinPeakHeight', thr, 'MinPeakDistance', refractory);
        locs = [locsPos; locsNeg];
        amps = [pksPos; -pksNeg];
        [locs, ord] = sort(locs);
        amps = amps(ord);
        [locs, amps] = dedupeNearbyPeaks(locs, amps, refractory);
    otherwise
        error('Unknown DetectionPolarity');
end
end

function [locsOut, ampsOut] = dedupeNearbyPeaks(locs, amps, refractory)
if isempty(locs)
    locsOut = locs; ampsOut = amps; return;
end
keep = true(size(locs));
i = 1;
while i < numel(locs)
    j = i + 1;
    group = i;
    while j <= numel(locs) && (locs(j) - locs(j-1) <= refractory)
        group(end+1) = j; %#ok<AGROW>
        j = j + 1;
    end
    [~, imax] = max(abs(amps(group)));
    loser = setdiff(group, group(imax));
    keep(loser) = false;
    i = j;
end
locsOut = locs(keep);
ampsOut = amps(keep);
end

function [centersOut, keepIdx] = enforceRefractoryStrongest(centers, strength, refractory)
if isempty(centers)
    centersOut = centers;
    keepIdx = false(size(centers));
    return;
end
[centersS, ord] = sort(centers);
strengthS = strength(ord);
keepS = true(size(centersS));
i = 1;
while i < numel(centersS)
    j = i + 1;
    group = i;
    while j <= numel(centersS) && (centersS(j) - centersS(j-1) <= refractory)
        group(end+1) = j; %#ok<AGROW>
        j = j + 1;
    end
    [~, imax] = max(strengthS(group));
    loser = setdiff(group, group(imax));
    keepS(loser) = false;
    i = j;
end
centersOut = centersS(keepS);
keepIdx = false(size(centers));
keepIdx(ord) = keepS;
end

function [bin_t, fr_raw, fr_smooth, validFrac] = computeMaskedFiringRate(spikeCenters, N, fs, binSec, smoothSec, invalidMask)
edges = 0:binSec:(N/fs);
if edges(end) < N/fs
    edges(end+1) = N/fs;
end

spikeTimes = (spikeCenters - 1) / fs;
fr_counts = histcounts(spikeTimes, edges);
nBins = numel(edges)-1;

validFrac = zeros(1, nBins);
fr_raw = nan(1, nBins);

for i = 1:nBins
    s1 = max(1, floor(edges(i)*fs)+1);
    s2 = min(N, floor(edges(i+1)*fs));
    if s2 < s1
        continue
    end
    validSamples = sum(~invalidMask(s1:s2));
    validFrac(i) = validSamples / (s2-s1+1);
    validSec = validSamples / fs;
    if validSec > 0
        fr_raw(i) = fr_counts(i) / validSec;
    end
end

bin_t = edges(1:end-1) + diff(edges)/2;
smoothBins = max(1, round(smoothSec / binSec));
fr_smooth = movmean(fr_raw, smoothBins, 'omitnan');
end

function burst = emptyBurstStruct(fr_t)
burst = struct();
burst.hasBurst = false;
burst.rate_per_min = NaN;
burst.peak_times = [];
burst.fr_env = nan(size(fr_t));
burst.fr_t = fr_t;
burst.interval_sec = [];
end

function burst = computeBurstRate(fr_t, fr_smooth, burstSmoothSec, peakMinDistanceSec, peakPromFrac)
burst = emptyBurstStruct(fr_t);

if isempty(fr_t) || isempty(fr_smooth) || all(~isfinite(fr_smooth))
    return;
end

if numel(fr_t) < 3
    return;
end

dt = median(diff(fr_t));
if ~isfinite(dt) || dt <= 0
    return;
end

smoothBins = max(1, round(burstSmoothSec / dt));
fr_env = movmean(fr_smooth, smoothBins, 'omitnan');
burst.fr_env = fr_env;

good = isfinite(fr_env);
if nnz(good) < 3
    return;
end

envUse = fr_env;
envUse(~good) = nanmedian(fr_env);

prom = peakPromFrac * max(envUse, [], 'omitnan');
if ~isfinite(prom) || prom <= 0
    prom = 0;
end

minDistBins = max(1, round(peakMinDistanceSec / dt));
[~, locs] = findpeaks(envUse, 'MinPeakDistance', minDistBins, 'MinPeakProminence', prom);

if numel(locs) < 2
    return;
end

peakTimes = fr_t(locs);
intervals = diff(peakTimes);

burst.hasBurst = true;
burst.peak_times = peakTimes;
burst.rate_per_min = mean(1 ./ intervals) * 60;
burst.interval_sec = intervals;
end

function [clusterID, clusterInfo] = sortSpikesSimple(waveforms, spikeTimes, peakAmps_uv, nClusters, nPCs, minClusterSize)
[~, score] = pca(waveforms);
nPCs = min(nPCs, size(score,2));
X = score(:,1:nPCs);

rng(1);
clusterID = kmeans(X, nClusters, 'Replicates', 10, 'MaxIter', 1000);

ids = unique(clusterID);
clusterInfo = struct([]);
newID = clusterID;
nextGood = 1;

for i = 1:numel(ids)
    cid = ids(i);
    idx = (clusterID == cid);
    if sum(idx) < minClusterSize
        newID(idx) = 0;
        continue;
    end

    mw = mean(waveforms(idx,:), 1);
    vw = max(mw) - min(mw);
    isi_ms = diff(spikeTimes(idx)) * 1000;
    if isempty(isi_ms)
        refr = NaN;
    else
        refr = mean(isi_ms < 1.5);
    end

    entry = struct();
    entry.id = nextGood;
    entry.orig_id = cid;
    entry.nSpikes = sum(idx);
    entry.meanWaveform = mw;
    entry.Vpp_uv = vw * 1e6;
    entry.refracViolFrac = refr;
    entry.meanAmp_uv = mean(peakAmps_uv(idx), 'omitnan');
    entry.spikeTimes = spikeTimes(idx);

    newID(idx) = nextGood;

    if isempty(clusterInfo)
        clusterInfo = entry;
    else
        clusterInfo(end+1) = entry; %#ok<AGROW>
    end
    nextGood = nextGood + 1;
end

clusterID = newID;
end

function stimTimesSec = detectStimTimes(stimVec, fs, stimThreshold)
stimVec = double(stimVec(:));
if isempty(stimThreshold)
    lo = prctile(stimVec, 10);
    hi = prctile(stimVec, 90);
    stimThreshold = (lo + hi) / 2;
end
dig = stimVec > stimThreshold;
rising = find(diff([0; dig]) == 1);
stimTimesSec = (rising - 1) / fs;
end

function ecap = computeECAPs(xf, fs, stimTimesSec, invalidMask, ...
    preMs, postMs, artifactPreMs, artifactPostMs, label)

ecap = struct();

if isempty(stimTimesSec)
    return;
end

preSamp = round(preMs * 1e-3 * fs);
postSamp = round(postMs * 1e-3 * fs);

artifactPreSamp = round(artifactPreMs * 1e-3 * fs);
artifactPostSamp = round(artifactPostMs * 1e-3 * fs);

stimIdx = round(stimTimesSec * fs) + 1;

valid = stimIdx > preSamp & stimIdx < numel(xf) - postSamp;
stimIdx = stimIdx(valid);
stimTimesSec = stimTimesSec(valid);

relSamp = -preSamp:postSamp;
t_ms = relSamp / fs * 1e3;

% This is the region where NaNs are allowed because they are stimulation artifact blanks.
artifactAllowed = relSamp >= -artifactPreSamp & relSamp <= artifactPostSamp;

snips = [];
usedTimes = [];

for i = 1:numel(stimIdx)
    idx = stimIdx(i)-preSamp : stimIdx(i)+postSamp;

    bad = invalidMask(idx);

    % Reject trial only if invalid samples occur outside the allowed artifact region.
    if any(bad(:) & ~artifactAllowed(:))
        continue;
    end

    w = xf(idx);

    % Keep blanked artifact portion as NaN so it does not affect mean ECAP.
    badRow = bad(:)';
    w(badRow & artifactAllowed) = NaN;

    if all(~isfinite(w))
        continue;
    end

    snips(end+1,:) = w; %#ok<AGROW>
    usedTimes(end+1,1) = stimTimesSec(i); %#ok<AGROW>
end

ecap.label = label;
ecap.t_ms = t_ms;
ecap.artifactWindow_ms = [-artifactPreMs artifactPostMs];
ecap.usedStimTimesSec = usedTimes;
ecap.nTriggers = size(snips,1);

if isempty(snips)
    ecap.mean = [];
    ecap.std = [];
    ecap.sem = [];
    ecap.traces = [];
    ecap.Vpp_uv = NaN;
    ecap.peakLatency_ms = NaN;
    ecap.troughLatency_ms = NaN;
    return;
end

ecap.traces = snips;
ecap.mean = mean(snips, 1, 'omitnan');
ecap.std = std(snips, 0, 1, 'omitnan');
ecap.sem = ecap.std ./ sqrt(sum(isfinite(snips), 1));

% Measure ECAP only after the artifact window.
measureMask = t_ms > artifactPostMs & t_ms <= postMs;

if any(measureMask) && any(isfinite(ecap.mean(measureMask)))
    m = ecap.mean(measureMask);
    tm = t_ms(measureMask);

    [mx, imax] = max(m);
    [mn, imin] = min(m);

    ecap.Vpp_uv = (mx - mn) * 1e6;
    ecap.peakLatency_ms = tm(imax);
    ecap.troughLatency_ms = tm(imin);
else
    ecap.Vpp_uv = NaN;
    ecap.peakLatency_ms = NaN;
    ecap.troughLatency_ms = NaN;
end
end

function makeIntegratedPlot(res, outDir)
showBurst = isfield(res, 'burst') && ~isempty(res.burst) && res.burst.hasBurst;

fig = figure('Color','w','Position',[100 100 1300 1100]);
tiledlayout(5,1,'Padding','compact','TileSpacing','compact');

sig_uv = res.filtered * 1e6;
t = res.t(:);

nexttile;
if isempty(res.fr_raw) || all(isnan(res.fr_raw))
    previewCenter = min(max(t)/2, max(t));
else
    [~, imax] = max(res.fr_smooth);
    previewCenter = res.fr_t(imax);
end
previewSec = 10;
t0 = max(0, previewCenter - previewSec/2);
t1 = min(t(end), t0 + previewSec);
idx0 = max(1, floor(t0*res.fs)+1);
idx1 = min(numel(t), floor(t1*res.fs)+1);

plot(t(idx0:idx1), sig_uv(idx0:idx1), 'b'); hold on;
spk = res.spikeCenters(res.spikeCenters >= idx0 & res.spikeCenters <= idx1);
if ~isempty(spk)
    plot(t(spk), sig_uv(spk), 'r.', 'MarkerSize', 8);
end
xlabel('Time (s)');
ylabel('Voltage (uV)');
title(sprintf('%s | n=%d | SNR=%.2f | Vpp=%.1f uV', ...
    res.label, res.nSpikes, res.SNR, res.Vpp_uv), 'Interpreter', 'none');

nexttile;
if ~isempty(res.waveforms)
    wf_uv = res.waveforms * 1e6;
    mean_uv = res.meanWaveform * 1e6;
    std_uv = res.stdWaveform * 1e6;
    t_ms = ((1:size(res.waveforms,2)) - ceil(size(res.waveforms,2)/2)) / res.fs * 1000;

    [~, ord] = sort(res.peakAmps_uv, 'descend');
    nShow = min(120, size(wf_uv,1));
    showIdx = ord(1:nShow);

    plot(t_ms, wf_uv(showIdx,:)', 'Color', [0.85 0.85 0.85]); hold on;
    plot(t_ms, mean_uv, 'k', 'LineWidth', 2);
    plot(t_ms, mean_uv + std_uv, '--', 'LineWidth', 1);
    plot(t_ms, mean_uv - std_uv, '--', 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    title(sprintf('Spike waveforms | refractory viol frac = %.4f', res.refracViolFrac));
else
    text(0.3,0.5,'No accepted spikes','Units','normalized');
    axis off;
end

nexttile;
if ~isempty(res.isi_ms)
    histogram(res.isi_ms, 0:0.5:20);
    xline(res.params.RefractoryMs, 'r--');
    xlabel('ISI (ms)');
    ylabel('Count');
    title('ISI histogram');
else
    text(0.3,0.5,'No ISI data','Units','normalized');
    axis off;
end

nexttile;
stem(res.spikeTimes, ones(size(res.spikeTimes)), '.', 'MarkerSize', 4); hold on;
yyaxis right
plot(res.fr_t, res.fr_smooth, 'k', 'LineWidth', 1.5);
ylabel('Smoothed FR (spikes/s)');
yyaxis left
ylabel('Spike events');
xlabel('Time (s)');
title('Spike raster and smoothed firing rate');

nexttile;
if showBurst
    plot(res.fr_t, res.fr_raw, '-', 'LineWidth', 1); hold on;
    plot(res.fr_t, res.fr_smooth, 'k', 'LineWidth', 2);
    plot(res.fr_t, res.burst.fr_env, '--', 'LineWidth', 1.5);
    for i = 1:numel(res.burst.peak_times)
        xline(res.burst.peak_times(i), ':');
    end
    xlabel('Time (s)');
    ylabel('Spikes/s');
    title(sprintf('Firing rate | burst rate = %.2f / min', res.burst.rate_per_min));
elseif isfield(res, 'clusterInfo') && ~isempty(res.clusterInfo)
    hold on;
    cmap = lines(numel(res.clusterInfo));
    for c = 1:numel(res.clusterInfo)
        plot(res.clusterInfo(c).fr_t, res.clusterInfo(c).fr_smooth, 'Color', cmap(c,:), 'LineWidth', 1.5);
    end
    xlabel('Time (s)');
    ylabel('Spikes/s');
    title('Per-cluster smoothed firing rates');
else
    yyaxis left
    plot(res.fr_t, res.fr_raw, '-', 'LineWidth', 1); hold on;
    plot(res.fr_t, res.fr_smooth, 'k', 'LineWidth', 2);
    ylabel('Spikes/s');
    yyaxis right
    plot(res.fr_t, res.validFrac, '--', 'LineWidth', 1);
    ylabel('Valid fraction');
    xlabel('Time (s)');
    title('Firing rate and valid-data fraction');
end
if isfield(res, 'ecap') && ~isempty(res.ecap) && ...
        isfield(res.ecap, 'mean') && ~isempty(res.ecap.mean)

    figE = figure('Color','w','Position',[100 100 900 500]);

    traces_uv = res.ecap.traces * 1e6;
    mean_uv = res.ecap.mean * 1e6;
    sem_uv = res.ecap.sem * 1e6;

    plot(res.ecap.t_ms, traces_uv', 'Color', [0.85 0.85 0.85]); hold on;
    plot(res.ecap.t_ms, mean_uv, 'k', 'LineWidth', 2);
    plot(res.ecap.t_ms, mean_uv + sem_uv, 'r--', 'LineWidth', 1);
    plot(res.ecap.t_ms, mean_uv - sem_uv, 'r--', 'LineWidth', 1);

    xline(0, 'b--');
    xline(res.ecap.artifactWindow_ms(1), 'm--');
    xline(res.ecap.artifactWindow_ms(2), 'm--');

    xlabel('Time from stimulation (ms)');
    ylabel('Voltage (uV)');
    title(sprintf('%s ECAP | n=%d | Vpp=%.2f uV', ...
        res.label, res.ecap.nTriggers, res.ecap.Vpp_uv), ...
        'Interpreter', 'none');

    grid on;

    exportgraphics(figE, fullfile(outDir, sprintf('%s_ECAP.png', res.label)));
    close(figE);
end
exportgraphics(fig, fullfile(outDir, sprintf('%s_integrated.png', res.label)));
close(fig);
end