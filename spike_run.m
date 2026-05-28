%% SETUP
close all; clear; clc;

rootDir = '/Users/andreaelizabethbiju/Library/CloudStorage/GoogleDrive-andreabiju@g.harvard.edu/Shared drives/BIONICs Lab Workspace/Project Folders/GEMS/Survivals/TDT/0505-0506-0507';

channelLabels = {'RVN','LVN'};
nerveChannels = [1 2];

%% FIND ALL USABLE BLANKMOTION FILES
allBlank = dir(fullfile(rootDir, '**', '*blankmotion.mat'));

if ~isempty(allBlank)
    keep = ~startsWith({allBlank.name}, '.');
    keep = keep & ~contains({allBlank.folder}, filesep + ".");
    allBlank = allBlank(keep);
end

fprintf('Found %d usable blankmotion files:\n', numel(allBlank));
for i = 1:numel(allBlank)
    fprintf('%3d  %s\n', i, fullfile(allBlank(i).folder, allBlank(i).name));
end

if isempty(allBlank)
    error('No usable *blankmotion.mat files found under %s', rootDir);
end

%% INITIALIZE SUMMARY STORAGE
summaryRows = struct( ...
    'dataset', {}, ...
    'blankedMatFile', {}, ...
    'outDir', {}, ...
    'channel', {}, ...
    'label', {}, ...
    'nSpikes', {}, ...
    'meanRate_spk_per_s', {}, ...
    'Vpp_uv', {}, ...
    'SNR', {}, ...
    'refracViolFrac', {}, ...
    'meanWidth_ms', {}, ...
    'medianWidth_ms', {}, ...
    'validDuration_s', {} );

allResults = {};
rowCounter = 1;

%% MAIN DETECTION LOOP
for f = 1:numel(allBlank)
    fname = allBlank(f).name;

    if startsWith(fname, '.')
        continue;
    end

    blankedMatFile = fullfile(allBlank(f).folder, fname);
    [~, baseName] = fileparts(blankedMatFile);
    datasetName = erase(baseName, '_blankmotion');

    outDir = fullfile(allBlank(f).folder, [datasetName '_detectsort']);

    fprintf('\n====================================================\n');
    fprintf('Running: %s\n', blankedMatFile);
    fprintf('Output : %s\n', outDir);
    fprintf('====================================================\n');

    try
        isStimFile = contains(datasetName, '_stim', 'IgnoreCase', true) && ...
             ~contains(datasetName, '_recovery', 'IgnoreCase', true);

        results = detectSortNerveSpikesECAP(blankedMatFile, outDir, ...
            'NerveChannels', nerveChannels, ...
            'ChannelLabels', channelLabels, ...
            'BandpassLow', 300, ...
            'BandpassHigh', 5000, ...
            'DetectionPolarity', 'neg', ...
            'ThreshSigma', 6, ...
            'MaxThreshSigma', 40, ...
            'PreMs', 0.6, ...
            'PostMs', 1.0, ...
            'RefractoryMs', 1.5, ...
            'EdgeBufferMs', 5, ...
            'MinAmpUV', 8, ...
            'MaxAmpUV', 150, ...
            'MinWidthMs', 0.15, ...
            'MaxWidthMs', 1.5, ...
            'FRBinSec', 5, ...
            'SmoothFRSec', 5, ...
            'DoSorting', true, ...
            'NumClusters', 3, ...
            'NumPCs', 3, ...
            'MinClusterSize', 10, ...
            'MinSpikesForBurst', 500, ...
            'MinMeanRateForBurst', 0.5, ...
            'UseBlankSegmentsAsStimTimes', isStimFile, ...
            'ECAPPreMs', 2, ...
            'ECAPPostMs', 10, ...
            'ECAPArtifactPreMs', 5, ...
            'ECAPArtifactPostMs', 5, ...
            'MakePlots', true);

        batchEntry = struct();
        batchEntry.dataset = datasetName;
        batchEntry.blankedMatFile = blankedMatFile;
        batchEntry.outDir = outDir;
        batchEntry.results = results;
        allResults{end+1} = batchEntry;

        for ch = 1:numel(results)
            res = results(ch);

            validMask = ~(res.invalidMask | res.edgeMask);
            validDuration_s = sum(validMask) / res.fs;

            if isempty(res.widths_ms)
                meanWidth = NaN;
                medianWidth = NaN;
            else
                meanWidth = mean(res.widths_ms, 'omitnan');
                medianWidth = median(res.widths_ms, 'omitnan');
            end

            summaryRows(rowCounter).dataset = datasetName;
            summaryRows(rowCounter).blankedMatFile = blankedMatFile;
            summaryRows(rowCounter).outDir = outDir;
            summaryRows(rowCounter).channel = res.channel;
            summaryRows(rowCounter).label = res.label;
            summaryRows(rowCounter).nSpikes = res.nSpikes;
            summaryRows(rowCounter).meanRate_spk_per_s = res.meanRate_spk_per_s;
            summaryRows(rowCounter).Vpp_uv = res.Vpp_uv;
            summaryRows(rowCounter).SNR = res.SNR;
            summaryRows(rowCounter).refracViolFrac = res.refracViolFrac;
            summaryRows(rowCounter).meanWidth_ms = meanWidth;
            summaryRows(rowCounter).medianWidth_ms = medianWidth;
            summaryRows(rowCounter).validDuration_s = validDuration_s;
            rowCounter = rowCounter + 1;
        end

    catch ME
        warning('Failed on %s\n%s', blankedMatFile, ME.getReport('basic'));
    end
end

if isempty(summaryRows)
    error('No datasets completed successfully.');
end

summaryTable = struct2table(summaryRows);

%% SAVE SUMMARY TABLES
summaryCsv = fullfile(rootDir, 'detectsort_summary_all_datasets.csv');
summaryMat = fullfile(rootDir, 'detectsort_summary_all_datasets.mat');

writetable(summaryTable, summaryCsv);
save(summaryMat, 'summaryTable', 'allResults');

fprintf('\nSaved summary CSV: %s\n', summaryCsv);
fprintf('Saved summary MAT: %s\n', summaryMat);

%% LOAD EXISTING RESULTS ONLY
% Run this section instead of the detection loop if you already have saved results.
% Comment it out if you just ran the full batch above.

% close all; clear; clc;
% rootDir = '/Volumes/PortableSSD/GEMS_recording_042026';
% load(fullfile(rootDir, 'detectsort_summary_all_datasets.mat'), ...
%     'summaryTable', 'allResults');

%% STANDARD SUMMARY PLOTS
% Safe to run by itself if summaryTable and allResults are already in workspace.

makeChannelOverallRateBarplots(summaryTable, rootDir, 1, 'RVN');
makeChannelOverallRateBarplots(summaryTable, rootDir, 2, 'LVN');

makeChannelStatPlots(summaryTable, rootDir, 1, 'RVN');
makeChannelStatPlots(summaryTable, rootDir, 2, 'LVN');

makeWaveformGrid(allResults, rootDir, 1, 'RVN');
makeWaveformGrid(allResults, rootDir, 2, 'LVN');

fprintf('\nFinished standard summary plots.\n');

%% BASELINE VS RECOVERY WINDOWED FIRING-RATE PLOTS
% Safe to run by itself if allResults is already in workspace.

makeBaselineRecoveryWindowPlots(allResults, rootDir, 1);
makeBaselineRecoveryWindowPlots(allResults, rootDir, 2);

fprintf('\nFinished baseline vs recovery windowed firing-rate plots.\n');

%% HELPER FUNCTIONS

function makeChannelOverallRateBarplots(summaryTable, rootDir, channelNumber, channelName)
    T = summaryTable(summaryTable.channel == channelNumber, :);
    if isempty(T)
        return;
    end

    [datasetNames, ia] = unique(T.dataset, 'stable');
    T = T(ia, :);

    fig = figure('Color','w','Position',[100 100 1500 550]);
    bar(T.meanRate_spk_per_s);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Mean firing rate (spikes/s)');
    title(sprintf('%s overall mean firing rate', channelName), 'Interpreter', 'none');
    grid on;

    exportgraphics(fig, fullfile(rootDir, sprintf('%s_mean_firing_rate_barplot.png', channelName)));
    close(fig);
end

function makeChannelStatPlots(summaryTable, rootDir, channelNumber, channelName)
    T = summaryTable(summaryTable.channel == channelNumber, :);
    if isempty(T)
        return;
    end

    [datasetNames, ia] = unique(T.dataset, 'stable');
    T = T(ia, :);

    fig = figure('Color','w','Position',[100 100 1500 900]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    bar(T.nSpikes);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Spike count');
    title(sprintf('%s spike count', channelName), 'Interpreter', 'none');
    grid on;

    nexttile;
    bar(T.Vpp_uv);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Vpp (uV)');
    title(sprintf('%s mean waveform Vpp', channelName), 'Interpreter', 'none');
    grid on;

    nexttile;
    bar(T.SNR);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('SNR');
    title(sprintf('%s SNR', channelName), 'Interpreter', 'none');
    grid on;

    nexttile;
    bar(T.meanWidth_ms);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Mean width (ms)');
    title(sprintf('%s mean spike width', channelName), 'Interpreter', 'none');
    grid on;

    exportgraphics(fig, fullfile(rootDir, sprintf('%s_other_stats.png', channelName)));
    close(fig);

    fig2 = figure('Color','w','Position',[100 100 1500 500]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    bar(T.refracViolFrac);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Refractory violation fraction');
    title(sprintf('%s refractory violation fraction', channelName), 'Interpreter', 'none');
    grid on;

    nexttile;
    bar(T.validDuration_s);
    xticks(1:numel(datasetNames));
    xticklabels(datasetNames);
    xtickangle(45);
    ylabel('Valid duration (s)');
    title(sprintf('%s valid analyzed duration', channelName), 'Interpreter', 'none');
    grid on;

    exportgraphics(fig2, fullfile(rootDir, sprintf('%s_refrac_and_valid_duration.png', channelName)));
    close(fig2);
end

function makeWaveformGrid(allResults, rootDir, channelIdx, channelName)
    if isempty(allResults)
        return;
    end

    nDatasets = numel(allResults);
    nCols = ceil(sqrt(nDatasets));
    nRows = ceil(nDatasets / nCols);

    fig = figure('Color','w','Position',[100 100 1700 1000]);
    tiledlayout(nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

    for i = 1:nDatasets
        nexttile;

        batchEntry = allResults{i};
        results = batchEntry.results;

        if channelIdx > numel(results)
            axis off;
            title(sprintf('%s (missing)', batchEntry.dataset), 'Interpreter', 'none');
            continue;
        end

        res = results(channelIdx);

        if isempty(res.waveforms)
            axis off;
            title(sprintf('%s | no spikes', batchEntry.dataset), 'Interpreter', 'none');
            continue;
        end

        wf_uv = res.waveforms * 1e6;
        mean_uv = res.meanWaveform * 1e6;

        if isempty(res.stdWaveform)
            std_uv = zeros(size(mean_uv));
        else
            std_uv = res.stdWaveform * 1e6;
        end

        t_ms = ((1:size(res.waveforms,2)) - ceil(size(res.waveforms,2)/2)) / res.fs * 1000;

        [~, ord] = sort(res.peakAmps_uv, 'descend');
        nShow = min(80, size(wf_uv,1));
        showIdx = ord(1:nShow);

        plot(t_ms, wf_uv(showIdx,:)', 'Color', [0.85 0.85 0.85]); hold on;
        plot(t_ms, mean_uv, 'k', 'LineWidth', 2);
        plot(t_ms, mean_uv + std_uv, '--', 'LineWidth', 1);
        plot(t_ms, mean_uv - std_uv, '--', 'LineWidth', 1);

        xlabel('ms');
        ylabel('uV');
        title(sprintf('%s\nn=%d | Vpp=%.1f uV', ...
            batchEntry.dataset, res.nSpikes, res.Vpp_uv), ...
            'Interpreter', 'none', 'FontSize', 9);
        grid on;
    end

    sgtitle(sprintf('%s spike waveforms across datasets', channelName), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(rootDir, sprintf('%s_waveform_grid.png', channelName)));
    close(fig);
end

function makeBaselineRecoveryWindowPlots(allResults, rootDir, ch)
    if isempty(allResults)
        return;
    end

    if ch == 1
        channelName = 'RVN';
    else
        channelName = 'LVN';
    end

    windowNames = {'Baseline','Rec 1 min','Rec 2 min','Rec 5 min','Rec 10 min','Rec full'};
    windowDurSec = [NaN, 60, 120, 300, 600, Inf];

    datasetMap = struct();

    for i = 1:numel(allResults)
        batchEntry = allResults{i};
        blankedMatFile = batchEntry.blankedMatFile;
        [~, fileBase] = fileparts(blankedMatFile);

        if contains(fileBase, '_stim_blankmotion')
            continue;
        end
        
        isBaseline = endsWith(fileBase, '_blankmotion') && ...
                     ~contains(fileBase, '_recovery_blankmotion') && ...
                     ~contains(fileBase, '_stim_blankmotion');
        
        isRecovery = contains(fileBase, '_recovery_blankmotion');        
        if ~(isBaseline || isRecovery)
            continue;
        end

        dsKey = makePairKey(fileBase);

        fn = matlab.lang.makeValidName(dsKey);

        if ~isfield(datasetMap, fn)
            entry = struct();
            entry.key = dsKey;
            entry.baseline = [];
            entry.recovery = [];
            datasetMap.(fn) = entry;
        end

        if isBaseline
            datasetMap.(fn).baseline = batchEntry;
        elseif isRecovery
            datasetMap.(fn).recovery = batchEntry;
        end
    end

    mapFields = fieldnames(datasetMap);
    matched = struct([]);
    mc = 1;

    for i = 1:numel(mapFields)
        E = datasetMap.(mapFields{i});
        if ~isempty(E.baseline) && ~isempty(E.recovery)
            matched(mc).key = E.key; %#ok<AGROW>
            matched(mc).baseline = E.baseline;
            matched(mc).recovery = E.recovery;
            mc = mc + 1;
        end
    end

    if isempty(matched)
        warning('No matched baseline/recovery pairs found for channel %d.', ch);
        return;
    end

    groupOrder = zeros(numel(matched),1);
    for i = 1:numel(matched)
        groupOrder(i) = groupOfKey(matched(i).key);
    end
    [~, ord] = sortrows([groupOrder, (1:numel(matched))']);
    matched = matched(ord);
    groupOrder = groupOrder(ord);

    nD = numel(matched);
    Y = nan(nD, numel(windowNames));
    datasetLabels = cell(nD,1);

    for i = 1:nD
        datasetLabels{i} = matched(i).key;

        bres = matched(i).baseline.results(ch);
        rres = matched(i).recovery.results(ch);

        Y(i,1) = computeMeanRateWithinWindow(bres, 0, Inf);
        for w = 2:numel(windowDurSec)
            Y(i,w) = computeMeanRateWithinWindow(rres, 0, windowDurSec(w));
        end
    end

    fig = figure('Color','w','Position',[100 100 1800 700]);
    bar(Y, 'grouped');
    xlabel('Dataset');
    ylabel('Mean firing rate (spikes/s)');
    title(sprintf('Channel %d (%s): Baseline vs Recovery Windowed Firing Rates', ...
        ch, channelName), 'Interpreter', 'none');

    xticks(1:nD);
    xticklabels(datasetLabels);
    xtickangle(45);
    legend(windowNames, 'Location', 'bestoutside');
    grid on;
    hold on;

    xlinePos = [];
    for i = 1:nD-1
        if groupOrder(i) ~= groupOrder(i+1)
            xlinePos(end+1) = i + 0.5; %#ok<AGROW>
        end
    end

    yl = ylim;
    for x = xlinePos
        plot([x x], yl, 'k--', 'LineWidth', 1);
    end
    ylim(yl);

    idx1 = find(groupOrder == 1);
    idx2 = find(groupOrder == 2);

    if ~isempty(idx1)
        text(mean(idx1), yl(2)*0.98, 'MOC', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold');
    end

    if ~isempty(idx2)
        text(mean(idx2), yl(2)*0.98, 'loll + fre', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold');
    end

    exportgraphics(fig, fullfile(rootDir, sprintf('Channel%d_baseline_recovery_windowed_firing_rate.png', ch)));
    close(fig);

    T = array2table(Y, 'VariableNames', matlab.lang.makeValidName(windowNames));
    T = addvars(T, datasetLabels, 'Before', 1, 'NewVariableNames', 'dataset');
    writetable(T, fullfile(rootDir, sprintf('Channel%d_baseline_recovery_windowed_firing_rate.csv', ch)));
end

function rate = computeMeanRateWithinWindow(res, tStart, tEnd)
    if isempty(res) || ~isfield(res, 'spikeTimes') || ~isfield(res, 'fs')
        rate = NaN;
        return;
    end

    t = res.t(:);
    fs = res.fs;

    if isempty(t)
        N = numel(res.filtered);
        t = (0:N-1)'/fs;
    end

    if isinf(tEnd)
        tEnd = t(end);
    else
        tEnd = min(tEnd, t(end));
    end

    if tEnd <= tStart
        rate = NaN;
        return;
    end

    sampleMask = (t >= tStart) & (t < tEnd);
    if ~any(sampleMask)
        rate = NaN;
        return;
    end

    validMask = ~(res.invalidMask(:) | res.edgeMask(:));
    validMask = validMask & sampleMask;

    validSec = sum(validMask) / fs;
    if validSec <= 0
        rate = NaN;
        return;
    end

    spk = res.spikeTimes(:);
    nSpikes = sum(spk >= tStart & spk < tEnd);

    rate = nSpikes / validSec;
end

function g = groupOfKey(key)
    k = lower(key);
    if contains(k, '_moc')
        g = 1;
    elseif contains(k, '_loll') || contains(k, '_fre')
        g = 2;
    else
        g = 3;
    end
end

function key = makePairKey(fileBase)
    key = fileBase;

    key = erase(key, '_blankmotion');
    key = erase(key, '_recovery');
    key = erase(key, '_stim');

    key = regexprep(key, '_bl_\d+$', '');
    key = regexprep(key, '_stim_rec_\d+$', '');
end