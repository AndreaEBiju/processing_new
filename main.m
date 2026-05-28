%% Batch Preprocessing of GI Multimodal Stimulation Data
% last updated: 28 Apr 2026

%% Setup
clc; clear; close all;
% pub = pubfig_setup( ...
%     'FontName','Times New Roman', ...
%     'BaseFontSize',16, ...
%     'LineWidth',1.5, ...
%     'MarkerSize',10, ...
%     'PlotsDir','plots', ...
%     'EnableLaTeX',true);

% fprintf('Saving plots in folder: %s\n', pub.dir);

% p = get(0, "MonitorPositions");
% if size(p,1)>1
%     f.Position = p(2, :); % second display
%     set(0, "DefaultFigurePosition", f.Position);
% end

% to extract tdt data
electrode_inds = [1:2,17:19];
chanlabels = {'RVN','LVN','ANT1','ANT2','ANT3'};
SDKPATH = '/Users/andreaelizabethbiju/Library/CloudStorage/GoogleDrive-andreabiju@g.harvard.edu/Shared drives/BIONICs Lab Workspace/Project Folders/GEMS/Survivals/Chronic Rat A - 11-22-25/mstim sweep 11-26-25/ME_STIM-251126_2/TDTMatlabSDK'; % or whatever path you extracted the SDK zip into
addpath(genpath(SDKPATH));

%% Get all subfolder names in current directory
d = dir(fullfile(pwd, '**'));
isub = [d.isdir];
folderPaths = fullfile({d(isub).folder}, {d(isub).name})';

% remove '.' and '..'
folderPaths = folderPaths(~endsWith(folderPaths, {'.', '..'}));

% extract just names if needed
folderNames = cellfun(@(p) string(p(max(strfind(p, filesep))+1:end)), folderPaths);

% filter by keywords
keepMask = contains(folderNames, ["bl","recovery","stim","Stim","stim_rec"]);
folderPaths = folderPaths(keepMask);
folderNames = folderNames(keepMask);

%% Process each folder and load data

for folderidx = 1:length(folderNames)
    folderPath = folderPaths{folderidx};
    condition = string(folderNames{folderidx});
    data = TDTbin2mat(folderPath);
    
    % get tdt data thats most relevant - CHECK IF DATA STREAMS MATCH
    signal=data.streams.Raww.data(electrode_inds,:)';
    stim = data.streams.BiPl.data';
    vib = data.streams.adc1.data';
    fs = data.streams.Raww.fs;
    fs_vib = data.streams.adc1.fs;
    times = (0:size(signal,1)-1)/fs;
    times_vib = (0:size(vib,1)-1)/fs_vib;
    save(strcat(folderPath,'/',condition,'_sig'),"signal",'fs',"times");
    save(strcat(folderPath,'/',condition,'_vib'),"vib",'fs_vib',"times_vib");
    save(strcat(folderPath,'/',condition,'_stim'),'stim','fs','times');
    
    %% notch filter data
    notch_order = 4;
    notch_freq = 60; %Hz
    notch_width = 1; %Hz, single sided
    signal = double(signal);
    Wn = [notch_freq - notch_width, notch_freq + notch_width] / (fs / 2); % Normalize frequency
    [z,p,k] = butter(notch_order, Wn, 'stop');
    [sos,g] = zp2sos(z, p, k);
    y = filtfilt(sos, g, signal);
    y = detrend(y);   % optional
    save(strcat(folderPath,'/',condition,'_notched'),'y','fs','times');
end

clear y fs times

%% remove artifacts
for folderidx = 21:23%1:length(folderNames)
    folderPath = folderPaths{folderidx};
    condition = string(folderNames{folderidx});

    fprintf('\nPreparing %s ...\n', condition);

    notchedFile = fullfile(folderPath, condition + "_notched.mat");
    vibFile     = fullfile(folderPath, condition + "_vib.mat");

    if ~isfile(notchedFile)
        error('Missing file: %s', notchedFile);
    end
    if ~isfile(vibFile)
        error('Missing file: %s', vibFile);
    end

    tic
    M = matfile(notchedFile);
    y = M.y;
    fs = M.fs;

    M2 = matfile(vibFile);
    vib = M2.vib;
    fs_vib = M2.fs_vib;
    fprintf('Loaded data for %s in %.2f s\n', condition, toc);

    if contains(condition,'stim_rec')
        [x_stim, x_recovery, stimMask, recoveryMask, stimEndTimeSec] = ...
            splitStimRecoveryManual(y, fs, fullfile(folderPath, condition), true);
    
        save(fullfile(folderPath, condition + "_alignedDig.mat"), ...
            "stimMask", "recoveryMask", "fs", "stimEndTimeSec");
    
        dataSplit.stim = x_stim;
        dataSplit.recovery = x_recovery;
    
        [stimSegments, recoverySegments, stimOut, recoveryOut, stimIdx, recoveryIdx, stimHRChanIdx, recoveryHRChanIdx] = ...
            browseMotionArtifacts(dataSplit, fs, 10, 'stim_rec', fullfile(folderPath, condition));
        if isempty(stimSegments) && isempty(recoverySegments) && isempty(stimIdx) && isempty(recoveryIdx)
            error('Motion artifact browsing was aborted by user in %s.', folderPath);
        end

    elseif (contains(condition,'bl') || contains(condition,'recovery')) && ~contains(condition,'stim_rec')

        [stimSegments, recoverySegments, stimOut, recoveryOut, stimIdx, recoveryIdx, stimHRChanIdx, recoveryHRChanIdx] = ...
            browseMotionArtifacts(y, fs, 10, 'single', fullfile(folderPath, condition));

        if isempty(stimSegments) && isempty(stimIdx) && isempty(stimHRChanIdx)
            error('Motion artifact browsing was aborted by user in %s.', folderPath);
        end
    end

    clear y vib times times_vib dataSplit x_stim x_recovery stimMask recoveryMask dig_aligned detectInfo
end

clearvars -except chanlabels folderNames folderPaths
clc

%% continuous run: HR/HRV from blankmotion files
for folderidx = 21:23%[7:17,19:20]%1:length(folderNames)

    folderPath = folderPaths{folderidx};
    condition = string(folderNames{folderidx});

    fprintf('\n==============================\n');
    fprintf('Folder: %s\n', condition);
    fprintf('==============================\n');

    files = dir(fullfile(folderPath, '*_blankmotion.mat'));
    files = files(~startsWith({files.name}, '._'));

    clear stim_y stim_t stim_segments stim_idx stim_hrChanIdx
    clear rec_y rec_t rec_segments rec_idx rec_hrChanIdx
    clear y t segments idx hrChanIdx fs
    clear data blankIdx results swResults

    if isempty(files)
        fprintf('No blankmotion files found. Skipping.\n');
        continue;
    end

    for k = 1:numel(files)

        fname = files(k).name;
        fpath = fullfile(files(k).folder, fname);

        fprintf('Reading: %s\n', fname);

        S = load(fpath);

        if ~isfield(S, 'yOut') || ~isfield(S, 'fs')
            warning('Skipping %s because it lacks yOut or fs.', fname);
            continue;
        end

        fs = S.fs;
        this_y = S.yOut;

        if isfield(S, 't')
            this_t = S.t;
        elseif isfield(S, 'times')
            this_t = S.times(:);
        else
            this_t = (0:size(this_y,1)-1)' / fs;
        end

        this_idx = getBlankIdxFromFile(S, size(this_y,1));
        this_segments = getSegmentsFromFile(S);
        this_hrChanIdx = getHRChanIdxFromFile(S);

        if contains(fname, '_stim_blankmotion')
            stim_y = this_y;
            stim_t = this_t;
            stim_segments = this_segments;
            stim_idx = this_idx;
            stim_hrChanIdx = this_hrChanIdx;

        elseif contains(fname, '_recovery_blankmotion')
            rec_y = this_y;
            rec_t = this_t;
            rec_segments = this_segments;
            rec_idx = this_idx;
            rec_hrChanIdx = this_hrChanIdx;

        else
            y = this_y;
            t = this_t;
            segments = this_segments;
            idx = this_idx;
            hrChanIdx = this_hrChanIdx;
        end
    end

    cutoff = 8;
    order = 4;

    lowPassOn = false;
    lowPassCutoff = 2;
    lowPassOrder = 4;
    window = 10;
    figtrue = true;

    if exist('stim_y','var') && exist('rec_y','var')

        data.stim = stim_y;
        data.recovery = rec_y;

        blankIdx.stim = stim_idx;
        blankIdx.recovery = rec_idx;

        hrChanIdx.stim = stim_hrChanIdx;
        hrChanIdx.recovery = rec_hrChanIdx;

        results = HR_BR_HRVAnalysis( ...
            data, fs, cutoff, order, folderPath, condition, hrChanIdx, ...
            'stim_rec', blankIdx, 0.25);

        swResults = slowWaveAnalysis(data, lowPassOn, lowPassCutoff, ...
        lowPassOrder, fs, window, figtrue, folderPath, condition, ...
        'stim_rec', blankIdx, 0.25);

    elseif exist('rec_y','var')

        condOut = char(condition + "_recovery");

        results = HR_BR_HRVAnalysis( ...
            rec_y, fs, cutoff, order, folderPath, condOut, rec_hrChanIdx, ...
            'single', rec_idx, 0.25);

        swResults = slowWaveAnalysis(data, lowPassOn, lowPassCutoff, ...
        lowPassOrder, fs, window, figtrue, folderPath, condition, ...
        'stim_rec', blankIdx, 0.25);

    elseif exist('y','var')

        results = HR_BR_HRVAnalysis( ...
            y, fs, cutoff, order, folderPath, char(condition), hrChanIdx, ...
            'single', idx, 0.25);

        swResults = slowWaveAnalysis( ... 
        y, lowPassOn, lowPassCutoff, lowPassOrder, fs, ... 
        window, figtrue, folderPath, condition, ... 
        'single', idx, 0.25);

    else
        warning('No usable blankmotion data loaded for %s. Skipping.', condition);
        continue;
    end
end

%% helper functions for mixed old/new blankmotion files

function idx = getBlankIdxFromFile(S, N)

    idx = zeros(0,2);

    if isfield(S, 'removedSegmentIdx')
        idx = S.removedSegmentIdx;
    end

    if isfield(S, 'removedSegmentIdx_motion')
        idx = [idx; S.removedSegmentIdx_motion];
    end

    if isfield(S, 'removedSegmentIdx_stim')
        idx = [idx; S.removedSegmentIdx_stim];
    end

    idx = normalizeIdxLocal(idx, N);
    idx = mergeSegmentsLocal(idx);
end

function segments = getSegmentsFromFile(S)

    if isfield(S, 'removedSegments')
        segments = S.removedSegments;
    elseif isfield(S, 'motionSegments')
        segments = S.motionSegments;
    elseif isfield(S, 'stimSegments')
        segments = S.stimSegments;
    elseif isfield(S, 'recoverySegments')
        segments = S.recoverySegments;
    else
        segments = zeros(0,2);
    end
end

function hrChanIdx = getHRChanIdxFromFile(S)

    if isfield(S, 'hrChanIdx') && ~isempty(S.hrChanIdx)
        hrChanIdx = S.hrChanIdx;
    elseif isfield(S, 'stimHRChanIdx') && ~isempty(S.stimHRChanIdx)
        hrChanIdx = S.stimHRChanIdx;
    elseif isfield(S, 'recoveryHRChanIdx') && ~isempty(S.recoveryHRChanIdx)
        hrChanIdx = S.recoveryHRChanIdx;
    else
        hrChanIdx = 1;
        warning('No HR channel index found. Defaulting to channel 1.');
    end
end

function idx = normalizeIdxLocal(idx, N)

    if isempty(idx)
        idx = zeros(0,2);
        return;
    end

    idx = round(idx);
    idx = idx(:,1:2);

    idx(:,1) = max(idx(:,1), 1);
    idx(:,2) = min(idx(:,2), N);

    bad = idx(:,2) < idx(:,1);
    idx(bad,:) = [];
end

function merged = mergeSegmentsLocal(idx)

    if isempty(idx)
        merged = zeros(0,2);
        return;
    end

    idx = sortrows(idx, 1);
    merged = idx(1,:);

    for i = 2:size(idx,1)
        if idx(i,1) <= merged(end,2) + 1
            merged(end,2) = max(merged(end,2), idx(i,2));
        else
            merged(end+1,:) = idx(i,:); %#ok<AGROW>
        end
    end
end