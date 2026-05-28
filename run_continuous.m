%% Standalone driver replicating the "continuous run" block of main_mod
% Runs HR/BR/HRV + slow wave analysis on every *_blankmotion.mat file pair
% found in subfolders of the current directory.

clc; clear; close all;

%% Discover folders (same logic as main_mod) ------------------------------
d = dir(fullfile(pwd, '**'));
isub = [d.isdir];
folderPaths = fullfile({d(isub).folder}, {d(isub).name})';
folderPaths = folderPaths(~endsWith(folderPaths, {'.', '..','detectsort'}));
folderNames = cellfun(@(p) string(p(max(strfind(p, filesep))+1:end)), folderPaths);
keepMask = contains(folderNames, ["bl","recovery","stim","Stim","stim_rec"]);
folderPaths = folderPaths(keepMask);
folderNames = folderNames(keepMask);

fprintf('Discovered %d candidate folder(s):\n', numel(folderNames));
for i = 1:numel(folderNames)
    fprintf('  [%d] %s\n', i, folderNames(i));
end

%% Continuous run loop ----------------------------------------------------
for folderidx = 1:numel(folderNames)

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
    clear dataHR dataSW blankIdx results swResults

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

    %% Parameters (same as main_mod) --------------------------------------
    cutoff       = 8;
    order        = 4;
    lowPassOn    = false;
    lowPassCutoff = 2;
    lowPassOrder = 4;
    window       = 10;
    figtrue      = true;

    if exist('stim_y','var') && exist('rec_y','var')

        % --- HR/BR/HRV: pass full 5-channel data, hrChanIdx picks the channel
        dataHR.stim     = stim_y;
        dataHR.recovery = rec_y;

        blankIdx.stim     = stim_idx;
        blankIdx.recovery = rec_idx;

        hrChan.stim     = stim_hrChanIdx;
        hrChan.recovery = rec_hrChanIdx;

        fprintf('\n>> HR/BR/HRV (stim_rec)\n');
        results = HR_BR_HRVAnalysis_new(dataHR, fs, cutoff, order, ...
            folderPath, condition, hrChan, 'stim_rec', blankIdx, ...
            0.75, 20, 1, figtrue);

        % --- Slow wave: pass stomach channels only (cols 3:5 per chanlabels)
        dataSW.stim     = stim_y(:, 3:5);
        dataSW.recovery = rec_y(:, 3:5);

        fprintf('\n>> Slow wave (stim_rec)\n');
        swResults = slowWaveAnalysis_new(dataSW, lowPassOn, lowPassCutoff, ...
            lowPassOrder, fs, window, figtrue, folderPath, condition, ...
            'stim_rec', blankIdx, 3);

    elseif exist('rec_y','var')

        condOut = char(condition + "_recovery");

        fprintf('\n>> HR/BR/HRV (single, recovery only)\n');
        results = HR_BR_HRVAnalysis_new( ...
            rec_y, fs, cutoff, order, folderPath, condOut, rec_hrChanIdx, ...
            'single', rec_idx, 0.75, 20, 1, figtrue);

        fprintf('\n>> Slow wave (single, recovery only)\n');
        swResults = slowWaveAnalysis_new(rec_y(:,3:5), lowPassOn, lowPassCutoff, ...
            lowPassOrder, fs, window, figtrue, folderPath, condOut, ...
            'single', rec_idx, 3);

    elseif exist('y','var')

        fprintf('\n>> HR/BR/HRV (single)\n');
        results = HR_BR_HRVAnalysis_new( ...
            y, fs, cutoff, order, folderPath, char(condition), hrChanIdx, ...
            'single', idx, 0.75, 20, 1, figtrue);

        fprintf('\n>> Slow wave (single)\n');
        swResults = slowWaveAnalysis_new(y(:,3:5), lowPassOn, lowPassCutoff, ...
            lowPassOrder, fs, window, figtrue, folderPath, char(condition), ...
            'single', idx, 3);

    else
        warning('No usable blankmotion data loaded for %s. Skipping.', condition);
        continue;
    end

    % Persist top-level summaries for later inspection
    try
        save(fullfile(folderPath, condition + "_HR_BR_HRV_results.mat"), 'results');
        save(fullfile(folderPath, condition + "_slowWave_results.mat"), 'swResults');
    catch ME
        warning('Could not save summary results: %s', ME.message);
    end
end

fprintf('\nDone.\n');


%% ---- helper functions (copied verbatim from main_mod) ----

function idx = getBlankIdxFromFile(S, N)

    idx = zeros(0,2);
    if isfield(S, 'removedSegmentIdx');        idx = S.removedSegmentIdx;                       end
    if isfield(S, 'removedSegmentIdx_motion'); idx = [idx; S.removedSegmentIdx_motion];         end
    if isfield(S, 'removedSegmentIdx_stim');   idx = [idx; S.removedSegmentIdx_stim];           end

    idx = normalizeIdxLocal(idx, N);
    idx = mergeSegmentsLocal(idx);
end

function segments = getSegmentsFromFile(S)
    if isfield(S, 'removedSegments');        segments = S.removedSegments;
    elseif isfield(S, 'motionSegments');     segments = S.motionSegments;
    elseif isfield(S, 'stimSegments');       segments = S.stimSegments;
    elseif isfield(S, 'recoverySegments');   segments = S.recoverySegments;
    else;                                    segments = zeros(0,2);
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
    if isempty(idx); idx = zeros(0,2); return; end
    idx = round(idx);
    idx = idx(:,1:2);
    idx(:,1) = max(idx(:,1), 1);
    idx(:,2) = min(idx(:,2), N);
    bad = idx(:,2) < idx(:,1);
    idx(bad,:) = [];
end

function merged = mergeSegmentsLocal(idx)
    if isempty(idx); merged = zeros(0,2); return; end
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
