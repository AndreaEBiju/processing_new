function batch_process()
% BATCH_PROCESS  Interactive batch driver for HR/HRV and slow wave analysis
%
% Workflow
% --------
% 1) Popup: stim/recovery separation time in seconds. Applied to every queued
%    file whose name contains "stim_rec" — the script discards the first
%    sepTimeSec seconds (= "stim" segment) and analyses ONLY the remainder
%    (= "recovery" segment). Files without "stim_rec" in their name are
%    analysed in full.
% 2) UI: build a queue of .mat files. Each row has an editable "Animal"
%    column auto-parsed from the filename as the first letter of the first
%    underscore-token (case-insensitive). Example:
%        E10_FRE_E10_stim_rec_0030_recovery_blankmotion.mat  ->  animal "e"
%    Multiple files for the same animal can keep the same identifier (group
%    them) or be edited to separate identifiers.
% 3) Popup per unique animal: integer channel index for HR/HRV/SampEn
%    detection. The same channel is used for every recording belonging to
%    that animal.
% 4) Parallel processing: parfor across files, pool size = floor(maxNumCompThreads/2).
%    Each worker loads its file, slices recovery if applicable, then calls
%    HR_BR_HRVAnalysis_new + slowWaveAnalysis_new and saves outputs alongside
%    the source file. Errors on individual files do not stop the batch.
%
% Accepted input file types
%   *_blankmotion.mat : yOut [Nx5], fs, removedSegmentIdx, hrChanIdx, t  (preferred)
%   *_notched.mat     : y [Nx5], fs, times                              (no blankIdx)
%   anything else with a 2-D numeric matrix and fs                       (best effort)
%
% Outputs (saved next to each source file)
%   <basename>_HRBR.mat, <basename>_HRVMeasures.mat       (HR/HRV)
%   <basename>_slowWaves.mat                              (slow wave)
%   <basename>_detection.png/fig, <basename>_HRVmetrics.png/fig
%   <basename>_slowwave.png/fig

    %% ====================== Tunable parameters ======================
    P = struct();
    % HR / HRV
    P.hr_cutoff        = 8;      % function arg; lowpass is commented out
    P.hr_order         = 4;
    P.hr_edgeBufferSec = 0.75;
    P.hr_winSec        = 20;     % HRV / beats-per-window
    P.hr_hrBrWinSec    = 60;     % HR / BR
    P.hr_stepSec       = 1;

    % Slow wave
    P.sw_lowPassOn     = true;
    P.sw_lowPassCutoff = 0.15;   % Hz
    P.sw_lowPassOrder  = 2;
    P.sw_edgeBufferSec = 15;     % covers ~13.3s settling time
    P.sw_smoothWindow  = 5;      % seconds gaussian smooth

    P.figtrue          = true;

    % Parallel
    P.poolSize         = max(1, floor(maxNumCompThreads / 2));

    %% ====================== Step 1: separation time ======================
    sepTimeSec = ui_ask_separation_time();
    if isempty(sepTimeSec); fprintf('Cancelled.\n'); return; end
    fprintf('Stim/recovery separation: %g s (applied to all *stim_rec* files)\n', sepTimeSec);

    %% ====================== Step 2: build queue ======================
    [files, animals] = ui_build_queue();
    if isempty(files); fprintf('Cancelled.\n'); return; end
    fprintf('Queue: %d file(s).\n', numel(files));

    %% ====================== Step 3: per-animal hrChanIdx ======================
    chanMap = ui_ask_animal_channels(animals);
    if isempty(chanMap); fprintf('Cancelled.\n'); return; end

    %% ====================== Step 4: parallel processing ======================
    run_in_parallel(files, animals, chanMap, sepTimeSec, P);
end


% ==========================================================================
% ============================ UI Step 1 ===================================
% ==========================================================================
function sepTimeSec = ui_ask_separation_time()
    prompt   = {'Stim/recovery separation time (seconds):'};
    dlgtitle = 'Step 1 — Stim/recovery separation';
    dims     = [1 60];
    default  = {'120'};
    answer   = inputdlg(prompt, dlgtitle, dims, default);
    if isempty(answer)
        sepTimeSec = [];
        return;
    end
    sepTimeSec = str2double(answer{1});
    if isnan(sepTimeSec) || sepTimeSec < 0
        errordlg('Separation time must be a non-negative number.', 'Bad input');
        sepTimeSec = [];
    end
end


% ==========================================================================
% ============================ UI Step 2 ===================================
% ==========================================================================
function [files, animals] = ui_build_queue()
% Build a queue of files with editable animal-name column.

    files   = {};
    animals = {};

    fig = uifigure('Name', 'Step 2 — Build processing queue', ...
        'Position', [200 200 900 500]);

    grid = uigridlayout(fig, [3 1]);
    grid.RowHeight   = {'fit', '1x', 'fit'};
    grid.ColumnWidth = {'1x'};

    % --- Top: buttons ---
    btnRow = uipanel(grid, 'BorderType', 'none');
    btnRow.Layout.Row    = 1;
    btnRow.Layout.Column = 1;
    btnLayout = uigridlayout(btnRow, [1 4]);
    btnLayout.ColumnWidth = {120, 120, '1x', 120};
    btnLayout.Padding     = [0 4 0 4];

    btnAdd = uibutton(btnLayout, 'Text', 'Add files...');
    btnRem = uibutton(btnLayout, 'Text', 'Remove selected');
    uilabel(btnLayout, 'Text', '');
    btnGo  = uibutton(btnLayout, 'Text', 'Continue', ...
        'BackgroundColor', [0.7 0.9 0.7]);

    % --- Middle: table ---
    tbl = uitable(grid, ...
        'ColumnName',     {'File', 'Animal'}, ...
        'ColumnEditable', [false, true], ...
        'ColumnWidth',    {'auto', 80}, ...
        'Data',           cell(0, 2), ...
        'CellSelectionCallback', @(s,e) setappdata(s, 'sel', e.Indices));
    tbl.Layout.Row    = 2;
    tbl.Layout.Column = 1;

    % --- Bottom: status ---
    status = uilabel(grid, 'Text', 'No files added.', ...
        'HorizontalAlignment', 'left');
    status.Layout.Row    = 3;
    status.Layout.Column = 1;

    % --- Callbacks ---
    btnAdd.ButtonPushedFcn = @(s,e) onAdd();
    btnRem.ButtonPushedFcn = @(s,e) onRemove();
    btnGo.ButtonPushedFcn  = @(s,e) onContinue();
    fig.CloseRequestFcn    = @(s,e) onCancel();

    cancelled = false;

    uiwait(fig);

    if cancelled || isempty(tbl.Data)
        files   = {};
        animals = {};
    else
        files   = tbl.Data(:, 1);
        animals = tbl.Data(:, 2);
    end
    if isvalid(fig); delete(fig); end

    % ---- nested ----
    function onAdd()
        [names, folder] = uigetfile({'*.mat', 'MAT files'}, ...
            'Add files to queue', 'MultiSelect', 'on');
        if isequal(names, 0); return; end
        if ischar(names); names = {names}; end
        D = tbl.Data;
        for k = 1:numel(names)
            fpath = fullfile(folder, names{k});
            ani   = parse_animal_id(names{k});
            D(end+1, :) = {fpath, ani}; %#ok<AGROW>
        end
        tbl.Data = D;
        status.Text = sprintf('%d file(s) queued.', size(tbl.Data, 1));
    end

    function onRemove()
        sel = getappdata(tbl, 'sel');
        if isempty(sel); return; end
        rows = unique(sel(:, 1));
        D = tbl.Data;
        D(rows, :) = [];
        tbl.Data = D;
        status.Text = sprintf('%d file(s) queued.', size(tbl.Data, 1));
    end

    function onContinue()
        if isempty(tbl.Data)
            status.Text = 'Queue is empty — add files or close to cancel.';
            return;
        end
        uiresume(fig);
    end

    function onCancel()
        cancelled = true;
        uiresume(fig);
    end
end


function ani = parse_animal_id(filename)
% Pick first letter of the SECOND underscore-token, lowercased.
% Filenames are expected as xx_<animal name>_yy.mat
%   e.g. "E10_FRE_E10_stim_rec_0030_recovery_blankmotion.mat"
%        tokens = {'E10','FRE','E10','stim','rec',...}; animal = 'f'
    [~, base, ~] = fileparts(filename);
    tokens = strsplit(base, '_');
    if numel(tokens) >= 2 && ~isempty(tokens{2})
        ani = lower(tokens{2}(1));
    elseif ~isempty(tokens) && ~isempty(tokens{1})
        ani = lower(tokens{1}(1));      % fallback if no second token
    else
        ani = '?';
    end
end


% ==========================================================================
% ============================ UI Step 3 ===================================
% ==========================================================================
function chanMap = ui_ask_animal_channels(animals)
% For every unique animal (case-insensitive), ask for an integer channel index.

    unique_animals = unique(lower(animals));
    nA = numel(unique_animals);

    prompts = cell(nA, 1);
    defaults = cell(nA, 1);
    for k = 1:nA
        prompts{k}  = sprintf('Animal "%s" — hrChanIdx (1–5):', unique_animals{k});
        defaults{k} = '1';
    end

    answer = inputdlg(prompts, 'Step 3 — Per-animal HR channel', [1 50], defaults);
    if isempty(answer)
        chanMap = struct();
        return;
    end

    chanMap = struct();
    for k = 1:nA
        ch = str2double(answer{k});
        if isnan(ch) || ch < 1 || ch ~= round(ch)
            errordlg(sprintf('Bad channel index for animal "%s"', unique_animals{k}), ...
                'Bad input');
            chanMap = struct();
            return;
        end
        % Field names must be valid identifiers — prepend "a_" for single-letter keys
        chanMap.(['a_' unique_animals{k}]) = ch;
    end
end


% ==========================================================================
% ============================ Step 4: parallel ===========================
% ==========================================================================
function run_in_parallel(files, animals, chanMap, sepTimeSec, P)
    n = numel(files);
    fprintf('\nRunning %d file(s) with poolSize = %d.\n', n, P.poolSize);

    % Start pool if not already running
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= P.poolSize
        if ~isempty(pool); delete(pool); end
        try
            parpool('local', P.poolSize);
        catch ME
            warning('batch_process:noPool', ...
                'Could not start parpool: %s. Falling back to serial.', ME.message);
        end
    end

    % Progress via DataQueue
    q = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(q, @(msg) progressTick(msg));

    statuses = cell(n, 1);

    parfor i = 1:n
        try
            statuses{i} = process_one_file(files{i}, animals{i}, chanMap, sepTimeSec, P);
            send(q, sprintf('OK  [%d/%d] %s', i, n, files{i}));
        catch ME
            statuses{i} = sprintf('ERROR: %s', ME.message);
            send(q, sprintf('ERR [%d/%d] %s  -- %s', i, n, files{i}, ME.message));
        end
    end

    % Summary
    nOK  = sum(cellfun(@(s) ~startsWith(s, 'ERROR'), statuses));
    nErr = n - nOK;
    fprintf('\n========== Batch summary ==========\n');
    fprintf('  OK    : %d\n', nOK);
    fprintf('  Error : %d\n', nErr);
    if nErr > 0
        fprintf('  Errors:\n');
        for i = 1:n
            if startsWith(statuses{i}, 'ERROR')
                fprintf('    %s : %s\n', files{i}, statuses{i});
            end
        end
    end

    % ---- nested ----
    function progressTick(msg)
        nDone = nDone + 1;
        fprintf('  [%3d/%3d] %s\n', nDone, n, msg);
    end
end


function status = process_one_file(fpath, animal, chanMap, sepTimeSec, P)
% Process a single file. Saves outputs alongside the source.
    status = 'OK';

    if ~isfile(fpath)
        error('File not found: %s', fpath);
    end

    [folderPath, base, ~] = fileparts(fpath);
    keyAnimal = ['a_' lower(animal(1))];

    if ~isfield(chanMap, keyAnimal)
        error('No HR channel mapped for animal "%s"', animal);
    end
    hrChanIdx = chanMap.(keyAnimal);

    % --- Load and dispatch on contents ---
    S = load(fpath);
    [signal, fs, blankIdx] = extract_signal(S);
    if isempty(signal)
        error('Could not find a signal matrix in %s', base);
    end
    N0 = size(signal, 1);

    % --- Special handling for files with "stim_rec" in the name ---
    %   1. Files already named *_stim_blankmotion : skip entirely
    %      (these ARE the stim portion; user does not want stim analysed).
    %   2. Files already named *_recovery_blankmotion : analyse as-is.
    %   3. Other stim_rec files : split into stim+recovery blankmotion mats,
    %      save both alongside the source, then analyse the recovery slice.
    if contains(base, '_stim_blankmotion')
        fprintf('Skipping stim-only file %s\n', base);
        status = 'SKIPPED (stim-only file)';
        return;
    end

    if contains(base, 'stim_rec') && ~contains(base, '_recovery_blankmotion')
        [signal, blankIdx, condition] = split_and_save(S, signal, fs, blankIdx, ...
            folderPath, base, sepTimeSec);
        if size(signal, 1) < 1
            error('Recovery segment is empty after splitting %s', base);
        end
    else
        % Use the base name directly; the analysis functions will write
        % output files prefixed with this string.
        condition = base;
    end

    if hrChanIdx > size(signal, 2)
        error('hrChanIdx %d exceeds %d available channels', ...
            hrChanIdx, size(signal, 2));
    end

    %% ----------------- HR / HRV -----------------
    HR_BR_HRVAnalysis_new( ...
        signal, fs, P.hr_cutoff, P.hr_order, folderPath, condition, ...
        hrChanIdx, blankIdx, P.hr_edgeBufferSec, ...
        P.hr_winSec, P.hr_stepSec, P.figtrue, P.hr_hrBrWinSec);

    %% ----------------- Slow wave -----------------
    % Slow wave function expects stomach channels. Per chanlabels
    % {'RVN','LVN','ANT1','ANT2','ANT3'} those are cols 3:5.
    if size(signal, 2) >= 5
        swData = signal(:, 3:5);
    elseif size(signal, 2) >= 3
        swData = signal(:, 1:3);
    else
        warning('process_one_file:tooFewChans', ...
            'File %s has only %d channels; skipping slow wave step.', ...
            base, size(signal, 2));
        close all force; %#ok<CLALL>
        return;
    end

    slowWaveAnalysis_new( ...
        swData, P.sw_lowPassOn, P.sw_lowPassCutoff, P.sw_lowPassOrder, ...
        fs, P.sw_smoothWindow, P.figtrue, folderPath, condition, ...
        blankIdx, P.sw_edgeBufferSec);

    % Close any figures left behind to free worker memory
    close all force; %#ok<CLALL>
end


function [recSignal, recBlankIdx, recCondition] = split_and_save( ...
    S, signal, fs, blankIdx, folderPath, base, sepTimeSec)
% For a "stim_rec" source file: write out a *_stim_blankmotion.mat and a
% *_recovery_blankmotion.mat alongside the source, then return the recovery
% slice for downstream analysis.

    N = size(signal, 1);
    cutSamp = max(1, round(sepTimeSec * fs) + 1);
    if cutSamp >= N
        error('sepTimeSec=%g s exceeds file duration (%.1f s) for %s', ...
            sepTimeSec, N / fs, base);
    end

    % --- Determine a "stem" by stripping any known suffix ---
    stem = regexprep(base, ...
        '_(notched|blankmotion|recovery_blankmotion|stim_blankmotion)$', '');

    % --- STIM slice [1 .. cutSamp-1] ---
    stimSignal = signal(1 : cutSamp - 1, :);
    stimBlankIdx = trim_blank(blankIdx, 1, cutSamp - 1, 0);

    % --- RECOVERY slice [cutSamp .. N] ---
    recSignal = signal(cutSamp : N, :);
    recBlankIdx = trim_blank(blankIdx, cutSamp, N, cutSamp - 1);

    % --- Save both alongside the source ---
    stimPath = fullfile(folderPath, [stem '_stim_blankmotion.mat']);
    recPath  = fullfile(folderPath, [stem '_recovery_blankmotion.mat']);

    splitInfo = struct('sepTimeSec', sepTimeSec, 'sourceFile', base, ...
        'cutSamp', cutSamp); %#ok<NASGU>

    save_blankmotion(stimPath, stimSignal, fs, stimBlankIdx, S, 'stim',     splitInfo);
    save_blankmotion(recPath,  recSignal,  fs, recBlankIdx,  S, 'recovery', splitInfo);

    recCondition = [stem '_recovery'];
end


function bk = trim_blank(blankIdx, lo, hi, offset)
% Restrict blankIdx to samples in [lo, hi] and shift indices by `offset`.
    if isempty(blankIdx); bk = zeros(0, 2); return; end
    bk = blankIdx;
    bk(bk(:, 2) < lo, :) = [];
    bk(bk(:, 1) > hi, :) = [];
    bk(:, 1) = max(bk(:, 1), lo) - offset;
    bk(:, 2) = min(bk(:, 2), hi) - offset;
end


function save_blankmotion(outPath, yOut, fs, blankIdx, srcS, region, splitInfo) %#ok<INUSD,INUSL>
% Save a blankmotion-style .mat alongside the source.
    N = size(yOut, 1);
    t = (0 : N - 1)' / fs;                                            %#ok<NASGU>
    removedSegmentIdx = blankIdx;                                     %#ok<NASGU>
    blankingApplied   = ~isempty(blankIdx);                           %#ok<NASGU>
    removedSegments   = [];                                           %#ok<NASGU>

    % Carry over hrChanIdx from source if present (animal-level channel
    % is decided separately in batch_process).
    if isfield(srcS, 'hrChanIdx') && ~isempty(srcS.hrChanIdx)
        hrChanIdx = srcS.hrChanIdx; %#ok<NASGU>
    else
        hrChanIdx = NaN; %#ok<NASGU>
    end

    save(outPath, 'yOut', 'fs', 't', 'removedSegmentIdx', 'blankingApplied', ...
        'removedSegments', 'hrChanIdx', 'splitInfo', 'region');
end


function [signal, fs, blankIdx] = extract_signal(S)
% Robustly pull (signal, fs, blankIdx) out of one of the supported .mat layouts.
    signal = []; fs = []; blankIdx = zeros(0, 2);

    if isfield(S, 'fs');     fs = S.fs;
    elseif isfield(S, 'Fs'); fs = S.Fs;
    end

    if isfield(S, 'yOut')
        signal = double(S.yOut);
    elseif isfield(S, 'y')
        signal = double(S.y);
    elseif isfield(S, 'signal')
        signal = double(S.signal);
    end

    if isfield(S, 'removedSegmentIdx');        blankIdx = [blankIdx; S.removedSegmentIdx]; end
    if isfield(S, 'removedSegmentIdx_motion'); blankIdx = [blankIdx; S.removedSegmentIdx_motion]; end
    if isfield(S, 'removedSegmentIdx_stim');   blankIdx = [blankIdx; S.removedSegmentIdx_stim]; end

    if ~isempty(blankIdx)
        N = size(signal, 1);
        blankIdx = round(blankIdx);
        blankIdx(:, 1) = max(blankIdx(:, 1), 1);
        blankIdx(:, 2) = min(blankIdx(:, 2), N);
        blankIdx(blankIdx(:, 2) < blankIdx(:, 1), :) = [];
        % Merge overlapping/adjacent segments
        if size(blankIdx, 1) > 1
            blankIdx = sortrows(blankIdx, 1);
            merged = blankIdx(1, :);
            for i = 2:size(blankIdx, 1)
                if blankIdx(i, 1) <= merged(end, 2) + 1
                    merged(end, 2) = max(merged(end, 2), blankIdx(i, 2));
                else
                    merged(end+1, :) = blankIdx(i, :); %#ok<AGROW>
                end
            end
            blankIdx = merged;
        end
    end
end
