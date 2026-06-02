function batch_spike_detect()
% BATCH_SPIKE_DETECT  Parallel nerve-spike detection driver.
%
% Workflow
% --------
% 1) Popup: stim/recovery separation time. Applied to every queued file
%    whose name contains "stim_rec" — the stim portion (0..sepTime) is
%    split out and written as *_stim_blankmotion.mat alongside the source,
%    and the recovery portion is written as *_recovery_blankmotion.mat.
%    Only the recovery slice is sent to spike detection (mirrors
%    batch_process.m).
% 2) UI: build a queue of .mat files to process (single "File" column).
% 3) Parallel processing: parfor across files. Each iteration:
%      - splits stim_rec source files if needed (or skips already-split
%        *_stim_blankmotion files);
%      - calls detectSortNerveSpikesECAP on the recovery / standalone file;
%      - lets detectSortNerveSpikesECAP draw its individual per-file figures
%        (MakePlots = true).
%    NO combined cross-dataset plots are produced.
% 4) After processing finishes, a per-file per-channel summary CSV + MAT
%    is written alongside the first queued file.
%
% Parameter overrides live in the `P` struct at the top of the function.

    %% ====================== Tunable parameters ======================
    P = struct();

    % Channels
    P.nerveChannels      = [1 2];
    P.channelLabels      = {'RVN','LVN'};

    % Bandpass + detection
    P.bandpassLow        = 300;
    P.bandpassHigh       = 5000;
    P.detectPolarity     = 'neg';
    P.threshSigma        = 6;
    P.maxThreshSigma     = 40;
    P.preMs              = 0.6;
    P.postMs             = 1.0;
    P.refractoryMs       = 1.5;
    P.edgeBufferMs       = 5;
    P.minAmpUV           = 8;
    P.maxAmpUV           = 150;
    P.minWidthMs         = 0.15;
    P.maxWidthMs         = 1.5;

    % Firing rate
    P.frBinSec           = 5;
    P.smoothFRSec        = 5;

    % Sorting / clustering
    P.doSorting          = true;
    P.numClusters        = 3;
    P.numPCs             = 3;
    P.minClusterSize     = 10;

    % Burst detection
    P.minSpikesForBurst    = 500;
    P.minMeanRateForBurst  = 0.5;

    % ECAP
    P.ecapPreMs          = 2;
    P.ecapPostMs         = 10;
    P.ecapArtifactPreMs  = 5;
    P.ecapArtifactPostMs = 5;

    P.makePlots          = true;

    % Parallel
    P.poolSize           = max(1, floor(maxNumCompThreads / 2));

    %% ====================== Step 1: separation time ======================
    sepTimeSec = ui_ask_separation_time();
    if isempty(sepTimeSec); fprintf('Cancelled.\n'); return; end
    fprintf('Stim/recovery separation: %g s (applied to all *stim_rec* files)\n', sepTimeSec);

    %% ====================== Step 2: build queue ======================
    files = ui_build_simple_queue();
    if isempty(files); fprintf('Cancelled.\n'); return; end
    fprintf('Queue: %d file(s).\n', numel(files));

    %% ====================== Step 3: parallel processing ======================
    [summaryRows, allResults, statuses] = run_spike_parallel(files, sepTimeSec, P);

    %% ====================== Step 4: write summary ======================
    if ~isempty(summaryRows)
        outRoot = fileparts(files{1});
        if isempty(outRoot); outRoot = pwd; end
        summaryTable = struct2table(summaryRows);
        summaryCsv = fullfile(outRoot, 'detectsort_summary.csv');
        summaryMat = fullfile(outRoot, 'detectsort_summary.mat');
        try
            writetable(summaryTable, summaryCsv);
            save(summaryMat, 'summaryTable', 'allResults');
            fprintf('\nSaved summary CSV: %s\n', summaryCsv);
            fprintf('Saved summary MAT: %s\n', summaryMat);
        catch ME
            warning('batch_spike_detect:saveSummary', ...
                'Could not save summary: %s', ME.message);
        end
    else
        fprintf('\nNo successful runs — no summary saved.\n');
    end

    %% ====================== Batch summary ======================
    n    = numel(files);
    nOK  = sum(cellfun(@(s) startsWith(s, 'OK'), statuses));
    nSkp = sum(cellfun(@(s) startsWith(s, 'SKIPPED'), statuses));
    nErr = n - nOK - nSkp;
    fprintf('\n========== Batch summary ==========\n');
    fprintf('  OK      : %d\n', nOK);
    fprintf('  Skipped : %d (stim-only files)\n', nSkp);
    fprintf('  Errors  : %d\n', nErr);
    if nErr > 0
        for i = 1:n
            if startsWith(statuses{i}, 'ERROR')
                fprintf('    %s : %s\n', files{i}, statuses{i});
            end
        end
    end
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
    if isempty(answer); sepTimeSec = []; return; end
    sepTimeSec = str2double(answer{1});
    if isnan(sepTimeSec) || sepTimeSec < 0
        errordlg('Separation time must be a non-negative number.', 'Bad input');
        sepTimeSec = [];
    end
end


% ==========================================================================
% ============================ UI Step 2 ===================================
% ==========================================================================
function files = ui_build_simple_queue()
% File-queue UI with a single "File" column. Same look as batch_process's
% queue UI minus the Animal column.

    files = {};

    fig = uifigure('Name', 'Step 2 — Build spike-detection queue', ...
        'Position', [200 200 900 500]);

    grid = uigridlayout(fig, [3 1]);
    grid.RowHeight   = {'fit', '1x', 'fit'};
    grid.ColumnWidth = {'1x'};

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

    tbl = uitable(grid, ...
        'ColumnName',     {'File'}, ...
        'ColumnEditable', false, ...
        'ColumnWidth',    {'auto'}, ...
        'Data',           cell(0, 1), ...
        'CellSelectionCallback', @(s,e) setappdata(s, 'sel', e.Indices));
    tbl.Layout.Row    = 2;
    tbl.Layout.Column = 1;

    status = uilabel(grid, 'Text', 'No files added.', 'HorizontalAlignment', 'left');
    status.Layout.Row    = 3;
    status.Layout.Column = 1;

    btnAdd.ButtonPushedFcn = @(s,e) onAdd();
    btnRem.ButtonPushedFcn = @(s,e) onRemove();
    btnGo.ButtonPushedFcn  = @(s,e) onContinue();
    fig.CloseRequestFcn    = @(s,e) onCancel();

    cancelled = false;

    uiwait(fig);

    if cancelled || isempty(tbl.Data)
        files = {};
    else
        files = tbl.Data(:, 1);
    end
    if isvalid(fig); delete(fig); end

    function onAdd()
        [names, folder] = uigetfile({'*.mat', 'MAT files'}, ...
            'Add files to queue', 'MultiSelect', 'on');
        if isequal(names, 0); return; end
        if ischar(names); names = {names}; end
        D = tbl.Data;
        for k = 1:numel(names)
            D(end+1, 1) = {fullfile(folder, names{k})}; %#ok<AGROW>
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


% ==========================================================================
% ============================ Step 3: parallel ============================
% ==========================================================================
function [summaryRows, allResults, statuses] = run_spike_parallel(files, sepTimeSec, P)
    n = numel(files);
    fprintf('\nRunning %d file(s) with poolSize = %d.\n', n, P.poolSize);

    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= P.poolSize
        if ~isempty(pool); delete(pool); end
        try
            parpool('local', P.poolSize);
        catch ME
            warning('batch_spike_detect:noPool', ...
                'Could not start parpool: %s. Falling back to serial.', ME.message);
        end
    end

    q = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(q, @(msg) progressTick(msg));

    summaryCells = cell(n, 1);
    resultCells  = cell(n, 1);
    statuses     = cell(n, 1);

    parfor i = 1:n
        try
            [summaryCells{i}, resultCells{i}, statusMsg] = ...
                process_one_spike_file(files{i}, sepTimeSec, P);
            statuses{i} = statusMsg;
            send(q, sprintf('%s [%d/%d] %s', statusMsg, i, n, files{i}));
        catch ME
            statuses{i} = sprintf('ERROR: %s', ME.message);
            summaryCells{i} = struct([]);
            resultCells{i}  = {};
            send(q, sprintf('ERR [%d/%d] %s  -- %s', i, n, files{i}, ME.message));
        end
    end

    % Concatenate per-file summaries (each is a struct array of rows)
    summaryRows = struct([]);
    for i = 1:n
        s = summaryCells{i};
        if ~isempty(s)
            if isempty(summaryRows)
                summaryRows = s;
            else
                summaryRows = [summaryRows, s]; %#ok<AGROW>
            end
        end
    end

    % Keep the heavy per-file results in a cell array (mirrors spike_run.m)
    allResults = {};
    for i = 1:n
        if ~isempty(resultCells{i})
            allResults{end+1} = resultCells{i}; %#ok<AGROW>
        end
    end

    function progressTick(msg)
        nDone = nDone + 1;
        fprintf('  [%3d/%3d] %s\n', nDone, n, msg);
    end
end


function [summaryRows, batchEntry, statusMsg] = process_one_spike_file(fpath, sepTimeSec, P)
% Process a single file: optionally split, then call detectSortNerveSpikesECAP.
% Returns:
%   summaryRows : struct array of summary rows (one per channel, 0 if skipped)
%   batchEntry  : struct mirroring spike_run.m's allResults entries
%   statusMsg   : "OK ..." / "SKIPPED ..." / "ERROR ..."

    summaryRows = struct([]);
    batchEntry  = struct([]);
    statusMsg   = 'OK';

    if ~isfile(fpath)
        error('File not found: %s', fpath);
    end

    [folderPath, base, ~] = fileparts(fpath);

    % --- Decide path forward based on filename pattern ---
    if contains(base, '_stim_blankmotion')
        % Stim-only file → skip; user does not want stim processed.
        statusMsg = 'SKIPPED (stim-only file)';
        return;
    end

    if contains(base, 'stim_rec') && ~contains(base, '_recovery_blankmotion')
        % Split first, then operate on the recovery file.
        recPath = split_and_save_recovery(fpath, sepTimeSec);
        runPath = recPath;
    else
        runPath = fpath;
    end

    % --- Compose outDir like spike_run.m did ---
    [runFolder, runBase] = fileparts(runPath);
    datasetName = erase(runBase, '_blankmotion');
    outDir = fullfile(runFolder, [datasetName '_detectsort']);

    isStimFile = contains(datasetName, '_stim', 'IgnoreCase', true) && ...
                 ~contains(datasetName, '_recovery', 'IgnoreCase', true);

    % --- Run detection ---
    results = detectSortNerveSpikesECAP(runPath, outDir, ...
        'NerveChannels',                P.nerveChannels, ...
        'ChannelLabels',                P.channelLabels, ...
        'BandpassLow',                  P.bandpassLow, ...
        'BandpassHigh',                 P.bandpassHigh, ...
        'DetectionPolarity',            P.detectPolarity, ...
        'ThreshSigma',                  P.threshSigma, ...
        'MaxThreshSigma',               P.maxThreshSigma, ...
        'PreMs',                        P.preMs, ...
        'PostMs',                       P.postMs, ...
        'RefractoryMs',                 P.refractoryMs, ...
        'EdgeBufferMs',                 P.edgeBufferMs, ...
        'MinAmpUV',                     P.minAmpUV, ...
        'MaxAmpUV',                     P.maxAmpUV, ...
        'MinWidthMs',                   P.minWidthMs, ...
        'MaxWidthMs',                   P.maxWidthMs, ...
        'FRBinSec',                     P.frBinSec, ...
        'SmoothFRSec',                  P.smoothFRSec, ...
        'DoSorting',                    P.doSorting, ...
        'NumClusters',                  P.numClusters, ...
        'NumPCs',                       P.numPCs, ...
        'MinClusterSize',               P.minClusterSize, ...
        'MinSpikesForBurst',            P.minSpikesForBurst, ...
        'MinMeanRateForBurst',          P.minMeanRateForBurst, ...
        'UseBlankSegmentsAsStimTimes',  isStimFile, ...
        'ECAPPreMs',                    P.ecapPreMs, ...
        'ECAPPostMs',                   P.ecapPostMs, ...
        'ECAPArtifactPreMs',            P.ecapArtifactPreMs, ...
        'ECAPArtifactPostMs',           P.ecapArtifactPostMs, ...
        'MakePlots',                    P.makePlots);

    % --- Build the heavy results record (matches spike_run.m schema) ---
    batchEntry = struct();
    batchEntry.dataset        = datasetName;
    batchEntry.blankedMatFile = runPath;
    batchEntry.outDir         = outDir;
    batchEntry.results        = results;

    % --- Build per-channel summary rows ---
    nCh = numel(results);
    summaryRows = repmat(emptySummaryRow(), 1, nCh);
    for ch = 1:nCh
        res = results(ch);
        validMask       = ~(res.invalidMask | res.edgeMask);
        validDuration_s = sum(validMask) / res.fs;
        if isempty(res.widths_ms)
            meanWidth   = NaN;
            medianWidth = NaN;
        else
            meanWidth   = mean(res.widths_ms, 'omitnan');
            medianWidth = median(res.widths_ms, 'omitnan');
        end

        summaryRows(ch).dataset            = datasetName;
        summaryRows(ch).blankedMatFile     = runPath;
        summaryRows(ch).outDir             = outDir;
        summaryRows(ch).channel            = res.channel;
        summaryRows(ch).label              = res.label;
        summaryRows(ch).nSpikes            = res.nSpikes;
        summaryRows(ch).meanRate_spk_per_s = res.meanRate_spk_per_s;
        summaryRows(ch).Vpp_uv             = res.Vpp_uv;
        summaryRows(ch).SNR                = res.SNR;
        summaryRows(ch).refracViolFrac     = res.refracViolFrac;
        summaryRows(ch).meanWidth_ms       = meanWidth;
        summaryRows(ch).medianWidth_ms     = medianWidth;
        summaryRows(ch).validDuration_s    = validDuration_s;
    end

    close all force; %#ok<CLALL>
end


function row = emptySummaryRow()
    row = struct( ...
        'dataset', '', 'blankedMatFile', '', 'outDir', '', ...
        'channel', NaN, 'label', '', 'nSpikes', NaN, ...
        'meanRate_spk_per_s', NaN, 'Vpp_uv', NaN, 'SNR', NaN, ...
        'refracViolFrac', NaN, 'meanWidth_ms', NaN, 'medianWidth_ms', NaN, ...
        'validDuration_s', NaN);
end


% ==========================================================================
% ============================ Splitting helpers ===========================
% ==========================================================================
function recPath = split_and_save_recovery(fpath, sepTimeSec)
% Mirrors batch_process.m: load source, slice at sepTimeSec, write
% *_stim_blankmotion.mat and *_recovery_blankmotion.mat alongside, return
% the recovery file path.

    [folderPath, base, ~] = fileparts(fpath);

    S = load(fpath);
    [signal, fs, blankIdx] = extract_signal(S);
    if isempty(signal)
        error('Could not find a signal matrix in %s', base);
    end
    N = size(signal, 1);

    cutSamp = max(1, round(sepTimeSec * fs) + 1);
    if cutSamp >= N
        error('sepTimeSec=%g s exceeds file duration (%.1f s) for %s', ...
            sepTimeSec, N / fs, base);
    end

    stem = regexprep(base, ...
        '_(notched|blankmotion|recovery_blankmotion|stim_blankmotion)$', '');

    stimSignal   = signal(1 : cutSamp - 1, :);
    stimBlankIdx = trim_blank(blankIdx, 1, cutSamp - 1, 0);

    recSignal    = signal(cutSamp : N, :);
    recBlankIdx  = trim_blank(blankIdx, cutSamp, N, cutSamp - 1);

    stimPath = fullfile(folderPath, [stem '_stim_blankmotion.mat']);
    recPath  = fullfile(folderPath, [stem '_recovery_blankmotion.mat']);

    splitInfo = struct('sepTimeSec', sepTimeSec, 'sourceFile', base, ...
        'cutSamp', cutSamp); %#ok<NASGU>

    save_blankmotion(stimPath, stimSignal, fs, stimBlankIdx, S, 'stim',     splitInfo);
    save_blankmotion(recPath,  recSignal,  fs, recBlankIdx,  S, 'recovery', splitInfo);
end


function [signal, fs, blankIdx] = extract_signal(S)
    signal = []; fs = []; blankIdx = zeros(0, 2);
    if isfield(S, 'fs');     fs = S.fs;
    elseif isfield(S, 'Fs'); fs = S.Fs;
    end
    if isfield(S, 'yOut');       signal = double(S.yOut);
    elseif isfield(S, 'y');      signal = double(S.y);
    elseif isfield(S, 'signal'); signal = double(S.signal);
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


function bk = trim_blank(blankIdx, lo, hi, offset)
    if isempty(blankIdx); bk = zeros(0, 2); return; end
    bk = blankIdx;
    bk(bk(:, 2) < lo, :) = [];
    bk(bk(:, 1) > hi, :) = [];
    bk(:, 1) = max(bk(:, 1), lo) - offset;
    bk(:, 2) = min(bk(:, 2), hi) - offset;
end


function save_blankmotion(outPath, yOut, fs, blankIdx, srcS, region, splitInfo) %#ok<INUSD,INUSL>
    N = size(yOut, 1);
    t = (0 : N - 1)' / fs;                                            %#ok<NASGU>
    removedSegmentIdx = blankIdx;                                     %#ok<NASGU>
    blankingApplied   = ~isempty(blankIdx);                           %#ok<NASGU>
    removedSegments   = [];                                           %#ok<NASGU>

    if isfield(srcS, 'hrChanIdx') && ~isempty(srcS.hrChanIdx)
        hrChanIdx = srcS.hrChanIdx; %#ok<NASGU>
    else
        hrChanIdx = NaN; %#ok<NASGU>
    end

    save(outPath, 'yOut', 'fs', 't', 'removedSegmentIdx', 'blankingApplied', ...
        'removedSegments', 'hrChanIdx', 'splitInfo', 'region');
end
