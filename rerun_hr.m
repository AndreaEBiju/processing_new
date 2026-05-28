%% Re-run only HR/BR/HRV on E10 data after adding heartCountSeries.

clc; clear; close all;

folderPath = fullfile(pwd, 'E10_FRE_E10_stim_rec_0030');
condition  = "E10_FRE_E10_stim_rec_0030";

% Use only the recovery blankmotion file (single-signal API).
recFile  = fullfile(folderPath, condition + "_recovery_blankmotion.mat");
Sr = load(recFile);

fs        = Sr.fs;
rec_y     = Sr.yOut;
rec_idx   = getBlankIdxFromFile(Sr, size(rec_y,1));
hrChanIdx = getHRChanIdxFromFile(Sr);

cutoff       = 8;
order        = 4;
winSec       = 20;     % HRV / per-window beat count window
hrBrWinSec   = 60;     % HR / BR window
stepSec      = 1;
edgeBufferSec = 0.75;
figtrue      = true;

condOut = char(condition + "_recovery");

fprintf('Re-running HR/BR/HRV on recovery (winSec=%ds, hrBrWinSec=%ds, centered)...\n', ...
    winSec, hrBrWinSec);
results = HR_BR_HRVAnalysis_new(rec_y, fs, cutoff, order, ...
    folderPath, condOut, hrChanIdx, rec_idx, ...
    edgeBufferSec, winSec, stepSec, figtrue, hrBrWinSec);

save(fullfile(folderPath, condition + "_HR_BR_HRV_results.mat"), 'results');
fprintf('Done.\n');


%% ---- helpers ----
function idx = getBlankIdxFromFile(S, N)
    idx = zeros(0,2);
    if isfield(S, 'removedSegmentIdx');        idx = S.removedSegmentIdx; end
    if isfield(S, 'removedSegmentIdx_motion'); idx = [idx; S.removedSegmentIdx_motion]; end
    if isfield(S, 'removedSegmentIdx_stim');   idx = [idx; S.removedSegmentIdx_stim]; end
    idx = normalizeIdxLocal(idx, N);
    idx = mergeSegmentsLocal(idx);
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
