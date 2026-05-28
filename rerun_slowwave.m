%% Re-run only slow wave analysis on the E10 dataset after gate fix.

clc; clear; close all;

folderPath = fullfile(pwd, 'E10_FRE_E10_stim_rec_0030');
condition  = "E10_FRE_E10_stim_rec_0030";

% Use only the recovery blankmotion file (single-signal API).
recFile  = fullfile(folderPath, condition + "_recovery_blankmotion.mat");
Sr = load(recFile);

fs       = Sr.fs;
rec_y    = Sr.yOut;
rec_idx  = getBlankIdxFromFile(Sr, size(rec_y,1));
swData   = rec_y(:, 3:5);  % stomach channels only

lowPassOn     = true;
lowPassOrder  = 2;       % settling time = 2/0.15 = 13.3 s
lowPassCutoff = 0.15;    % Hz — passes 3-6 cpm slow waves, cuts EMG bursts
edgeBufferSec = 15;      % covers the 13.3 s settling time of order=2 filter
window        = 5;       % s gaussian smooth — suppresses transition-band ripple
figtrue       = true;

condOut = char(condition + "_recovery");

fprintf('Re-running slow wave on recovery (lowpass ON, cutoff=%g Hz, edgeBuf=%gs)...\n', ...
    lowPassCutoff, edgeBufferSec);
swResults = slowWaveAnalysis_new(swData, lowPassOn, lowPassCutoff, ...
    lowPassOrder, fs, window, figtrue, folderPath, condOut, ...
    rec_idx, edgeBufferSec);

save(fullfile(folderPath, condition + "_slowWave_results.mat"), 'swResults');
fprintf('Done.\n');


%% ---- helper (same as in run_continuous.m) ----
function idx = getBlankIdxFromFile(S, N)
    idx = zeros(0,2);
    if isfield(S, 'removedSegmentIdx');        idx = S.removedSegmentIdx; end
    if isfield(S, 'removedSegmentIdx_motion'); idx = [idx; S.removedSegmentIdx_motion]; end
    if isfield(S, 'removedSegmentIdx_stim');   idx = [idx; S.removedSegmentIdx_stim]; end
    idx = normalizeIdxLocal(idx, N);
    idx = mergeSegmentsLocal(idx);
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
