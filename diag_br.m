%% Diagnose empty BR stretches in 200-400s.

folderPath = fullfile(pwd, 'E10_FRE_E10_stim_rec_0030');
condition  = "E10_FRE_E10_stim_rec_0030";

S = load(fullfile(folderPath, condition + "_recovery_HRBR.mat"));

br = S.breathRateSeries;
mt = S.metrics_t;
fs = 24414.0625;

inWin = mt >= 200 & mt <= 400;
fprintf('200-400s range: %d metric_t points\n', sum(inWin));
fprintf('  BR NaN: %d (%.1f%%)\n', sum(isnan(br(inWin))), 100*mean(isnan(br(inWin))));
fprintf('  BR finite: %d (mean = %.2f bpm)\n', sum(~isnan(br(inWin))), mean(br(inWin),'omitnan'));

% --- Inspect breath peak times in 200-400s ---
brSamp = S.br_locs_true;
brT    = brSamp / fs;
sel    = brT >= 200 & brT <= 400;
brSel  = brT(sel);
fprintf('\nBreath troughs in 200-400s: %d\n', numel(brSel));
if numel(brSel) > 1
    isiSec  = diff(brSel);
    isiRate = 60 ./ isiSec;
    fprintf('  ISI range: %.2f - %.2f s\n', min(isiSec), max(isiSec));
    fprintf('  ISI rate range: %.1f - %.1f bpm\n', min(isiRate), max(isiRate));
    fprintf('  intervals < 40 bpm: %d (%.1f%%)\n', ...
        sum(isiRate < 40), 100*mean(isiRate < 40));
    fprintf('  intervals > 170 bpm: %d (%.1f%%)\n', ...
        sum(isiRate > 170), 100*mean(isiRate > 170));
end

% --- Sample 10 NaN BR windows in 200-400s and check why ---
nanIdx = find(inWin & isnan(br));
sampleIdx = nanIdx(round(linspace(1, numel(nanIdx), min(10, numel(nanIdx)))));
fprintf('\n--- Sample of NaN BR windows in 200-400s ---\n');
fprintf('Constants: brMinStretchSec=%g, brMinPeaksInStretch=%d, minBreathRate=%g, maxBreathRate=%g, hrBrWinSec=%g\n', ...
    S.brMinStretchSec, S.brMinPeaksInStretch, S.minBreathRate_bpm, S.maxBreathRate_bpm, S.hrBrWinSec);

invalid = S.invalidMask;
edge    = S.edgeMask;
clean   = ~invalid & ~edge;
N       = numel(invalid);
hrBrWinSec = S.hrBrWinSec;
halfW   = hrBrWinSec / 2;

for k = 1:numel(sampleIdx)
    tc   = mt(sampleIdx(k));
    if tc < halfW || tc > (N-1)/fs - halfW
        fprintf('t=%.1fs: edge (window does not fit)\n', tc);
        continue;
    end
    idxCtr = max(1, min(N, round(tc*fs)+1));
    if invalid(idxCtr)
        fprintf('t=%.1fs: center inside invalidMask\n', tc);
        continue;
    end
    idx0 = max(1, round((tc-halfW)*fs)+1);
    idx1 = min(N, round((tc+halfW)*fs)+1);
    winClean = clean(idx0:idx1);
    [stretchLen, stretchStart] = longestContiguousRun(winClean);
    stretchSec = stretchLen / fs;

    absStart = idx0 + stretchStart - 1;
    absEnd   = min(N, absStart + stretchLen - 1);
    peaksInStretch = brSamp(brSamp >= absStart & brSamp <= absEnd);
    nP = numel(peaksInStretch);

    if stretchSec < S.brMinStretchSec
        fprintf('t=%.1fs: longest clean stretch %.1fs < %g (likely small blanks fragmenting window)\n', ...
            tc, stretchSec, S.brMinStretchSec);
        continue;
    end
    if nP < S.brMinPeaksInStretch
        fprintf('t=%.1fs: only %d peaks in stretch (need %d)\n', tc, nP, S.brMinPeaksInStretch);
        continue;
    end
    isiSec  = diff(peaksInStretch) / fs;
    isiRate = 60 ./ isiSec;
    nOOR    = sum(isiRate < S.minBreathRate_bpm | isiRate > S.maxBreathRate_bpm);
    if nOOR > 0
        fprintf('t=%.1fs: %d peaks, %d/%d intervals out of [%g %g] bpm (range %.1f-%.1f)\n', ...
            tc, nP, nOOR, numel(isiRate), S.minBreathRate_bpm, S.maxBreathRate_bpm, ...
            min(isiRate), max(isiRate));
    else
        fprintf('t=%.1fs: stretch ok, %d peaks, all ISI in range — should NOT be NaN\n', tc, nP);
    end
end

%% --- helper: longest contiguous true run ---
function [longestLen, longestStart] = longestContiguousRun(mask)
    if isempty(mask); longestLen = 0; longestStart = 1; return; end
    mask = mask(:)';
    d    = diff([false, mask, false]);
    s    = find(d == 1);
    e    = find(d == -1) - 1;
    lens = e - s + 1;
    [longestLen, j] = max([0, lens]);
    if longestLen == 0; longestStart = 1; else; longestStart = s(j-1); end
end
