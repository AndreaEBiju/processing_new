function hr = resolve_hrbr(neuralPath, phase, conv)
% RESOLVE_HRBR  Find the heartbeat (HRBR) file for a neural file, robustly.
%
%   hr = resolve_hrbr(neuralPath, phase, conv)
%
% Tries, in order:
%   1. the exact constructed name  <stem><HRBRsuffix>.mat
%   2. a prefix glob               <stem>*HRBR*.mat
%   3. folder-wide match: among *HRBR*.mat in the same folder, the one of the
%      right phase (recovery HRBR contains 'recovery'; baseline does not) with
%      the longest shared leading prefix with the neural name. This catches the
%      case where the HRBR name is SHORTER than the neural (e.g. it drops the
%      'stim_rec_1406_notched' middle and is just <animal_condition>_..._recovery_HRBR).
% Returns '' if nothing plausible is found.

    [folder, base] = fileparts(neuralPath);
    isRec = strcmpi(phase, 'recovery');
    if isRec; sufN = conv.sufRecNeural; sufH = conv.sufRecHRBR;
    else;     sufN = conv.sufBaseNeural; sufH = conv.sufBaseHRBR; end

    if ~isempty(sufN) && endsWith(base, sufN); stem = base(1:end-numel(sufN)); else; stem = base; end

    % 1) exact
    hr = fullfile(folder, [stem sufH '.mat']);
    if exist(hr,'file'); return; end

    % 2) prefix glob (HRBR starts with the full stem)
    alt = dir(fullfile(folder, [stem '*HRBR*.mat']));
    alt = alt(~[alt.isdir]);
    if ~isempty(alt); hr = fullfile(alt(1).folder, alt(1).name); return; end

    % 3) folder-wide HRBR, phase-filtered, longest shared prefix with `base`
    cand = dir(fullfile(folder, '*HRBR*.mat'));
    cand = cand(~[cand.isdir]);
    names = {cand.name};
    isRecName = contains(names, 'recovery', 'IgnoreCase', true);
    keep = isRec == isRecName;                  % recovery<->recovery, baseline<->non-recovery
    cand = cand(keep);
    if isempty(cand); hr = ''; return; end

    best = ''; bestL = 0;
    for i = 1:numel(cand)
        L = lcp_len(base, cand(i).name);
        if L > bestL; bestL = L; best = fullfile(cand(i).folder, cand(i).name); end
    end
    if bestL >= 6; hr = best; else; hr = ''; end   % require a meaningful shared prefix
end

function L = lcp_len(a, b)
    n = min(numel(a), numel(b)); L = 0;
    while L < n && a(L+1) == b(L+1); L = L + 1; end
end
