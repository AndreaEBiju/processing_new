function results = run_dfa_batch(folders, opts)
% RUN_DFA_BATCH  Run gap-aware DFA -- Step 7 (nerve) + dfaRR_gapAware
% (cardiac) -- across every file discovered under `folders`, reusing this
% repo's existing bulk-scan/cache conventions wherever one already exists.
%
%   results = run_dfa_batch(folders, opts)
%
% folders : cell array (or single char/string) of root folders to scan.
%           Same convention as run_pipeline_bulk.m / bulk_scan_files.m --
%           searched recursively.
%
% opts (all optional):
%   .cacheFile      path to an EXISTING run_pipeline_bulk.m cache (default:
%                   bulk_cache.mat in pwd). Read to reuse its saved
%                   file-naming .conv conventions (baseline/recovery
%                   suffixes, RVN/LVN column, heartbeat var name) so you
%                   don't retype them. If that cache has no 'conv' yet,
%                   bulk_conventions_ui prompts ONCE (a UI dialog -- this
%                   is the one part of this script that isn't headless,
%                   same as every other bulk_* tool here) and persists the
%                   answer into the SAME file for next time.
%   .forceRefresh   ignore this script's own per-file DFA caches and
%                   recompute everything (default false)
%   .saveFigs       save Step 7's diagnostic PNGs per channel (default true)
%   .figsDir        override where Step 7 figures land (default: a
%                   'dfa_figs' folder next to each neural file)
%   .P              pipeline_params() struct (default: pipeline_params())
%   .lowRateFlagPct percentile (of this batch's firing rates) used to flag
%                   alphaFull~0.5 + low-rate channels as possible
%                   noise-floor artifacts rather than real findings
%                   (default 25)
%
% Returns results.nerve (one row per channel per neural file) and
% results.cardiac (one row per _HRVMeasures.mat file found), plus prints a
% summary to the console.
%
% ---------------------------------------------------------------------
% This covers, automatically, across every file, the three follow-up
% checks that came up after the Step 7 build:
%
%   1. SAVE-PATH VERIFICATION. After Step 7 computes D.dfa for a file, this
%      calls pipeline_save_summary and RELOADS the saved .mat to confirm
%      S.dfa_alpha1(k)/dfa_alpha2(k) on disk match D.dfa(k) in memory --
%      not just trusting the in-memory struct. Recorded per-channel as
%      .summaryMatch; the summary should report 0 mismatches.
%
%   2. THE "Ch2-STYLE" ANOMALY CHECK. Any channel whose alphaFull lands
%      within 0.05 of 0.5 (near-white-noise scaling) AND whose mean firing
%      rate falls in the bottom opts.lowRateFlagPct percentile of this
%      whole batch gets .lowRateSuspect = true. This does NOT decide
%      whether it's a real finding or a sparse-rate artifact -- per the
%      plan's own framing, that's a human judgment call -- it just
%      surfaces the candidates instead of leaving you to eyeball every
%      channel's numbers by hand.
%
%   3. THE CARDIAC REAL-DATA CHECK. Every _HRVMeasures.mat found gets its
%      dfa_alpha1/dfa_alpha2 checked against the ~0.5-1.3 plausible range
%      for rat HRV (per DFA_IMPLEMENTATION_PLAN.md T15); anything outside
%      it is flagged rather than silently accepted.
%
% ---------------------------------------------------------------------
% SCOPE NOTES (read before trusting "all files" to mean everything):
%
% - CARDIAC: if a _HRVMeasures.mat predates the DFA wiring (has
%   RR_intervals/RR_times but no dfa_alpha1 yet), this script computes DFA
%   directly from the cached RR series and re-saves it into that same
%   file -- it does NOT re-run HR_BR_HRVAnalysis_new.m's full heartbeat
%   detection from the raw signal (that needs the original fs/cutoff/
%   order/channel index, which this script has no reliable way to infer
%   for an arbitrary file). If a file has NEITHER dfa_alpha1 NOR
%   RR_intervals/RR_times, it's reported as 'needs-full-rerun' and
%   skipped -- run HR_BR_HRVAnalysis_new.m yourself for those.
%
% - NERVE: confirmed (via direct inspection) that no existing cache in
%   this repo persists Step 6's fr_validFrac/D.validMask -- the existing
%   <base>_vengmetrics.mat sidecar (written by run_pipeline_bulk.m) only
%   stores fr_hz/fr_t, not the valid-fraction Step 7 needs. So Step 7
%   cannot be reconstructed from that cache alone; this script re-runs
%   Steps 1-6 (process_dataset) the FIRST time it sees a given neural
%   file, then writes its OWN sidecar cache (<base>_dfa_nerve.mat,
%   deliberately separate from run_pipeline_bulk.m's cache/schema so
%   nothing here can corrupt that other tool's data) so a repeat run of
%   THIS script skips the expensive part.

    if nargin < 1 || isempty(folders)
        error('run_dfa_batch:noFolders', 'Pass a cell array of root folders to scan.');
    end
    if ischar(folders) || isstring(folders); folders = {char(folders)}; end
    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts,'cacheFile');      opts.cacheFile = fullfile(pwd,'bulk_cache.mat'); end
    if ~isfield(opts,'forceRefresh');   opts.forceRefresh = false; end
    if ~isfield(opts,'saveFigs');       opts.saveFigs = true; end
    if ~isfield(opts,'figsDir');        opts.figsDir = ''; end
    if ~isfield(opts,'P') || isempty(opts.P); opts.P = pipeline_params(); end
    if ~isfield(opts,'lowRateFlagPct'); opts.lowRateFlagPct = 25; end
    P = opts.P;

    fprintf('=== run_dfa_batch: nerve (Step 7) ===\n');
    nerveRows = run_nerve_dfa(folders, opts, P);

    fprintf('\n=== run_dfa_batch: cardiac ===\n');
    cardiacRows = run_cardiac_dfa(folders, opts, P);

    results = struct('nerve', nerveRows, 'cardiac', cardiacRows);
    print_summary(results);
end

% ========================================================================
% NERVE (Step 7)
% ========================================================================
function rows = run_nerve_dfa(folders, opts, P)
    cf = opts.cacheFile;
    conv = bulk_cache_get(cf, 'conv');
    if ~isstruct(conv) || isempty(fieldnames(conv))
        conv = bulk_conventions_ui(cf);
        if isempty(conv)
            warning('run_dfa_batch:noConv', 'Conventions cancelled; skipping nerve DFA.');
            rows = repmat(empty_nerve_row(), 0, 1);
            return;
        end
    end

    files = bulk_scan_files(folders, conv);
    loadOpts = struct('neuralCols',[conv.rvnCol conv.lvnCol], 'labels',{{'RVN','LVN'}}, 'rVar',conv.heartVar);

    rows = repmat(empty_nerve_row(), 0, 1);
    SCHEMA = dfa_nerve_schema();
    for f = 1:numel(files)
        F = files(f);
        fprintf('[nerve %d/%d] %s (%s / %s)\n', f, numel(files), F.stem, F.condition, F.phase);
        [nfDir, nfBase] = fileparts(F.neural);
        cachePath = fullfile(nfDir, [nfBase '_dfa_nerve.mat']);
        mt = file_mtime_local(F.neural);

        fileRows = [];
        if ~opts.forceRefresh
            fileRows = try_load_nerve_cache(cachePath, mt, SCHEMA);
        end
        if isempty(fileRows)
            try
                fileRows = compute_nerve_dfa(F, P, loadOpts, opts);
                rows_ = fileRows; srcMtime = mt; schema = SCHEMA; %#ok<NASGU>
                save(cachePath, 'rows_', 'srcMtime', 'schema');
            catch ME
                warning('run_dfa_batch:nerveFail', '  skipped %s (%s)', F.stem, ME.message);
                continue;
            end
        else
            fprintf('  (cached)\n');
        end
        rows = [rows; fileRows]; %#ok<AGROW>
    end

    % Anomaly flag #2 needs the whole batch's rate distribution first, so it
    % runs as a pass over the assembled rows, not per-file.
    if ~isempty(rows)
        rates = [rows.meanRate_hz];
        finiteRates = rates(isfinite(rates));
        if ~isempty(finiteRates)
            thr = prctile(finiteRates, opts.lowRateFlagPct);
            for i = 1:numel(rows)
                nearHalf = isfinite(rows(i).alphaFull) && abs(rows(i).alphaFull - 0.5) < 0.05;
                lowRate  = isfinite(rows(i).meanRate_hz) && rows(i).meanRate_hz <= thr;
                rows(i).lowRateSuspect = nearHalf && lowRate;
            end
        end
    end
end

function rows = compute_nerve_dfa(F, P, loadOpts, opts)
    prevVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    cu = onCleanup(@() set(0, 'DefaultFigureVisible', prevVis)); %#ok<NASGU>

    D = bulk_load_one(F.neural, F.hrbr, loadOpts);
    D = process_dataset(D, P, false);          % Steps 1-6, no per-file figures
    D = step7_dfa_report(D, P, opts.saveFigs);  % plotMode controls Step 7's own diagnostic figure

    if opts.saveFigs
        figDir = opts.figsDir;
        if isempty(figDir); figDir = fullfile(fileparts(F.neural), 'dfa_figs'); end
        try
            save_all_figures(figDir, [F.stem '_step7'], {'png'});
        catch ME
            warning('run_dfa_batch:figsave', '  figure save failed for %s (%s)', F.stem, ME.message);
        end
    end
    close all force; %#ok<CLALL>

    % Check #1: reload the saved summary rather than trusting D in memory.
    summaryMatch = false(1, numel(D.neuralChannels));
    try
        outPath = pipeline_save_summary(D, P, F.phase);
        S = load(outPath, 'summary');
        for k = 1:numel(D.neuralChannels)
            summaryMatch(k) = isequaln(S.summary.dfa_alpha1(k), D.dfa(k).alpha1) && ...
                               isequaln(S.summary.dfa_alpha2(k), D.dfa(k).alpha2) && ...
                               isequaln(S.summary.dfa_alphaFull(k), D.dfa(k).alphaFull);
        end
    catch ME
        warning('run_dfa_batch:summaryCheck', '  pipeline_save_summary check failed for %s (%s)', F.stem, ME.message);
    end

    rows = repmat(empty_nerve_row(), numel(D.neuralChannels), 1);
    for k = 1:numel(D.neuralChannels)
        r = empty_nerve_row();
        r.stem = F.stem; r.animal = F.animal; r.condition = F.condition; r.phase = F.phase;
        r.label = D.channelLabels{k};
        r.alpha1 = D.dfa(k).alpha1; r.alpha2 = D.dfa(k).alpha2; r.alphaFull = D.dfa(k).alphaFull;
        r.R2_1 = D.dfa(k).R2_1; r.R2_2 = D.dfa(k).R2_2; r.R2_full = D.dfa(k).R2_full;
        r.nCross = D.dfa(k).nCross; r.nExcluded = numel(D.dfa(k).excludedScales);
        r.meanRate_hz = D.metrics(k).meanRate_hz;
        r.summaryMatch = summaryMatch(k);
        r.lowRateSuspect = false; % filled in by the batch-wide pass in run_nerve_dfa
        rows(k) = r;
    end
end

function r = empty_nerve_row()
    r = struct('stem','', 'animal','', 'condition','', 'phase','', 'label','', ...
        'alpha1',NaN, 'alpha2',NaN, 'alphaFull',NaN, 'R2_1',NaN, 'R2_2',NaN, 'R2_full',NaN, ...
        'nCross',NaN, 'nExcluded',0, 'meanRate_hz',NaN, 'summaryMatch',false, 'lowRateSuspect',false);
end

function s = dfa_nerve_schema(); s = 1; end % bump if empty_nerve_row's fields change

function rows = try_load_nerve_cache(pf, mt, schema)
    rows = [];
    if ~exist(pf, 'file'); return; end
    try
        S = load(pf);
        if isfield(S,'schema') && isequal(S.schema,schema) && isfield(S,'srcMtime') ...
                && isequal(S.srcMtime,mt) && isfield(S,'rows_')
            rows = S.rows_;
        end
    catch
    end
end

function mt = file_mtime_local(p)
    di = dir(p); if isempty(di); mt = NaN; else; mt = di.datenum; end
end

% ========================================================================
% CARDIAC (dfaRR_gapAware)
% ========================================================================
function rows = run_cardiac_dfa(folders, opts, P)
    hrvHits  = scan_pattern(folders, '*_HRVMeasures.mat');
    hrbrHits = scan_pattern(folders, '*_HRBR.mat');
    nHrvBefore = numel(hrvHits); nHrbrBefore = numel(hrbrHits);
    hrvHits  = filter_canonical(hrvHits,  '_HRVMeasures.mat');
    hrbrHits = filter_canonical(hrbrHits, '_HRBR.mat');
    fprintf('[cardiac] %d _HRVMeasures.mat (%d after naming-convention filter), %d _HRBR.mat (%d after filter).\n', ...
        nHrvBefore, numel(hrvHits), nHrbrBefore, numel(hrbrHits));

    % _HRBR.mat and _HRVMeasures.mat are saved from the SAME conditionLabel +
    % folderpath by HR_BR_HRVAnalysis_new.m (just different suffix), so their
    % stems match exactly for the same trial -- group by that shared key so a
    % trial whose _HRVMeasures.mat lacks RR data can still fall back to its
    % _HRBR.mat sibling (which already has heartlocs+invalidMask from the
    % peak detection that already ran -- no raw-signal re-run needed).
    keysHrv  = arrayfun(@(h) trial_key(h.mat, '_HRVMeasures.mat'), hrvHits,  'UniformOutput', false);
    keysHrbr = arrayfun(@(h) trial_key(h.mat, '_HRBR.mat'),        hrbrHits, 'UniformOutput', false);
    allKeys = unique([keysHrv(:); keysHrbr(:)]);

    rows = repmat(empty_cardiac_row(), 0, 1);
    for t = 1:numel(allKeys)
        key = allKeys{t};
        hrvPath = ''; hrbrPath = '';
        j = find(strcmp(keysHrv, key), 1);  if ~isempty(j); hrvPath  = hrvHits(j).mat;  end
        j = find(strcmp(keysHrbr, key), 1); if ~isempty(j); hrbrPath = hrbrHits(j).mat; end

        fprintf('[cardiac %d/%d] %s\n', t, numel(allKeys), key);
        r = empty_cardiac_row();
        r.file = ternary_str(~isempty(hrvPath), hrvPath, hrbrPath);
        try
            S = struct();
            if ~isempty(hrvPath)
                S = load(hrvPath);
            elseif exist([key '_dfa_from_hrbr.mat'], 'file')
                % A prior run already derived DFA from this HRBR-only trial
                % and left it in the sidecar -- reuse it instead of
                % recomputing from scratch every time.
                S = load([key '_dfa_from_hrbr.mat']);
                r.file = [key '_dfa_from_hrbr.mat'];
            end
            if isfield(S, 'dfa_alpha1')
                r.source = 'cached';
                r.dfa_alpha1 = S.dfa_alpha1; r.dfa_alpha2 = S.dfa_alpha2; r.dfa_alphaFull = S.dfa_alphaFull;
            elseif isfield(S,'RR_intervals') && isfield(S,'RR_times') && numel(S.RR_intervals) >= 2
                dfaOut = dfaRR_gapAware(S.RR_intervals, S.RR_times, P);
                r.source = 'recomputed-from-RR';
                [r, S1] = fill_dfa_row(r, dfaOut); %#ok<NASGU>
                save(hrvPath, '-struct', 'S1', '-append');
            elseif isfield(S, 'alpha1') && isfield(S, 'alpha2')
                % Pre-existing OLD, non-gap-aware DFA (dfaRR.m's output, saved
                % under bare alpha1/alpha2/nVals/F -- from before the gap-aware
                % rewrite). No RR_intervals/RR_times cached, so it can't be
                % upgraded to dfa_alpha1 here -- but this is NOT "nothing":
                % flagged distinctly so it's never confused with or pooled
                % alongside the new gap-aware dfa_alpha1/dfa_alpha2.
                r.source = 'old-ungapped-only';
                r.dfa_alpha1 = S.alpha1; r.dfa_alpha2 = S.alpha2;
            elseif ~isempty(hrbrPath)
                H = load(hrbrPath);
                if isfield(H,'heartlocs') && isfield(H,'invalidMask') && isfield(H,'t') && numel(H.heartlocs) >= 2
                    fs = 1 / median(diff(H.t(:)));
                    [RR_intervals, RR_times] = computeValidRRIntervals(H.heartlocs, H.invalidMask, fs);
                    if numel(RR_intervals) >= 2
                        dfaOut = dfaRR_gapAware(RR_intervals, RR_times, P);
                        r.source = 'recomputed-from-HRBR';
                        [r, S1] = fill_dfa_row(r, dfaOut);
                        S1.RR_intervals = RR_intervals; S1.RR_times = RR_times; %#ok<STRNU>
                        if ~isempty(hrvPath)
                            % A real _HRVMeasures.mat exists (just missing RR/dfa
                            % data) -- augment it in place, same as the
                            % recomputed-from-RR path.
                            save(hrvPath, '-struct', 'S1', '-append');
                        else
                            % No _HRVMeasures.mat at all for this trial -- do NOT
                            % fabricate one under that name (it would misleadingly
                            % imply the full hrv/rmssd/pnn5/... suite exists when
                            % only DFA was computed here). Save a distinctly-named
                            % sidecar instead.
                            sidecarPath = [key '_dfa_from_hrbr.mat'];
                            save(sidecarPath, '-struct', 'S1');
                            r.file = sidecarPath;
                        end
                    else
                        r.source = 'needs-full-rerun';
                    end
                else
                    r.source = 'needs-full-rerun';
                end
            else
                r.source = 'needs-full-rerun';
            end
        catch ME
            warning('run_dfa_batch:cardiacFail', '  failed for %s (%s)', key, ME.message);
            r.source = 'error';
        end
        if any(strcmp(r.source, {'cached','recomputed-from-RR','old-ungapped-only','recomputed-from-HRBR'}))
            r.plausible = isfinite(r.dfa_alpha1) && r.dfa_alpha1 >= 0.5 && r.dfa_alpha1 <= 1.3 && ...
                (isnan(r.dfa_alpha2) || (r.dfa_alpha2 >= 0.5 && r.dfa_alpha2 <= 1.3));
        end
        rows(end+1) = r; %#ok<AGROW>
    end
end

function [r, S1] = fill_dfa_row(r, dfaOut)
    r.dfa_alpha1 = dfaOut.alpha1; r.dfa_alpha2 = dfaOut.alpha2; r.dfa_alphaFull = dfaOut.alphaFull;
    S1 = struct('dfa_alpha1', dfaOut.alpha1, 'dfa_alpha2', dfaOut.alpha2, 'dfa_alphaFull', dfaOut.alphaFull, ...
        'dfa_R2_1', dfaOut.R2_1, 'dfa_R2_2', dfaOut.R2_2, 'dfa_nCross', dfaOut.nCross, ...
        'dfa_nWindows', dfaOut.nWindows, 'dfa_excludedScales', dfaOut.excludedScales);
end

function s = ternary_str(c, a, b)
    if c; s = a; else; s = b; end
end

function k = trial_key(p, suf)
    k = p;
    if endsWith(k, suf); k = k(1:end-numel(suf)); end
end

function hits = scan_pattern(folders, pattern)
    hits = struct('mat', {});
    for fi = 1:numel(folders)
        d = dir(fullfile(folders{fi}, '**', pattern));
        d = d(~[d.isdir]);
        for k = 1:numel(d)
            hits(end+1).mat = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
    end
end

function r = empty_cardiac_row()
    r = struct('file','', 'source','', 'dfa_alpha1',NaN, 'dfa_alpha2',NaN, 'dfa_alphaFull',NaN, 'plausible',false);
end

function out = filter_canonical(hits, suffix)
% Prefer the canonical "<stem>_v0.2.X_[recovery_]blankmotion<suffix>"
% naming (matches the neural-file convention; suffix is '_HRVMeasures.mat'
% or '_HRBR.mat') over older/bare-named duplicates of the SAME recording
% in the same folder, per the user's confirmation that the canonical
% version is always the one to trust.
% Grouped per-FOLDER (not per-trial-stem) -- real folders were observed to
% mix "_notched_v0.2.1_blankmotion", "_notched_v0.2.2_blankmotion", and a
% bare "*_HRVMeasures.mat" all referring to the same underlying recording,
% but the "_notched" segment breaks simple stem-prefix matching. Folder-
% level grouping is the pragmatic tradeoff: if a folder ever legitimately
% holds multiple DISTINCT recordings saved under different versions, this
% could over-exclude an older-but-still-valid one -- flagged here so it's
% not a silent assumption.
    if isempty(hits); out = hits; return; end
    n = numel(hits);
    folders = cell(1, n);
    versions = nan(1, n);
    isCanon = false(1, n);
    suffixTag = regexprep(suffix, '\.mat$', ''); % e.g. '_HRVMeasures' or '_HRBR'
    pat = ['_v0\.2\.(\d+)_(?:recovery_)?blankmotion' regexptranslate('escape', suffixTag) '$'];
    for i = 1:n
        [d, b] = fileparts(hits(i).mat);
        folders{i} = d;
        tok = regexp(b, pat, 'tokens', 'once');
        if ~isempty(tok)
            isCanon(i) = true;
            versions(i) = str2double(tok{1});
        end
    end
    uFolders = unique(folders);
    keep = false(1, n);
    for fi = 1:numel(uFolders)
        idx = find(strcmp(folders, uFolders{fi}));
        canonIdx = idx(isCanon(idx));
        if ~isempty(canonIdx)
            maxV = max(versions(canonIdx));
            keep(canonIdx(versions(canonIdx) == maxV)) = true;
        else
            keep(idx) = true; % no canonical option in this folder -- use what's there
        end
    end
    out = hits(keep);
end

% ========================================================================
% SUMMARY
% ========================================================================
function print_summary(results)
    fprintf('\n================ run_dfa_batch summary ================\n');
    n = results.nerve;
    fprintf('Nerve: %d channel-file row(s).\n', numel(n));
    if ~isempty(n)
        nMismatch = sum(~[n.summaryMatch]);
        fprintf('  %d/%d summary-save mismatches (should be 0)\n', nMismatch, numel(n));
        if nMismatch > 0
            for i = find(~[n.summaryMatch])
                fprintf('    MISMATCH: %s / %s / %s (%s)\n', n(i).animal, n(i).condition, n(i).phase, n(i).label);
            end
        end
        suspectIdx = find([n.lowRateSuspect]);
        fprintf('  %d channel(s) flagged low-rate/near-0.5 alphaFull (check by hand):\n', numel(suspectIdx));
        for i = suspectIdx
            fprintf('    %s / %s / %s (%s): alphaFull=%.3f, meanRate=%.3f Hz\n', ...
                n(i).animal, n(i).condition, n(i).phase, n(i).label, n(i).alphaFull, n(i).meanRate_hz);
        end
    end

    c = results.cardiac;
    fprintf('Cardiac: %d file(s).\n', numel(c));
    if ~isempty(c)
        srcs = {c.source};
        uniqueSrcs = unique(srcs);
        for s = uniqueSrcs
            fprintf('  %d file(s): %s\n', sum(strcmp(srcs, s{1})), s{1});
        end
        badIdx = find(~[c.plausible] & (strcmp(srcs,'cached') | strcmp(srcs,'recomputed-from-RR') | strcmp(srcs,'recomputed-from-HRBR')));
        if ~isempty(badIdx)
            fprintf('  %d file(s) with IMPLAUSIBLE dfa_alpha (outside ~0.5-1.3 for rat HRV):\n', numel(badIdx));
            for i = badIdx
                fprintf('    %s: alpha1=%.3f alpha2=%.3f\n', c(i).file, c(i).dfa_alpha1, c(i).dfa_alpha2);
            end
        end
        oldIdx = find(strcmp(srcs, 'old-ungapped-only'));
        if ~isempty(oldIdx)
            fprintf(['  %d file(s) have OLD non-gap-aware DFA only (dfaRR.m-era alpha1/alpha2,\n' ...
                     '  no cached RR data to upgrade) -- do NOT pool these with dfa_alpha1 elsewhere:\n'], numel(oldIdx));
            for i = oldIdx
                fprintf('    %s: alpha1=%.3f alpha2=%.3f (plausible=%d)\n', ...
                    c(i).file, c(i).dfa_alpha1, c(i).dfa_alpha2, c(i).plausible);
            end
        end
        needRerun = find(strcmp(srcs, 'needs-full-rerun'));
        if ~isempty(needRerun)
            fprintf('  %d file(s) need a full HR_BR_HRVAnalysis_new.m re-run (no RR data cached):\n', numel(needRerun));
            for i = needRerun; fprintf('    %s\n', c(i).file); end
        end
    end
    fprintf('=========================================================\n');
end
