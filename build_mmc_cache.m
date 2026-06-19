function build_mmc_cache(opts)
% BUILD_MMC_CACHE  Batch driver: for every queued _blankmotion neural file, find
% its paired _HRBR file (reusing the bulk conventions + scan), extract the
% gastric motor "MMC" activity, and write <base>_mmc.mat next to each file.
%
%   build_mmc_cache(opts)
%
% REQUIRED opts (the variable names are study-specific):
%   .gastricVar   name of the N x 3 gastric array in the _blankmotion file
%   .ecgVar       name of the ECG trace in the _HRBR file
% OPTIONAL opts: any field accepted by extract_mmc (band, W, S, k, sigmaWin,
%   refractory, cardiacBlankMs, storeFs, delayW/delayStep/delayMaxLag, fsVar,
%   ecgFsVar, ...). Reasonable defaults are used otherwise.
%
% Uses the same conventions/queue UI as run_pipeline_bulk, so the suffixes and
% the neural<->HRBR pairing are identical to the rest of the pipeline. Files
% already carrying a _mmc.mat are reprocessed (delete to force a clean rebuild).

    if nargin < 1; opts = struct(); end
    assert(isfield(opts,'gastricVar') && ~isempty(opts.gastricVar), ...
        'build_mmc_cache: set opts.gastricVar (the N x 3 gastric variable name).');
    assert(isfield(opts,'ecgVar') && ~isempty(opts.ecgVar), ...
        'build_mmc_cache: set opts.ecgVar (the ECG variable name in the _HRBR file).');
    for fn = {'bulk_conventions_ui','bulk_queue_ui','extract_mmc'}
        assert(exist(fn{1},'file')==2, 'build_mmc_cache: %s not on path.', fn{1});
    end

    cf = fullfile(pwd,'mmc_cache_conv.mat');     % remembers conventions between runs
    conv = bulk_conventions_ui(cf);
    if isempty(conv); fprintf('Cancelled at conventions.\n'); return; end
    files = bulk_queue_ui(conv, cf);             % Add folder -> recursive scan -> Continue
    if isempty(files); fprintf('No files queued.\n'); return; end

    n = numel(files); ok = 0;
    for f = 1:n
        F = files(f);
        if isempty(F.hrbr) || exist(F.hrbr,'file')~=2
            warning('mmc:nohrbr','%s: no paired _HRBR file -- skipped.', F.stem); continue;
        end
        try
            fprintf('[mmc %d/%d] %s  (%s / %s)\n', f, n, F.stem, F.condition, F.phase);
            extract_mmc(F.neural, F.hrbr, opts);
            ok = ok + 1;
        catch ME
            warning('mmc:fail','  skipped %s (%s)', F.stem, ME.message);
        end
    end
    fprintf('[mmc] wrote %d / %d _mmc.mat files.\n', ok, n);
end
