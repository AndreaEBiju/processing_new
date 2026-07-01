function build_mmc_cache(opts)
% BUILD_MMC_CACHE  Batch driver: for every queued _blankmotion neural file, find
% its paired _HRBR file (reusing the bulk conventions + scan), extract the
% gastric motor "MMC" activity, and write <base>_mmc.mat next to each file.
%
%   build_mmc_cache(opts)
%
% REQUIRED opts (study-specific):
%   .gastricCols  3 column indices of the stomach channels within the combined
%                 (nerve+stomach) array in the _blankmotion file, e.g. [5 6 7]
%   cardiac source -- one of:  .rpeakVar (precomputed R-peaks in _HRBR, preferred)
%                              or .ecgVar (raw ECG trace in _HRBR)
% OPTIONAL: .dataVar, .rpeakUnits ('seconds'|'samples'), .rpeakFs (see extract_mmc)
% OPTIONAL: .numWorkers  parallel pool size for the file loop ([] = default pool,
%   N = N workers, 0 = serial). Files run in parallel; each worker holds one file's
%   arrays, so lower this if RAM is tight.
% OPTIONAL opts: any field accepted by extract_mmc (band, W, S, k, sigmaWin,
%   spikeRefractory, burstRefractory, cardiacBlankMs, delayW/delayStep/delayMaxLag,
%   fsVar, ...). Reasonable defaults are used otherwise.
%
% Uses the same conventions/queue UI as run_pipeline_bulk, so the suffixes and
% the neural<->HRBR pairing are identical to the rest of the pipeline. Files
% already carrying a _mmc.mat are reprocessed (delete to force a clean rebuild).

    if nargin < 1; opts = struct(); end
    assert(isfield(opts,'gastricCols') && ~isempty(opts.gastricCols), ...
        'build_mmc_cache: set opts.gastricCols (3 stomach channel indices, e.g. [5 6 7]).');
    assert((isfield(opts,'rpeakVar') && ~isempty(opts.rpeakVar)) || ...
           (isfield(opts,'ecgVar')   && ~isempty(opts.ecgVar)), ...
        'build_mmc_cache: set opts.rpeakVar (precomputed R-peaks) or opts.ecgVar (raw ECG) in the _HRBR file.');
    for fn = {'bulk_conventions_ui','bulk_queue_ui','extract_mmc'}
        assert(exist(fn{1},'file')==2, 'build_mmc_cache: %s not on path.', fn{1});
    end

    cf = fullfile(pwd,'mmc_cache_conv.mat');     % remembers conventions between runs
    conv = bulk_conventions_ui(cf);
    if isempty(conv); fprintf('Cancelled at conventions.\n'); return; end
    files = bulk_queue_ui(conv, cf);             % Add folder -> recursive scan -> Continue
    if isempty(files); fprintf('No files queued.\n'); return; end

    % ---- parallel over files (they are fully independent; each writes its own
    %      _mmc.mat). Concurrency is memory-bound -- each worker loads one file's
    %      large arrays -- so cap workers with opts.numWorkers if RAM is tight.
    %      opts.numWorkers: [] -> default pool | N -> N workers | 0 -> serial.
    addpath(fileparts(mfilename('fullpath')));   % ensure workers inherit extract_mmc on their path
    M = Inf;
    if isfield(opts,'numWorkers') && ~isempty(opts.numWorkers); M = opts.numWorkers; end
    if M ~= 0 && exist('parpool','file')==2 && isempty(gcp('nocreate'))
        try
            if isfinite(M) && M > 0; parpool(M); else; parpool; end
        catch ME
            warning('mmc:nopool','could not start parallel pool (%s); continuing.', ME.message);
        end
    end

    n = numel(files); ok = 0;
    parfor (f = 1:n, M)
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
