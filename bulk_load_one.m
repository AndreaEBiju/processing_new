function D = bulk_load_one(neuralPath, rpeakPath, opts)
% BULK_LOAD_ONE  Headless loader for one VENG dataset (bulk mode).
%
%   D = bulk_load_one(neuralPath, rpeakPath, opts)
%
% Programmatic equivalent of step0_load_data's doLoad (no GUI), so the same
% D struct can be fed to process_dataset in a batch loop. Mirrors the file
% schema written by the upstream HR/blankmotion pipeline:
%
%   neuralPath : *_blankmotion.mat  with variables
%                  yOut  [N x nCol] (samples x channels; RVN,LVN,ANT1..3)
%                  fs    scalar sample rate (Hz)
%                  t     [N x 1] time vector            (optional)
%                  removedSegmentIdx [m x 2]            (optional)
%   rpeakPath  : *_HRBR.mat with the R-peak SAMPLE INDICES (default var
%                'heartlocs'); used for cardiac blanking. May be '' (then no
%                R-peaks -> step1a blanks nothing, with a warning downstream).
%
%   opts (struct, all optional; sensible defaults):
%     .yVar          neural matrix variable name        (default 'yOut')
%     .fsVar         sample-rate variable name           (default 'fs')
%     .tVar          time variable name or ''            (default 't')
%     .segVar        removed-segment var name or ''      (default 'removedSegmentIdx')
%     .rVar          R-peak variable name in rpeakPath   (default 'heartlocs')
%     .neuralCols    columns of yOut to treat as neural  (default [1 2])
%     .labels        cellstr labels for those columns    (default {'Ch1','Ch2'} ...)
%     .versionX      informational version number        (default NaN)
%
% Returns a struct D with exactly the fields the step functions consume
% (see step0_load_data), or errors if a required variable is missing.

    if nargin < 2; rpeakPath = ''; end
    if nargin < 3 || isempty(opts); opts = struct(); end

    % ---- defaults ----
    def = struct('yVar','yOut','fsVar','fs','tVar','t', ...
                 'segVar','removedSegmentIdx','rVar','heartlocs', ...
                 'neuralCols',[1 2],'labels',{{}},'versionX',NaN);
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
            opts.(f{i}) = def.(f{i});
        end
    end
    if ~isfile(neuralPath); error('bulk_load_one:noNeural','Neural file not found: %s', neuralPath); end

    % ---- load only the variables we need from the neural file ----
    want = {opts.yVar, opts.fsVar};
    if ~isempty(opts.tVar);   want{end+1} = opts.tVar;   end
    if ~isempty(opts.segVar); want{end+1} = opts.segVar; end
    Sn = load(neuralPath, want{:});
    if ~isfield(Sn, opts.yVar);  error('bulk_load_one:noY','Variable "%s" not in %s', opts.yVar, neuralPath); end
    if ~isfield(Sn, opts.fsVar); error('bulk_load_one:noFs','Variable "%s" not in %s', opts.fsVar, neuralPath); end

    y  = double(Sn.(opts.yVar));
    fs = double(Sn.(opts.fsVar)); fs = fs(1);

    % orient to samples x channels (more rows than cols)
    if size(y,1) < size(y,2); y = y.'; end
    N = size(y,1);

    % time vector
    if ~isempty(opts.tVar) && isfield(Sn, opts.tVar) && ~isempty(Sn.(opts.tVar))
        t = Sn.(opts.tVar); t = t(:);
    else
        t = (0:N-1)' / fs;
    end

    % removed-segment indices
    if ~isempty(opts.segVar) && isfield(Sn, opts.segVar) && ~isempty(Sn.(opts.segVar))
        removedSegmentIdx = Sn.(opts.segVar);
    else
        removedSegmentIdx = zeros(0,2);
    end

    % ---- neural channel selection ----
    neuralCols = opts.neuralCols(:)';
    if max(neuralCols) > size(y,2)
        error('bulk_load_one:cols','neuralCols [%s] exceed matrix columns (%d) in %s', ...
            num2str(neuralCols), size(y,2), neuralPath);
    end
    labels = opts.labels;
    if isempty(labels)
        labels = arrayfun(@(c) sprintf('Ch%d', c), 1:numel(neuralCols), 'UniformOutput', false);
    end
    if numel(labels) ~= numel(neuralCols)
        error('bulk_load_one:labels','Need one label per neural column.');
    end

    % ---- R-peaks (sample indices) ----
    rpeakSamples = [];
    if ~isempty(rpeakPath) && isfile(rpeakPath)
        Sr = load(rpeakPath, opts.rVar);
        if isfield(Sr, opts.rVar) && ~isempty(Sr.(opts.rVar))
            rpeakSamples = round(double(Sr.(opts.rVar)(:)));
            rpeakSamples = rpeakSamples(rpeakSamples >= 1 & rpeakSamples <= N);
        else
            warning('bulk_load_one:noR','R-peak var "%s" not in %s; cardiac blanking will be empty.', ...
                opts.rVar, rpeakPath);
        end
    else
        warning('bulk_load_one:noRfile','No R-peak file for %s; cardiac blanking will be empty.', neuralPath);
    end
    rpeakTimes = (rpeakSamples - 1) / fs;

    % ---- assemble D (same shape as step0_load_data's Dout) ----
    D = struct();
    D.neuralFile        = neuralPath;
    D.rpeakFile         = rpeakPath;
    D.slowWaveFile      = '';
    D.versionX          = opts.versionX;
    D.fs                = fs;
    D.y                 = y;                  % samples x channels (raw)
    D.t                 = t;
    D.removedSegmentIdx = removedSegmentIdx;
    D.neuralChannels    = neuralCols;         % column indices into y
    D.channelLabels     = labels;
    D.rpeakSamples      = rpeakSamples;
    D.rpeakTimes        = rpeakTimes;
    D.slowWave          = [];
    D.meta              = struct('loadedAt', datetime('now'), ...
                                 'neuralVar', opts.yVar, 'fsVar', opts.fsVar, ...
                                 'mode', 'bulk');
end
