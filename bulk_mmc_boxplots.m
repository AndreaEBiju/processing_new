function figs = bulk_mmc_boxplots(queueFile, saveDir, opts)
% BULK_MMC_BOXPLOTS  Box + violin of the baseline-normalized gastric MMC metrics,
% grouped by condition, using the SAME styling as the nerve plots
% (bulk_plot_boxviolin -> your boxViolinPlot.m + whiten_figure).
%
%   bulk_mmc_boxplots                          % gemsplots_queue.mat in pwd, display
%   figs = bulk_mmc_boxplots(queueFile, saveDir, opts)
%
% One figure per metric (rate/burst/amp) x channel label. Each violin is that
% condition's normalized recovery distribution (per-window, pooled across
% animals), y = (window - mean_baseline)/mean_baseline, dotted line at 0.
%
% opts.groups      : cell of condition groups (defineGroupsUI form).  [queue's own]
% opts.yMax        : symmetric y cap about 0 on every panel.           [auto]
% opts.channelMode : 'mean' (default, label 'G') | 'separate' (G1,G2,G3).
% opts.metrics     : subset of {'rate','burst','amp'}.                 [all three]

    if nargin < 1 || isempty(queueFile); queueFile = fullfile(pwd,'gemsplots_queue.mat'); end
    if nargin < 2; saveDir = ''; end
    if nargin < 3; opts = struct(); end
    assert(exist('compile_mmc','file')==2, 'bulk_mmc_boxplots: compile_mmc.m not on path.');
    assert(exist('bulk_plot_boxviolin','file')==2, 'bulk_mmc_boxplots: bulk_plot_boxviolin.m not on path.');
    formats = {'png','svg','fig'};

    mets = {'rate','burst','amp'};
    if isfield(opts,'metrics')&&~isempty(opts.metrics); mets = intersect(mets, lower(opts.metrics),'stable'); end
    chMode = 'mean'; if isfield(opts,'channelMode')&&~isempty(opts.channelMode); chMode = lower(opts.channelMode); end
    switch chMode; case 'separate'; chans = {'G1','G2','G3'}; otherwise; chans = {'G'}; end
    yMax = []; if isfield(opts,'yMax'); yMax = opts.yMax; end

    groups = get_groups(queueFile, opts);

    [~, norm] = compile_mmc(queueFile, opts);

    figs = gobjects(0);
    for k = 1:numel(mets)
        for ci = 1:numel(chans)
            f = bulk_plot_boxviolin(norm, groups, mets{k}, chans{ci}, yMax);
            figs(end+1) = f; %#ok<AGROW>
            if ~isempty(saveDir)
                nm = sprintf('mmc_box_%s_%s', mets{k}, chans{ci});
                if exist('save_one_figure','file')==2
                    save_one_figure(f, saveDir, nm, formats);
                else
                    if ~exist(saveDir,'dir'); mkdir(saveDir); end
                    exportgraphics(f, fullfile(saveDir,[nm '.png']), 'Resolution',150, 'BackgroundColor','white');
                end
            end
        end
    end
    if ~isempty(saveDir); fprintf('[bulk_mmc_boxplots] saved %d figures to %s\n', numel(figs), saveDir); end
end

% ======================================================================
function groups = get_groups(queueFile, opts)
% Prefer explicit opts.groups, then the queue's own saved grouping, else fall
% back to one group holding every distinct condition.
    if isfield(opts,'groups') && ~isempty(opts.groups); groups = opts.groups; return; end
    groups = {};
    try
        Q = load(queueFile,'groups');
        if isfield(Q,'groups') && ~isempty(Q.groups); groups = Q.groups; return; end
    catch; end
    % fallback: single group of all unique conditions
    try
        Q = load(queueFile,'condition');
        if isfield(Q,'condition'); groups = { unique(Q.condition,'stable') }; end
    catch; end
    if isempty(groups); error('bulk_mmc_boxplots: no groups given and none found in the queue.'); end
end
