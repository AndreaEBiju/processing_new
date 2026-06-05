function bulk_save_per_group(out, saveDir, formats)
% BULK_SAVE_PER_GROUP  Render and save the combined plots with ONE IMAGE PER GROUP.
%
%   bulk_save_per_group(out, saveDir)
%   bulk_save_per_group(out, saveDir, {'png','svg','fig'})
%
% For each group it writes, per channel (RVN, LVN):
%   full_<metric>_<chan>_group<g>      full-duration box-violin
%   windowed_<metric>_<chan>_group<g>  1/2/5/10min/full box-violin
%   fano_<chan>_group<g>               Fano-slope box
% and once per group:
%   totalrate_group<g>                 total (LVN+RVN)/s windowed box-violin
% Synergy heatmaps are not group-based (3x3), so they are saved whole:
%   synergy_<metric>_<chan>
%
% `out` is the struct returned by run_pipeline_bulk (needs out.raw, out.norm,
% out.groups). Figures are drawn hidden and closed after saving, so this does
% not disturb anything already on screen.

    if nargin < 3 || isempty(formats); formats = {'png','svg','fig'}; end
    assert(isstruct(out) && isfield(out,'groups') && ~isempty(out.groups), ...
        'bulk_save_per_group: out.groups is empty -- run run_pipeline_bulk first.');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end

    g = out.groups; raw = out.raw; nrm = out.norm;
    chs = {'RVN','LVN'}; metrics = {'rate','vpp','fwhm','cv2'};

    prev = get(0,'DefaultFigureVisible'); set(0,'DefaultFigureVisible','off');
    cu = onCleanup(@() set(0,'DefaultFigureVisible',prev)); %#ok<NASGU>

    for gi = 1:numel(g)
        gg = g(gi);                       % 1x1 cell -> single-group figure
        tag = sprintf('group%d', gi);
        for ci = 1:numel(chs)
            ch = chs{ci};
            for mi = 1:numel(metrics)
                m = metrics{mi};
                f = bulk_plot_boxviolin(nrm, gg, m, ch);
                save_one_figure(f, saveDir, sprintf('full_%s_%s_%s', m, ch, tag), formats); close(f);
                f = bulk_plot_windowed(raw, gg, m, ch);
                save_one_figure(f, saveDir, sprintf('windowed_%s_%s_%s', m, ch, tag), formats); close(f);
            end
            f = bulk_plot_fano_box(nrm, gg, ch);
            save_one_figure(f, saveDir, sprintf('fano_%s_%s', ch, tag), formats); close(f);
        end
        f = bulk_plot_windowed_total(raw, gg);
        save_one_figure(f, saveDir, sprintf('totalrate_%s', tag), formats); close(f);
    end

    % synergy + mean-%-change heatmaps are not grouped -- save one per metric/channel
    for ci = 1:numel(chs)
        for mi = 1:numel(metrics)
            f = bulk_plot_synergy(nrm, metrics{mi}, chs{ci});
            if ~isempty(f) && isgraphics(f)
                save_one_figure(f, saveDir, sprintf('synergy_%s_%s', metrics{mi}, chs{ci}), formats); close(f);
            end
            f = bulk_plot_me_heatmap(raw, metrics{mi}, chs{ci});
            if ~isempty(f) && isgraphics(f)
                save_one_figure(f, saveDir, sprintf('heatmap_%s_%s', metrics{mi}, chs{ci}), formats); close(f);
            end
        end
    end

    fprintf('[bulk] saved per-group figures to %s\n', saveDir);
end
