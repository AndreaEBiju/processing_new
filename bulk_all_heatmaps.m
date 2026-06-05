function bulk_all_heatmaps(out, saveDir, formats)
% BULK_ALL_HEATMAPS  All requested M x E mean-%-change heatmaps in one call.
%
%   bulk_all_heatmaps(out, saveDir)
%
% Nerve metrics (from out.raw, per channel): nerve firing (rate), Vpp, FWHM, CV2.
% HR/HRV metrics (via bulk_hrv_heatmaps -> your loaders): HR, breathing rate,
% HRV, RMSSD, pNN5, sample entropy, SD1, SD2.
%
% Each cell = mean over animals of the per-animal percent change vs baseline.
% saveDir optional (png/svg/fig); requires processing_new on the path.

    if nargin < 2; saveDir = ''; end
    if nargin < 3 || isempty(formats); formats = {'png','svg','fig'}; end
    if ~isempty(saveDir) && ~exist(saveDir,'dir'); mkdir(saveDir); end

    nerve = {'rate','vpp','fwhm','cv2'}; chs = {'RVN','LVN'};
    for ci = 1:numel(chs)
        for mi = 1:numel(nerve)
            f = bulk_plot_me_heatmap(out.raw, nerve{mi}, chs{ci});
            if ~isempty(f) && isgraphics(f) && ~isempty(saveDir) && exist('save_one_figure','file')==2
                save_one_figure(f, saveDir, sprintf('heatmap_%s_%s', nerve{mi}, chs{ci}), formats);
            end
        end
    end

    f = bulk_plot_me_heatmap_total(out.raw);   % total (LVN+RVN)/s
    if ~isempty(f) && isgraphics(f) && ~isempty(saveDir) && exist('save_one_figure','file')==2
        save_one_figure(f, saveDir, 'heatmap_totalrate', formats);
    end

    bulk_hrv_heatmaps(saveDir, formats);   % HR/HRV (opens your buildFileQueue)
end
