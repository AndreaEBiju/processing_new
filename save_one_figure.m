function save_one_figure(fig, outDir, name, formats)
% SAVE_ONE_FIGURE  Save a single figure handle to outDir as name.<ext> for each
% requested format (default png, svg, fig).

    if nargin < 4 || isempty(formats); formats = {'png','svg','fig'}; end
    if ~exist(outDir,'dir'); mkdir(outDir); end
    safe = regexprep(name, '[^\w\-.]', '_');
    base = fullfile(outDir, safe);
    whiten_figure(fig);                          % white bg + black text + font 20 (permanent)
    for i = 1:numel(formats)
        try
            switch lower(formats{i})
                case 'fig'
                    savefig(fig, [base '.fig']);
                case 'png'
                    exportgraphics(fig, [base '.png'], 'Resolution', 200, 'BackgroundColor', 'white');
                case 'svg'
                    try
                        exportgraphics(fig, [base '.svg'], 'ContentType', 'vector', 'BackgroundColor', 'white');
                    catch
                        print(fig, [base '.svg'], '-dsvg', '-vector');
                    end
            end
        catch ME
            warning('save_one_figure:fmt', '%s save failed (%s)', formats{i}, ME.message);
        end
    end
end
