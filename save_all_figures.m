function save_all_figures(outDir, prefix, formats)
% SAVE_ALL_FIGURES  Save every open figure to PNG, FIG and SVG.
%
%   save_all_figures(outDir, prefix)
%   save_all_figures(outDir, prefix, {'png','svg','fig'})
%
% Each open figure is written as <prefix>_NN_<figureName>.{png,svg,fig} into
% outDir (created if needed). NN is the order opened. Figures may be hidden;
% they are still rendered correctly.

    if nargin < 2 || isempty(prefix);  prefix  = 'fig'; end
    if nargin < 3 || isempty(formats); formats = {'png','svg','fig'}; end
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    prefix = sanitize(prefix);
    figs = findall(0, 'Type', 'figure');
    figs = flipud(figs(:));                 % oldest first
    if isempty(figs); fprintf('No open figures to save.\n'); return; end

    for i = 1:numel(figs)
        f = figs(i);
        nm = f.Name; if isempty(nm); nm = sprintf('figure%d', f.Number); end
        base = fullfile(outDir, sprintf('%s_%02d_%s', prefix, i, sanitize(nm)));

        whiten_figure(f);                       % white bg + black text + font 20 (permanent)
        if any(strcmpi(formats, 'png'))
            try, exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
            catch, try, saveas(f, [base '.png']); catch, end, end
        end
        if any(strcmpi(formats, 'svg'))
            try, exportgraphics(f, [base '.svg'], 'ContentType', 'vector', 'BackgroundColor','white');
            catch, try, saveas(f, [base '.svg']); catch, end, end
        end
        if any(strcmpi(formats, 'fig'))
            try, savefig(f, [base '.fig']); catch, end
        end
    end
    fprintf('Saved %d figure(s) [%s] to %s\n', numel(figs), strjoin(formats, '/'), outDir);
end

function s = sanitize(s)
    s = regexprep(char(s), '[^\w\-]+', '_');   % keep word chars and hyphen
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
    if isempty(s); s = 'x'; end
end
