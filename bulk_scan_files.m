function files = bulk_scan_files(folders, conv)
% BULK_SCAN_FILES  Recursively find baseline/recovery neural files under the
% given folder(s) and pair each with its HRBR file, using the explicit suffixes
% in conv (from bulk_conventions_ui).
%
%   files : struct array with fields
%     neural, hrbr   full paths (hrbr = '' if the sibling is missing)
%     stem           filename with the neural suffix removed
%     animal, condition, phase
%
% Pairing rule: a file named  <STEM><sufXNeural>.mat  is paired with
% <STEM><sufXHRBR>.mat in the same folder (X = baseline or recovery).

    if ischar(folders) || isstring(folders); folders = {char(folders)}; end
    files = struct('neural',{},'hrbr',{},'stem',{},'animal',{},'condition',{},'phase',{});

    % recovery suffix is more specific -> scan it FIRST; the baseline pass then
    % excludes anything ending with the recovery suffix (suffix nesting) AND any
    % 'stim_rec' recording (those are stim/recovery files, never baselines --
    % baselines are the 'bl' recordings).
    %        neuralSuf            hrbrSuf            phase        exclEndsWith          exclContains
    specs = {conv.sufRecNeural,  conv.sufRecHRBR,  'recovery',  '',                   ''; ...
             conv.sufBaseNeural, conv.sufBaseHRBR, 'baseline',  conv.sufRecNeural,    'stim_rec'};

    for fi = 1:numel(folders)
        root = folders{fi};
        for si = 1:size(specs,1)
            sufN = specs{si,1}; sufH = specs{si,2}; phase = specs{si,3};
            exclSuf = specs{si,4}; exclSub = specs{si,5};
            hits = dir(fullfile(root,'**',['*' sufN '.mat']));
            hits = hits(~[hits.isdir]);
            for k = 1:numel(hits)
                name = hits(k).name;
                if ~endsWith(name, [sufN '.mat']); continue; end
                if contains(name,'HRBR','IgnoreCase',true); continue; end
                % don't let a recovery file be picked up by the baseline pass
                if ~isempty(exclSuf) && endsWith(name, [exclSuf '.mat']); continue; end
                % don't let a stim_rec recording be picked up as a baseline
                if ~isempty(exclSub) && contains(name, exclSub, 'IgnoreCase',true); continue; end
                stem = name(1:end-numel([sufN '.mat']));
                neural = fullfile(hits(k).folder, name);
                hrbr   = resolve_hrbr(neural, phase, conv);   %#ok<NASGU> robust HRBR finder
                [ani, cond] = parse_stem(stem);
                files(end+1) = struct('neural',neural,'hrbr',hrbr,'stem',stem, ...
                    'animal',ani,'condition',cond,'phase',phase); %#ok<AGROW>
            end
        end
    end

    % de-duplicate by neural path (a file could match both scans if suffixes nest)
    if ~isempty(files)
        [~, ia] = unique({files.neural}, 'stable'); files = files(ia);
    end
    fprintf('[scan] %d file(s): %d baseline, %d recovery.\n', numel(files), ...
        sum(strcmpi({files.phase},'baseline')), sum(strcmpi({files.phase},'recovery')));
    coverage_report(files);
end

% ----------------------------------------------------------------------
function coverage_report(files)
% list what matched per phase and flag (animal,condition) missing a pair
    fprintf('  --- matched files ---\n');
    for i = 1:numel(files)
        [~, nm] = fileparts(files(i).neural);
        hr = ''; if isempty(files(i).hrbr); hr = '  (!! no HRBR)'; end
        fprintf('    [%-8s] %s / %-8s  %s%s\n', files(i).phase, files(i).animal, ...
            files(i).condition, nm, hr);
    end
    keys = arrayfun(@(f) sprintf('%s|%s', f.animal, f.condition), files, 'UniformOutput', false);
    uk = unique(keys, 'stable');
    fprintf('  --- pairing check ---\n');
    for i = 1:numel(uk)
        sel = strcmp(keys, uk{i});
        hasB = any(sel & strcmpi({files.phase},'baseline'));
        hasR = any(sel & strcmpi({files.phase},'recovery'));
        if ~(hasB && hasR)
            miss = 'baseline'; if hasB; miss = 'recovery'; end
            fprintf('    %s : MISSING %s -> will be dropped at compile\n', uk{i}, miss);
        end
    end
end

% ----------------------------------------------------------------------
function [ani, cond] = parse_stem(stem)
    tok = strsplit(stem,'_');
    if numel(tok)>=2 && ~isempty(tok{2}); ani = lower(tok{2}(1));
    elseif ~isempty(tok) && ~isempty(tok{1}); ani = lower(tok{1}(1));
    else; ani = '?'; end
    cond = '';
    for k = 1:numel(tok)
        if ~isempty(regexp(upper(tok{k}),'^[ME]\d+([ME]\d+)?$','once')); cond = upper(tok{k}); break; end
    end
    if isempty(cond) && ~isempty(tok); cond = tok{1}; end
end
