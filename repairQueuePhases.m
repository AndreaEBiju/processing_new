function repairQueuePhases(cacheFile)
%REPAIRQUEUEPHASES  Re-derive the phase column of an existing queue cache
% from each file's basename, applying the current buildFileQueue rule
% ('recovery' if the name contains 'recovery' OR 'stim_rec', else
% 'baseline'). Use this after a parser change to avoid having to rebuild
% the queue from scratch.
%
%   repairQueuePhases()                use the default cache path
%   repairQueuePhases('/path/to.mat')  use a specific cache
%
% Prints before/after counts so you can confirm the change.

    if nargin < 1
        cacheFile = fullfile(fileparts(mfilename('fullpath')), ...
            'gemsplots_queue.mat');
    end
    if ~exist(cacheFile, 'file')
        fprintf('No cache file at %s\n', cacheFile);
        return;
    end
    s = load(cacheFile);
    if ~isfield(s,'files') || isempty(s.files)
        fprintf('Cache is empty — nothing to repair.\n');
        return;
    end

    n = numel(s.files);
    if ~isfield(s,'phase') || numel(s.phase) ~= n
        s.phase = repmat({'baseline'}, n, 1);
    end

    fprintf('Before: %s\n', tabPhase(s.phase));
    nChanged = 0;
    for i = 1:n
        [~, base, ~] = fileparts(s.files{i});
        if contains(base, 'recovery', 'IgnoreCase', true) || ...
           contains(base, 'stim_rec', 'IgnoreCase', true)
            newPh = 'recovery';
        elseif contains(base, 'baseline','IgnoreCase', true) || ...
               ~isempty(regexp(base, '(^|_)bl(_|$)', 'once'))
            newPh = 'baseline';
        else
            newPh = 'baseline';
        end
        if ~strcmpi(s.phase{i}, newPh)
            s.phase{i} = newPh;
            nChanged = nChanged + 1;
        end
    end
    fprintf('After:  %s\n', tabPhase(s.phase));
    fprintf('Updated %d / %d rows.\n', nChanged, n);

    save(cacheFile, '-struct', 's');
end

function s = tabPhase(phases)
    u = unique(phases);
    parts = cell(numel(u),1);
    for k = 1:numel(u)
        n = sum(strcmpi(phases, u{k}));
        parts{k} = sprintf('%s=%d', u{k}, n);
    end
    s = strjoin(parts, ', ');
end
