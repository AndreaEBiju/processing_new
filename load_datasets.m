function Dsets = load_datasets(P, conditionNames)
% LOAD_DATASETS  Load several datasets (e.g. baseline + recovery) in one go.
%
%   Dsets = load_datasets(P)                    % prompts for condition names
%   Dsets = load_datasets(P, {'baseline','recovery'})
%
% Loops the step-0 loader once per condition (each can live in a different
% folder), tags each D with D.condition, and returns a cell array of D
% structs. Returns {} if cancelled.

    if nargin < 1 || isempty(P); P = pipeline_params(); end
    if nargin < 2 || isempty(conditionNames)
        ans_ = inputdlg('Condition names (comma-separated):', 'Conditions', ...
            [1 60], {'baseline,recovery'});
        if isempty(ans_); Dsets = {}; return; end
        conditionNames = strtrim(strsplit(ans_{1}, ','));
        conditionNames = conditionNames(~cellfun(@isempty, conditionNames));
    end

    Dsets = {};
    for i = 1:numel(conditionNames)
        cn = conditionNames{i};
        uiwait(msgbox(sprintf('Next: select the files for condition  "%s".', cn), ...
            'Load dataset', 'modal'));
        D = step0_load_data(P, sprintf('  [%s]', cn));
        if isempty(D)
            fprintf('Cancelled while loading "%s".\n', cn);
            Dsets = {}; return;
        end
        D.condition = cn;
        Dsets{end+1} = D; %#ok<AGROW>
    end
    fprintf('Loaded %d dataset(s): %s\n', numel(Dsets), strjoin(conditionNames, ', '));
end
