close all force;
% RUN_PIPELINE_PAIR  Load baseline + recovery together, process both, compare.
%
% One-shot workflow for a condition pair:
%   1. load_datasets prompts for condition names, then for each you pick its
%      3 files (neural / R-peak / slow-wave) — folders can differ.
%   2. process_dataset runs steps 1 -> 5c on each.
%   3. summaries are saved and the two comparison figures are drawn
%      automatically (no file picker).
%
% Run the whole file, or section by section.

%% Load both datasets
P = pipeline_params();
plotMode = false;   % set true to see every per-step figure for each dataset
                    % (false relies on the comparison plots + step5c verdicts)

Dsets = load_datasets(P, {'baseline', 'recovery'});
if isempty(Dsets); fprintf('Cancelled.\n'); return; end

%% Process each and save summaries
summaryPaths = {};
for i = 1:numel(Dsets)
    fprintf('\n========== Processing condition: %s ==========\n', Dsets{i}.condition);
    Dsets{i} = process_dataset(Dsets{i}, P, plotMode);
    summaryPaths{end+1} = pipeline_save_summary(Dsets{i}, P, Dsets{i}.condition); %#ok<SAGROW>
end

%% Compare (uses the saved summaries directly)
compare_conditions(summaryPaths);      % activity: excess RMS + rate
compare_distributions(summaryPaths);   % fiber question: FWHM / Vpp distributions + KS

% Per-condition modality verdicts are printed by step5c during processing
% (CONTINUUM vs MULTIMODAL). Dsets{i} retains the full D for each condition.
