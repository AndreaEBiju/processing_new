diary(fullfile(pwd,'real_layer12_fit_log.txt'));
diary on;

S = load(fullfile(pwd,'real_rawRows.mat'));
rawRows = S.out.raw;

fprintf('Loaded %d real rawRows.\n', numel(rawRows));

saveDir = fullfile(pwd,'layer12_real_figs');
if ~exist(saveDir,'dir'); mkdir(saveDir); end

fprintf('\n================ PER-CHANNEL (RVN/LVN) FIRST PASS ================\n');
resultsPerChannel = run_layer12_first_pass(rawRows, saveDir);

fprintf('\n================ TOTAL (RVN+LVN) FIRST PASS ================\n');
resultsTotal = run_layer12_total_channel(rawRows, saveDir);

save(fullfile(pwd,'real_layer12_results.mat'), 'resultsPerChannel', 'resultsTotal');
fprintf('\nSaved results to real_layer12_results.mat\n');
diary off;
