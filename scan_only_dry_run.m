diary(fullfile(pwd,'scan_dry_run_log.txt'));
diary on;

conv = struct('sufBaseNeural','_v0.2.2_blankmotion', ...
              'sufRecNeural','_v0.2.2_recovery_blankmotion', ...
              'sufBaseHRBR','_v0.2.2_blankmotion_HRBR', ...
              'sufRecHRBR','_v0.2.2_recovery_HRBR', ...
              'rvnCol',1, 'lvnCol',2, 'heartVar','heartlocs');

survivalsRoot = '/Users/andreaelizabethbiju/Library/CloudStorage/GoogleDrive-andreabiju@g.harvard.edu/Shared drives/BIONICs Lab Workspace/Project Folders/GEMS/Survivals';
dateFolders = {'05062026','05082026','05092026','05102026','05112026','05122026', ...
               '05132026','05142026','05152026','05202026','05252026','06192026'};
folders = cellfun(@(d) fullfile(survivalsRoot, d), dateFolders, 'UniformOutput', false);
folders{end+1} = fullfile(survivalsRoot, 'TDT', '0505-0506-0507');

fprintf('Scanning %d date folders...\n', numel(folders));
files = bulk_scan_files(folders, conv);

fprintf('\n--- vengmetrics cache presence check ---\n');
nHave = 0; nMiss = 0;
for i = 1:numel(files)
    [d,b] = fileparts(files(i).neural);
    pf = fullfile(d, [b '_vengmetrics.mat']);
    if exist(pf,'file'); nHave = nHave+1; else; nMiss = nMiss+1; fprintf('  MISSING cache: %s\n', files(i).stem); end
end
fprintf('%d have per-file cache, %d missing (would need fresh processing).\n', nHave, nMiss);

save(fullfile(pwd,'scan_dry_run_files.mat'), 'files', 'conv');
diary off;
