diary(fullfile(pwd,'real_harvest_log.txt'));
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

allFiles = bulk_scan_files(folders, conv);

% restrict to exactly the manifest's 90 canonical paths (gemsplots_queue.mat)
manifestPaths = strtrim(strsplit(fileread('/tmp/manifest_paths.txt'), sprintf('\n')));
manifestPaths = manifestPaths(~cellfun(@isempty, manifestPaths));
keep = ismember({allFiles.neural}, manifestPaths);
files = allFiles(keep);
fprintf('[manifest-scoped] %d of %d scanned files kept (manifest has %d paths); %d scanned-but-not-in-manifest dropped.\n', ...
    sum(keep), numel(allFiles), numel(manifestPaths), sum(~keep));
missingFromScan = setdiff(manifestPaths, {allFiles.neural});
if ~isempty(missingFromScan)
    fprintf('WARNING: %d manifest paths not found by the scan:\n', numel(missingFromScan));
    for i = 1:numel(missingFromScan); fprintf('  %s\n', missingFromScan{i}); end
end

opts = struct('conv', conv, 'files', files, 'perFileFigures', false, ...
              'cacheFile', fullfile(pwd,'bulk_cache.mat'));

out = run_pipeline_bulk([], opts);

if isempty(out)
    fprintf('run_pipeline_bulk returned empty.\n');
else
    fprintf('rawRows: %d rows. Channels: %s. Conditions: %s. Animals: %s\n', numel(out.raw), ...
        strjoin(unique({out.raw.label}), ', '), strjoin(unique({out.raw.condition}), ', '), ...
        strjoin(unique({out.raw.animal}), ', '));
    save(fullfile(pwd,'real_rawRows.mat'), 'out');
    fprintf('Saved to real_rawRows.mat\n');
end
diary off;
