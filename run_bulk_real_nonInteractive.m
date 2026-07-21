conv = struct('sufBaseNeural','_v0.2.2_blankmotion', ...
              'sufRecNeural','_v0.2.2_recovery_blankmotion', ...
              'sufBaseHRBR','_v0.2.2_blankmotion_HRBR', ...
              'sufRecHRBR','_v0.2.2_recovery_HRBR', ...
              'rvnCol',1, 'lvnCol',2, 'heartVar','heartlocs');

root = '/Users/andreaelizabethbiju/Library/CloudStorage/GoogleDrive-andreabiju@g.harvard.edu/Shared drives/BIONICs Lab Workspace/Project Folders/GEMS/Survivals';

opts = struct('conv', conv, 'folders', {{root}}, 'perFileFigures', false, ...
              'cacheFile', fullfile(pwd,'bulk_cache.mat'));

out = run_pipeline_bulk([], opts);

if isempty(out)
    fprintf('run_pipeline_bulk returned empty.\n');
else
    fprintf('rawRows: %d rows. Channels: %s. Conditions: %s\n', numel(out.raw), ...
        strjoin(unique({out.raw.label}), ', '), strjoin(unique({out.raw.condition}), ', '));
    save(fullfile(pwd,'real_rawRows.mat'), 'out');
    fprintf('Saved to real_rawRows.mat\n');
end
