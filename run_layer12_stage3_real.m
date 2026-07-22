diary(fullfile(pwd,'real_stage3_log.txt'));
diary on;

S = load(fullfile(pwd,'real_layer12_results.mat'));
R = load(fullfile(pwd,'real_rawRows.mat'));
rawRows = R.out.raw;

pc = S.resultsPerChannel.perChannelDetail;
channels = {'RVN','LVN'};
thetaShapeAll = struct();
for k = 1:numel(channels)
    ch = channels{k};
    fprintf('\n================ Stage 3 (CV2/FWHM term-share fit): %s ================\n', ch);
    thetaFull = pc.(ch).thetaFull;
    thetaShapeAll.(ch) = fit_layer12_stage3_shape(rawRows, ch, thetaFull);
end

fprintf('\n================ TOTAL channel Stage 3 ================\n');
Rtotal = S.resultsTotal;
% run_layer12_total_channel.m builds its own TOTAL rows internally; replicate
% that here since fit_layer12_stage3_shape needs TOTAL-labeled rawRows, not
% RVN/LVN ones. Sum RVN+LVN mean.rate/cv2/fwhm per matched trial (same
% grouping convention as build_total_rate_table in bulk_mixed_models.m).
totalRows = build_total_rows_local(rawRows);
thetaShapeAll.TOTAL = fit_layer12_stage3_shape(totalRows, 'TOTAL', Rtotal.thetaFull);

save(fullfile(pwd,'real_stage3_results.mat'), 'thetaShapeAll');
fprintf('\nSaved to real_stage3_results.mat\n');
diary off;

function totalRows = build_total_rows_local(rawRows)
    isRVN = strcmpi({rawRows.label},'RVN'); isLVN = strcmpi({rawRows.label},'LVN');
    rvn = rawRows(isRVN); lvn = rawRows(isLVN);
    keys = arrayfun(@(r) sprintf('%s|%s|%s|%s', r.animal, r.condition, r.phase, r.stem), rvn, 'UniformOutput', false);
    keysL = arrayfun(@(r) sprintf('%s|%s|%s|%s', r.animal, r.condition, r.phase, r.stem), lvn, 'UniformOutput', false);
    totalRows = rvn; totalRows = totalRows(1:0);
    for i = 1:numel(rvn)
        j = find(strcmp(keysL, keys{i}), 1);
        if isempty(j); continue; end
        r = rvn(i);
        r.label = 'TOTAL';
        r.mean.rate = rvn(i).mean.rate + lvn(j).mean.rate;
        r.mean.cv2  = mean([rvn(i).mean.cv2, lvn(j).mean.cv2], 'omitnan');
        r.mean.fwhm = mean([rvn(i).mean.fwhm, lvn(j).mean.fwhm], 'omitnan');
        totalRows(end+1) = r; %#ok<AGROW>
    end
    fprintf('[total] %d TOTAL rows built for Stage 3 (paired RVN+LVN; rate summed, CV2/FWHM averaged).\n', numel(totalRows));
end
