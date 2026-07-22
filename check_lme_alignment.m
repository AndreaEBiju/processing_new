S = load('real_stage3_results.mat');
R = load('real_rawRows.mat');
rawRows = R.out.raw;
thetaShape = S.thetaShapeAll.RVN;

isRec = strcmpi({rawRows.phase}, 'recovery');
isCh  = strcmpi({rawRows.label}, 'RVN');
uM = nan(numel(rawRows),1); uE = nan(numel(rawRows),1);
for i = 1:numel(rawRows)
    c = upper(strtrim(rawRows(i).condition));
    m = regexp(c,'M(\d+)','tokens','once'); if isempty(m); mm=0; else; mm=str2double(m{1}); end
    e = regexp(c,'E(\d+)','tokens','once'); if isempty(e); ee=0; else; ee=str2double(e{1}); end
    uM(i) = mm; uE(i) = ee;
end
isCombined = (uM > 0) & (uE > 0);
sel = find(isRec(:) & isCh(:) & isCombined(:));

uCells = unique([uM(sel) uE(sel)], 'rows');
fprintf('%-12s %8s %8s %8s %8s %8s\n', 'condition', 'obsCV2', 'hatCV2', 'obsFWHM', 'hatFWHM', 'p_2');
for c = 1:size(uCells,1)
    m = uCells(c,1); e = uCells(c,2);
    idx = sel(uM(sel)==m & uE(sel)==e);
    obsCV2 = mean(arrayfun(@(i) rawRows(i).mean.cv2, idx));
    obsFWHM = mean(arrayfun(@(i) rawRows(i).mean.fwhm, idx));
    hatCV2 = model_layer12_equations('cv2_hat', m, e, thetaShape);
    hatFWHM = model_layer12_equations('fwhm_hat', m, e, thetaShape);
    [~,~,p2] = model_layer12_equations('term_shares', m, e, thetaShape);
    marker = '';
    if m==100 && e==10; marker = '  <-- LME: FWHM(RVN) significant +65.2%%'; end
    if m==50 && e==100; marker = '  <-- LME: CV2(RVN) significant +38.1%%'; end
    fprintf('M%-3dE%-6d %8.4f %8.4f %8.4f %8.4f %8.4f%s\n', m, e, obsCV2, hatCV2, obsFWHM, hatFWHM, p2, marker);
end
