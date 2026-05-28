function localSegs = fullSegmentsToLocalSegments(fullSegs, mask)

    if isempty(fullSegs)
        localSegs = zeros(0,2);
        return;
    end

    mask = mask(:);
    fullToLocal = zeros(numel(mask),1);

    localSamples = find(mask);
    fullToLocal(localSamples) = 1:numel(localSamples);

    localSegs = zeros(0,2);

    for k = 1:size(fullSegs,1)
        s1 = max(1, fullSegs(k,1));
        s2 = min(numel(mask), fullSegs(k,2));

        overlap = s1:s2;
        overlap = overlap(mask(overlap));

        if isempty(overlap)
            continue;
        end

        localSegs(end+1,:) = ...
            [fullToLocal(overlap(1)), fullToLocal(overlap(end))]; %#ok<AGROW>
    end

    localSegs = mergeSegments(localSegs);
end