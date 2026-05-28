function merged = mergeSegments(idx)
    if isempty(idx)
        merged = zeros(0,2);
        return;
    end

    idx = sortrows(idx, 1);
    merged = idx(1,:);

    for i = 2:size(idx,1)
        if idx(i,1) <= merged(end,2) + 1
            merged(end,2) = max(merged(end,2), idx(i,2));
        else
            merged(end+1,:) = idx(i,:); %#ok<AGROW>
        end
    end
end