function idx = normalizeIdx(idx, N)
    if isempty(idx)
        idx = zeros(0,2);
        return;
    end

    idx = round(idx);
    idx(:,1) = max(idx(:,1), 1);
    idx(:,2) = min(idx(:,2), N);

    bad = idx(:,2) < idx(:,1);
    idx(bad,:) = [];
end
