function v = bulk_cache_get(cf, field)
% BULK_CACHE_GET  Read one field from the bulk cache (a single struct 'C' stored
% in cf). Returns [] if the file, the struct, or the field is missing.

    v = [];
    if isempty(cf) || ~exist(cf,'file'); return; end
    try
        vars = who('-file', cf);
        if ~ismember('C', vars); return; end
        S = load(cf, 'C');
        if isstruct(S.C) && isfield(S.C, field); v = S.C.(field); end
    catch
    end
end
