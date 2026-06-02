function bulk_cache_set(cf, field, value)
% BULK_CACHE_SET  Write one field into the bulk cache (a single struct 'C' in
% cf), preserving the other fields. Atomic-ish: load -> update -> save whole.

    if isempty(cf); return; end
    C = struct();
    if exist(cf,'file')
        try
            vars = who('-file', cf);
            if ismember('C', vars); S = load(cf,'C'); if isstruct(S.C); C = S.C; end; end
        catch
        end
    end
    C.(field) = value;
    try
        save(cf, 'C');
    catch ME
        warning('bulk:cacheset', 'cache write failed (%s): %s', field, ME.message);
    end
end
