function normRows = bulk_compile(rawRows)
% BULK_COMPILE  Pair baseline/recovery per (animal, condition, channel) and
% baseline-normalize, per the agreed rule: every recovery sample x becomes
% (x - mean_baseline)/mean_baseline. Produces one normalized row per
% (animal, condition, channel) holding:
%   .dist.(rate|vpp|fwhm|cv2)   normalized recovery DISTRIBUTIONS (per-sample)
%   .scalar.(rate|excess|vpp|fwhm|cv2|fano)  normalized scalar effects
%
% Baseline mean for each metric = that metric's mean in the matched baseline
% file; for fano it is the baseline Fano slope.

    normRows = struct('animal',{},'condition',{},'label',{}, ...
        'dist',{},'scalar',{});
    if isempty(rawRows); return; end

    keyRec = find(strcmpi({rawRows.phase},'recovery'));
    seen = {};
    for r = keyRec
        R = rawRows(r);
        key = sprintf('%s|%s|%s', R.animal, R.condition, R.label);
        if any(strcmp(seen,key)); continue; end
        seen{end+1} = key; %#ok<AGROW>

        recIdx  = find(match(rawRows, R.animal, R.condition, R.label, 'recovery'));
        baseIdx = find(match(rawRows, R.animal, R.condition, R.label, 'baseline'));
        if isempty(baseIdx); continue; end
        B = rawRows(baseIdx(1));

        nd = struct(); ns = struct();
        for m = {'rate','vpp','fwhm','cv2'}
            mm = m{1};
            mu = B.mean.(mm);
            % pool recovery distributions across any repeat trials, normalize
            d = [];
            for ri = recIdx(:)'
                d = [d; rawRows(ri).dist.(mm)(:)]; %#ok<AGROW>
            end
            nd.(mm) = norm_vec(d, mu);
            % normalized scalar = (mean recovery - mu)/mu (pool trial means)
            rm = mean(arrayfun(@(ri) rawRows(ri).mean.(mm), recIdx), 'omitnan');
            ns.(mm) = norm_scalar(rm, mu);
        end
        % excess (scalar only)
        muE = B.mean.excess;
        rmE = mean(arrayfun(@(ri) rawRows(ri).mean.excess, recIdx), 'omitnan');
        ns.excess = norm_scalar(rmE, muE);
        % fano slope (scalar only)
        muF = B.fanoSlope;
        rmF = mean(arrayfun(@(ri) rawRows(ri).fanoSlope, recIdx), 'omitnan');
        ns.fano = norm_scalar(rmF, muF);

        normRows(end+1) = struct('animal',R.animal,'condition',R.condition, ...
            'label',R.label,'dist',nd,'scalar',ns); %#ok<AGROW>
    end
    fprintf('[compile] %d normalized (animal x condition x channel) rows.\n', numel(normRows));
end

% ----------------------------------------------------------------------
function tf = match(rows, ani, cond, lab, phase)
    tf = strcmpi({rows.animal},ani) & strcmpi({rows.condition},cond) ...
        & strcmpi({rows.label},lab) & strcmpi({rows.phase},phase);
end

function y = norm_vec(x, mu)
    if ~isfinite(mu) || mu==0; y = []; return; end
    y = (x(isfinite(x)) - mu) / mu;
end

function v = norm_scalar(x, mu)
    if ~isfinite(mu) || mu==0 || ~isfinite(x); v = NaN; return; end
    v = (x - mu)/mu;
end
