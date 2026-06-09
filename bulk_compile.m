function normRows = bulk_compile(rawRows)
% BULK_COMPILE  Baseline-normalize each recovery TRIAL and emit ONE normalized
% row per (recovery file, channel). Repeat trials of the same (animal,
% condition) are kept SEPARATE -- they are independent observations and are NOT
% pooled/averaged together.
%
% Each recovery is paired to ITS OWN baseline by shared file stem (the per-trial
% id that bulk_scan_files computes by stripping the phase suffix; a trial's
% baseline and recovery share it). Falls back to (animal,condition,channel) when
% no stem is available (e.g. an older cache), which is exact in the nominal
% one-trial-per-cell case.
%
% Each output row holds, per the agreed rule (recovery x -> (x-mean_base)/mean_base):
%   .dist.(rate|vpp|fwhm|cv2)   normalized recovery DISTRIBUTION (per-sample)
%   .scalar.(rate|excess|vpp|fwhm|cv2|fano)  normalized scalar effects
% Baseline reference for each metric = that metric's mean in the matched
% baseline file; for fano it is the baseline Fano slope.

    normRows = struct('animal',{},'condition',{},'label',{},'dist',{},'scalar',{});
    if isempty(rawRows); return; end
    hasStem = isfield(rawRows,'stem');

    recI = find(strcmpi({rawRows.phase},'recovery'));
    for r = recI(:)'
        R = rawRows(r);

        % this trial's own baseline: same stem (trial) + same channel
        bsel = [];
        if hasStem && ~isempty(R.stem)
            bsel = find(strcmpi({rawRows.phase},'baseline') ...
                      & strcmp({rawRows.stem}, R.stem) ...
                      & strcmpi({rawRows.label}, R.label));
        end
        if isempty(bsel)   % fallback: pair by (animal,condition,channel)
            bsel = find(match(rawRows, R.animal, R.condition, R.label, 'baseline'));
        end
        if isempty(bsel); continue; end
        B = rawRows(bsel(1));

        nd = struct(); ns = struct();
        for m = {'rate','vpp','fwhm','cv2'}
            mm = m{1}; mu = B.mean.(mm);
            nd.(mm) = norm_vec(R.dist.(mm)(:), mu);    % this trial's recovery samples
            ns.(mm) = norm_scalar(R.mean.(mm), mu);    % this trial's scalar effect
        end
        ns.excess = norm_scalar(R.mean.excess, B.mean.excess);
        ns.fano   = norm_scalar(R.fanoSlope,   B.fanoSlope);

        normRows(end+1) = struct('animal',R.animal,'condition',R.condition, ...
            'label',R.label,'dist',nd,'scalar',ns); %#ok<AGROW>
    end
    fprintf('[compile] %d normalized rows (one per recovery trial x channel; repeats kept separate).\n', ...
        numel(normRows));
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
