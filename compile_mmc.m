function [rows, norm] = compile_mmc(queueFile, opts)
% COMPILE_MMC  Read every _mmc.mat listed in the bulk queue and build the same
% raw + baseline-normalized row structures the nerve pipeline uses (bulk_compile),
% so the MMC metrics drop straight into the existing box/violin + M x E machinery.
%
%   [rows, norm] = compile_mmc                     % uses gemsplots_queue.mat in pwd
%   [rows, norm] = compile_mmc(queueFile, opts)
%
% METRICS (per gastric channel, then a channel-mean 'G'):
%   rate  = individual muscle-firing rate  (mmc.firing.avgRate / firing.rate windows)
%   burst = grouped burst-episode rate      (mmc.burst.avgRate  / burst.rate windows)
%   amp   = mean firing peak amplitude       (mmc.firing.peakAmp windows)
%
% ROWS (raw, one per recording x channel-label):
%   .animal .condition .phase('baseline'|'recovery') .stem(trial id) .label
%   .mean.(rate|burst|amp)  scalar per recording
%   .dist.(rate|burst|amp)  the per-window series (finite windows) for violins
% Channel labels: 'G1','G2','G3' and 'G' (mean across the three channels).
%
% NORM (baseline-normalized, one per recovery x label), mirroring bulk_compile:
%   each recovery paired to ITS OWN baseline by trial stem (fallback:
%   animal+condition+label); repeats kept SEPARATE; x -> (x-mean_base)/mean_base.
%   .animal .condition .label .dist.(m) .scalar.(m)
%
% opts: .queueVar (force which queue variable holds the file list)  [auto]

    if nargin < 1 || isempty(queueFile); queueFile = fullfile(pwd,'gemsplots_queue.mat'); end
    if nargin < 2; opts = struct(); end
    assert(exist(queueFile,'file')==2, 'compile_mmc: queue file not found: %s', queueFile);

    Q = load(queueFile);
    [files, animal, condition, phase] = read_queue(Q, opts);
    n = numel(files);
    mets = {'rate','burst','amp'};

    rows = struct('animal',{},'condition',{},'phase',{},'stem',{}, ...
                  'label',{},'mean',{},'dist',{});
    nMiss = 0;
    for i = 1:n
        p = resolve_mmc(files{i});
        if isempty(p); nMiss = nMiss + 1; continue; end
        S = load(p); if ~isfield(S,'mmc'); continue; end
        m = S.mmc;

        [~,b] = fileparts(files{i});
        ph   = norm_phase(phase{i}, b);
        stem = trial_stem(b);
        ani  = char(string(animal{i}));
        cond = char(string(condition{i}));

        % per-channel series + scalars
        ser.rate = cell(1,3); ser.burst = cell(1,3); ser.amp = cell(1,3);
        sca.rate = nan(1,3);  sca.burst = nan(1,3);  sca.amp = nan(1,3);
        for ch = 1:3
            rr = col(m.firing.rate, ch);  ser.rate{ch}  = rr(isfinite(rr));
            bb = col(m.burst.rate,  ch);  ser.burst{ch} = bb(isfinite(bb));
            aa = col(m.firing.peakAmp,ch);ser.amp{ch}   = aa(isfinite(aa));
            sca.rate(ch)  = scalar_avg(m.firing,'avgRate',ch, ser.rate{ch});
            sca.burst(ch) = scalar_avg(m.burst, 'avgRate',ch, ser.burst{ch});
            sca.amp(ch)   = mean(ser.amp{ch},'omitnan');
        end

        % emit G1,G2,G3 rows
        for ch = 1:3
            rows(end+1) = mk_row(ani,cond,ph,stem,sprintf('G%d',ch), ...
                struct('rate',sca.rate(ch),'burst',sca.burst(ch),'amp',sca.amp(ch)), ...
                struct('rate',ser.rate{ch},'burst',ser.burst{ch},'amp',ser.amp{ch})); %#ok<AGROW>
        end
        % channel-mean 'G' row: scalar = mean across channels; dist = pooled windows
        gMean = struct('rate',mean(sca.rate,'omitnan'), ...
                       'burst',mean(sca.burst,'omitnan'),'amp',mean(sca.amp,'omitnan'));
        gDist = struct('rate',vertcat(ser.rate{:}),'burst',vertcat(ser.burst{:}), ...
                       'amp',vertcat(ser.amp{:}));
        rows(end+1) = mk_row(ani,cond,ph,stem,'G',gMean,gDist); %#ok<AGROW>
    end
    if nMiss>0; warning('compile_mmc:miss','%d of %d files had no _mmc.mat cache.', nMiss, n); end

    norm = normalize_rows(rows, mets);
    fprintf('[compile_mmc] %d raw rows -> %d normalized recovery rows (repeats separate).\n', ...
        numel(rows), numel(norm));
end

% ======================================================================
function norm = normalize_rows(rows, mets)
% Pair each recovery to its own baseline (stem+label, fallback animal+cond+label)
% and normalize x -> (x-mean_base)/mean_base, per bulk_compile.
    norm = struct('animal',{},'condition',{},'label',{},'dist',{},'scalar',{});
    if isempty(rows); return; end
    recI = find(strcmpi({rows.phase},'recovery'));
    for r = recI(:)'
        R = rows(r);
        bsel = find(strcmpi({rows.phase},'baseline') & strcmp({rows.stem},R.stem) ...
                  & strcmpi({rows.label},R.label));
        if isempty(bsel)
            bsel = find(strcmpi({rows.phase},'baseline') & strcmpi({rows.animal},R.animal) ...
                      & strcmpi({rows.condition},R.condition) & strcmpi({rows.label},R.label));
        end
        if isempty(bsel); continue; end
        B = rows(bsel(1));
        nd = struct(); ns = struct();
        for mm = mets
            k = mm{1}; mu = B.mean.(k);
            nd.(k) = norm_vec(R.dist.(k), mu);
            ns.(k) = norm_scalar(R.mean.(k), mu);
        end
        norm(end+1) = struct('animal',R.animal,'condition',R.condition, ...
            'label',R.label,'dist',nd,'scalar',ns); %#ok<AGROW>
    end
end

% ======================================================================
function [files, animal, condition, phase] = read_queue(Q, opts)
% Queue is saved as parallel cells files/animal/condition/phase.
    if isfield(opts,'queueVar') && ~isempty(opts.queueVar) && isfield(Q,opts.queueVar)
        files = Q.(opts.queueVar);
    elseif isfield(Q,'files'); files = Q.files;
    else
        fn = fieldnames(Q); files = [];
        for i=1:numel(fn); v=Q.(fn{i}); if iscell(v)&&numel(v)>=1&&ischar(v{1}); files=v; break; end; end
        assert(~isempty(files),'compile_mmc: could not find the file-list cell in the queue.');
    end
    get = @(nm) getcell(Q, nm, numel(files));
    animal    = get('animal');
    condition = get('condition');
    phase     = get('phase');
end

function c = getcell(Q, nm, n)
    if isfield(Q,nm) && iscell(Q.(nm)) && numel(Q.(nm))==n; c = Q.(nm);
    else; c = repmat({''},n,1); end
end

function p = resolve_mmc(entry)
% _mmc.mat sits next to the source file; find it by stem glob.
    [d,b,e] = fileparts(entry); if isempty(d); d = pwd; end
    stem = regexprep([b e],'(_blankmotion|_HRBR)?\.mat$','');
    c = dir(fullfile(d,[stem '*_mmc.mat']));
    if isempty(c); p = ''; else; p = fullfile(c(1).folder, c(1).name); end
end

function ph = norm_phase(pIn, basename)
% Canonicalize phase to 'baseline'/'recovery' from the queue field or filename.
    s = lower(char(string(pIn)));
    if contains(s,'rec');  ph='recovery'; return; end
    if contains(s,'base')||strcmp(s,'bl'); ph='baseline'; return; end
    bl = lower(basename);
    if contains(bl,'stim_rec')||contains(bl,'_rec_'); ph='recovery';
    elseif contains(bl,'_bl_')||contains(bl,'base');   ph='baseline';
    else; ph='recovery'; end   % default: treat unknown as recovery (kept, may fail to pair)
end

function stem = trial_stem(basename)
% Trial id shared by a baseline/recovery pair = everything before the phase token.
    stem = regexprep(basename,'_(stim_rec|stim_bl|recovery|baseline|rec|bl)_.*$','','ignorecase');
end

function v = col(A, ch)
    if isempty(A) || size(A,2) < ch; v = []; else; v = double(A(:,ch)); end
end

function s = scalar_avg(lvl, fld, ch, ser)
% Prefer the stored whole-record avgRate; fall back to the mean of the windows.
    if isfield(lvl,fld) && numel(lvl.(fld))>=ch && isfinite(lvl.(fld)(ch)); s = lvl.(fld)(ch);
    else; s = mean(ser,'omitnan'); end
end

function r = mk_row(ani,cond,ph,stem,label,meanS,distS)
    r = struct('animal',ani,'condition',cond,'phase',ph,'stem',stem, ...
               'label',label,'mean',meanS,'dist',distS);
end

function y = norm_vec(x, mu)
    if ~isfinite(mu) || mu==0; y = []; return; end
    y = (x(isfinite(x)) - mu)/mu;
end

function v = norm_scalar(x, mu)
    if ~isfinite(mu) || mu==0 || ~isfinite(x); v = NaN; return; end
    v = (x - mu)/mu;
end
