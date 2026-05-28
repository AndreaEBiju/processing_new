function diagnoseMetricCoverage(cacheFile)
%DIAGNOSEMETRICCOVERAGE  Show which queued sources have which output files.
%
%   diagnoseMetricCoverage()
%     Reads gemsplots_queue.mat next to this function. Prints, for every
%     queued source file × every metric spec, whether the corresponding
%     output file exists on disk (using the same path resolution that
%     loadMetric does — both <stem>_<suffix> and the <stem>_recovery_<suffix>
%     variant for stim_rec sources). Reports per-metric and per-source
%     coverage so you can spot missing analysis runs at a glance.
%
%   diagnoseMetricCoverage(cacheFile)
%     Same, with an explicit cache path.
%
% Uses only dir() (no load), so it's fast even on Google Drive.

    if nargin < 1
        cacheFile = fullfile(fileparts(mfilename('fullpath')), ...
            'gemsplots_queue.mat');
    end
    if ~exist(cacheFile, 'file')
        fprintf('No cache file found at: %s\n', cacheFile);
        return;
    end
    c = load(cacheFile);
    if ~isfield(c,'files') || isempty(c.files)
        fprintf('Cache has no queued files. Run plotAveragedMetrics once first.\n');
        return;
    end
    if ~isfield(c,'metrics') || isempty(c.metrics)
        fprintf('Cache has no metric specs. Run plotAveragedMetrics once to define them.\n');
        return;
    end

    files       = c.files;
    metricSpecs = c.metrics;
    nFiles      = numel(files);
    nMetrics    = numel(metricSpecs);

    % ----- existence matrix -----
    fileFound = false(nFiles, nMetrics);
    candTried = cell(nFiles, nMetrics);
    for i = 1:nFiles
        srcFile = files{i};
        [d, base, ~] = fileparts(srcFile);
        stem = regexprep(base, '_blankmotion$', '');

        for k = 1:nMetrics
            spec = metricSpecs(k);
            suffix = char(spec.suffix);
            if ~endsWith(suffix, '.mat'), suffix = [suffix '.mat']; end

            candidates = {fullfile(d, [base suffix]), ...
                          fullfile(d, [stem suffix])};
            if contains(stem, 'stim_rec', 'IgnoreCase', true) && ...
               ~endsWith(stem, '_recovery', 'IgnoreCase', true) && ...
               ~endsWith(stem, '_stim',     'IgnoreCase', true)
                candidates{end+1} = fullfile(d, [stem '_recovery' suffix]); %#ok<AGROW>
            end

            for ci = 1:numel(candidates)
                if exist(candidates{ci}, 'file')
                    fileFound(i,k) = true;
                    candTried{i,k} = candidates{ci};
                    break;
                end
            end
            if ~fileFound(i,k)
                candTried{i,k} = candidates{end};
            end
        end
    end

    % ----- per-metric summary -----
    fprintf('\n=== Per-metric coverage ===\n');
    fprintf('  %-32s %12s   %s\n', 'Metric', 'found/total', 'pct');
    fprintf('  %s\n', repmat('-', 1, 60));
    for k = 1:nMetrics
        s = metricSpecs(k);
        n = sum(fileFound(:,k));
        pct = 100*n/nFiles;
        bar = repmat('|', 1, max(0,round(pct/2)));
        fprintf('  %-32s %5d / %4d   %5.1f%%  %s\n', ...
            shortLabel(s.label, 32), n, nFiles, pct, bar);
    end

    % ----- per-source matrix (compact) -----
    fprintf('\n=== Per-source matrix (X=present, .=missing) ===\n');
    hdr = repmat(' ', 1, nMetrics);
    for k = 1:nMetrics, hdr(k) = num2str(mod(k,10)); end
    fprintf('  %s  %s\n', hdr, '(columns = metric index, see legend below)');
    for i = 1:nFiles
        row = repmat('.', 1, nMetrics);
        row(fileFound(i,:)) = 'X';
        [~, base, ~] = fileparts(files{i});
        fprintf('  %s  %s\n', row, shortLabel(base, 80));
    end
    fprintf('\nMetric legend:\n');
    for k = 1:nMetrics
        fprintf('  %d = %s  (suffix=%s, field=%s)\n', mod(k,10), ...
            metricSpecs(k).label, metricSpecs(k).suffix, metricSpecs(k).field);
    end

    % ----- show some missing examples -----
    missingPairs = find(~fileFound);
    if ~isempty(missingPairs)
        fprintf('\n=== First 5 missing output files it looked for ===\n');
        for n = 1:min(5, numel(missingPairs))
            idx = missingPairs(n);
            [i, k] = ind2sub([nFiles, nMetrics], idx);
            fprintf('  Source : %s\n', files{i});
            fprintf('  Tried  : %s\n', candTried{i,k});
            fprintf('  Spec   : %s   (suffix=%s)\n', ...
                metricSpecs(k).label, metricSpecs(k).suffix);
            fprintf('\n');
        end
    else
        fprintf('\nAll output files present.\n');
    end

    % ----- sources with zero coverage -----
    zeroSrc = find(~any(fileFound, 2));
    if ~isempty(zeroSrc)
        fprintf('=== Sources with NO output files found (%d / %d) ===\n', ...
            numel(zeroSrc), nFiles);
        for kk = 1:min(10, numel(zeroSrc))
            fprintf('  %s\n', files{zeroSrc(kk)});
        end
        if numel(zeroSrc) > 10
            fprintf('  ... (+%d more)\n', numel(zeroSrc) - 10);
        end
    end
end


function s = shortLabel(str, n)
    str = char(str);
    if numel(str) > n, s = ['...' str(end-(n-3)+1:end)]; else, s = str; end
end
