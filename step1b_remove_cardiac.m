function D = step1b_remove_cardiac(D, P, plotMode)
% STEP1B_REMOVE_CARDIAC  R-peak-aligned QRS template subtraction.
%
%   D = step1b_remove_cardiac(D, P, plotMode)
%
% The sharp cardiac QRS survives the 100 Hz high-pass as a large biphasic
% transient and dominates detection. This step removes it BEFORE the noise
% estimate and detection, using the supplied R-peaks:
%
%   1. For each R-peak, re-align to the local |signal| maximum within
%      +/- P.cardiacSearchMs (handles a fixed lag / jitter between the
%      cardiac file and the neural channel).
%   2. Extract a window [-P.cardiacPreMs, +P.cardiacPostMs] around each
%      fiducial; for each beat the template is the median of the nearest
%      P.cardiacTemplateBeats beats (a rolling template that tracks slow
%      drift in QRS shape over the session), robust to the occasional real
%      spike riding on a beat.
%   3. Subtract a least-squares-scaled copy of that template from each beat,
%      so respiration-modulated beats of differing amplitude subtract cleanly
%      while non-stereotyped neural events (incl. baroreceptor afferents) are
%      preserved.
%
% Beats whose window amplitude is far below the median (the dead recording
% error) are skipped, so the flat region is left untouched and step 2's
% validity mask still flags it.
%
% Requires D.filtered (step 1) and D.rpeakSamples (step 0).
%
% Replaces D.filtered with the cleaned signal and keeps:
%   D.filteredRaw   pre-subtraction bandpassed signal (samples x nNeural)
%   D.cardiac       1 x nNeural struct: template, t_ms, nBeatsUsed,
%                   Vpp_uv, fidTimesUsed
%
% plotMode (single-file): per channel, the QRS template over beats, and a
% representative window before vs after subtraction with R-peaks marked.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    if ~isfield(D, 'filtered')
        error('step1b_remove_cardiac:noFiltered', 'Run step1_bandpass first.');
    end
    if ~isfield(D, 'rpeakSamples') || isempty(D.rpeakSamples)
        warning('step1b_remove_cardiac:noRpeaks', ...
            'No R-peaks available; skipping cardiac subtraction.');
        D.filteredRaw = D.filtered;
        return;
    end

    fs   = D.fs;
    N    = size(D.filtered, 1);
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;

    preS    = round(P.cardiacPreMs    * 1e-3 * fs);
    postS   = round(P.cardiacPostMs   * 1e-3 * fs);
    searchS = round(P.cardiacSearchMs * 1e-3 * fs);
    winLen  = preS + postS + 1;
    t_ms    = (-preS:postS) / fs * 1e3;

    rs = round(D.rpeakSamples(:));
    rs = rs(rs >= 1 & rs <= N);
    nB = numel(rs);

    D.filteredRaw = D.filtered;
    cleaned = D.filtered;
    D.cardiac = repmat(struct('template', [], 't_ms', t_ms, 'nBeatsUsed', 0, ...
        'templateBeats', 0, 'Vpp_uv', NaN, 'fidTimesUsed', []), 1, nCh);

    fprintf('[step1b] Cardiac subtraction: %d R-peaks, window [-%g, +%g] ms.\n', ...
        nB, P.cardiacPreMs, P.cardiacPostMs);

    for k = 1:nCh
        xf = D.filtered(:, k);
        xc = xf;

        % --- pass 1: per-beat fiducial, window, amplitude ---
        fid     = nan(nB, 1);
        beatMax = nan(nB, 1);
        W       = nan(nB, winLen);
        for i = 1:nB
            r  = rs(i);
            lo = max(1, r - searchS);
            hi = min(N, r + searchS);
            seg = xf(lo:hi);
            if all(~isfinite(seg)); continue; end
            [~, m] = max(abs(seg));
            f = lo + m - 1;
            wlo = f - preS; whi = f + postS;
            if wlo < 1 || whi > N; continue; end
            w = xf(wlo:whi);
            if any(~isfinite(w)); continue; end
            fid(i)     = f;
            W(i, :)    = w(:).';
            beatMax(i) = max(abs(w));
        end

        medMax = median(beatMax, 'omitnan');
        good   = isfinite(fid) & isfinite(beatMax) & ...
                 (beatMax > P.cardiacMinActiveFrac * medMax);

        if nnz(good) < 5
            warning('step1b:fewBeats', ...
                'ch %d (%s): only %d usable beats; leaving signal unchanged.', ...
                ch(k), lab{k}, nnz(good));
            D.cardiac(k).template = nan(1, winLen);
            continue;
        end

        % --- global template (for display only) ---
        tmplGlobal = median(W(good, :), 1, 'omitnan').';   % winLen x 1

        % --- pass 2: subtract a LOCAL (rolling) template from each good beat ---
        % For each beat, the template is the median of the nearest K good beats
        % (in time), so a QRS shape that drifts over the session is tracked and
        % subtracted cleanly rather than leaving a residual.
        gi = find(good);
        nG = numel(gi);
        K  = max(5, min(P.cardiacTemplateBeats, nG));
        half = floor(K / 2);
        for j = 1:nG
            lo = max(1, j - half);
            hi = min(nG, lo + K - 1);
            lo = max(1, hi - K + 1);                 % keep a full K-beat window
            tmpl = median(W(gi(lo:hi), :), 1, 'omitnan').';
            den  = tmpl.' * tmpl;
            f = fid(gi(j));
            wlo = f - preS; whi = f + postS;
            w = xf(wlo:whi);
            if den > 0
                scale = (w(:).' * tmpl) / den;
            else
                scale = 1;
            end
            xc(wlo:whi) = w - scale * tmpl;
        end

        cleaned(:, k) = xc;
        D.cardiac(k).template      = tmplGlobal.';
        D.cardiac(k).nBeatsUsed    = nG;
        D.cardiac(k).templateBeats = K;
        D.cardiac(k).Vpp_uv        = (max(tmplGlobal) - min(tmplGlobal)) * 1e6;
        D.cardiac(k).fidTimesUsed  = (fid(good) - 1) / fs;

        fprintf('        ch %d (%s): %d beats | rolling template of %d beats | Vpp = %.1f uV\n', ...
            ch(k), lab{k}, nG, K, D.cardiac(k).Vpp_uv);
    end

    D.filtered = cleaned;

    if plotMode
        for k = 1:nCh
            plot_cardiac(D, k, t_ms);
        end
    end
end

% ========================================================================
function plot_cardiac(D, k, t_ms)
    fs   = D.fs;
    t    = D.t;
    raw  = D.filteredRaw(:, k);
    cln  = D.filtered(:, k);
    cd   = D.cardiac(k);
    lab  = D.channelLabels{k};

    fig = figure('Color', 'w', 'Name', sprintf('Step 1b — %s', lab), ...
        'Position', [150 120 1150 760]);
    tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 1b cardiac subtraction  |  %s  |  %d beats', ...
        lab, cd.nBeatsUsed), 'Interpreter', 'none');

    % ---- template + beat snippets ----
    nexttile;
    if isempty(cd.fidTimesUsed) || all(isnan(cd.template))
        text(0.4, 0.5, 'No usable beats', 'Units', 'normalized'); axis off;
    else
        preS = round((0 - t_ms(1)) / 1e3 * fs);
        postS = numel(t_ms) - preS - 1;
        fids = round(cd.fidTimesUsed * fs) + 1;
        nShow = min(150, numel(fids));
        sel = round(linspace(1, numel(fids), nShow));
        hold on;
        for ii = sel
            f = fids(ii); wlo = f - preS; whi = f + postS;
            if wlo >= 1 && whi <= numel(raw)
                plot(t_ms, raw(wlo:whi) * 1e6, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.3);
            end
        end
        plot(t_ms, cd.template * 1e6, 'k', 'LineWidth', 2);
        grid on; xlabel('Time from R-fiducial (ms)'); ylabel('Amplitude (\muV)');
        title(sprintf('QRS template (Vpp = %.0f \\muV) over %d beats', cd.Vpp_uv, nShow));
        xlim([t_ms(1) t_ms(end)]);
    end

    % ---- before vs after over a representative window ----
    nexttile;
    if isempty(cd.fidTimesUsed)
        c = t(round(numel(t)/2));
    else
        c = median(cd.fidTimesUsed);
    end
    lo = max(0, c - 2.5); hi = min(t(end), lo + 5);
    i0 = max(1, floor(lo*fs)+1); i1 = min(numel(t), floor(hi*fs)+1);
    plot(t(i0:i1), raw(i0:i1)*1e6, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, ...
        'DisplayName', 'before'); hold on;
    plot(t(i0:i1), cln(i0:i1)*1e6, 'b', 'LineWidth', 0.5, 'DisplayName', 'after');
    if isfield(D, 'rpeakTimes') && ~isempty(D.rpeakTimes)
        rp = D.rpeakTimes(D.rpeakTimes >= t(i0) & D.rpeakTimes <= t(i1));
        yl = ylim;
        for rr = 1:numel(rp)
            plot([rp(rr) rp(rr)], yl, ':', 'Color', [0 0.45 1], ...
                'LineWidth', 0.7, 'HandleVisibility', 'off');
        end
        if ~isempty(rp)
            plot(nan, nan, ':', 'Color', [0 0.45 1], 'DisplayName', 'R-peaks');
        end
    end
    grid on; xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    title('Representative window: before vs after cardiac subtraction');
    legend('show', 'Location', 'northeast');
    xlim([t(i0) t(i1)]);
end
