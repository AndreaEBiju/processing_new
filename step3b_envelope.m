function D = step3b_envelope(D, P, plotMode)
% STEP3B_ENVELOPE  Threshold-free spike-band activity envelope (RMS).
%
%   D = step3b_envelope(D, P, plotMode)
%
% Computes the root-mean-square of the cardiac-cleaned spike-band signal in
% time bins (P.envBinSec) over valid samples. This is the PRIMARY activity
% readout for this low-SNR, C-fiber-dominated recording: it captures the
% whole population, including sub-threshold events that discrete detection
% (step 3) misses, and does not depend on a detection threshold.
%
% A guard window of +/- P.envCardiacGuardMs around each R-peak is excluded so
% that any residual cardiac left after step 1b does not inflate the activity
% estimate. (Mean RMS, not MAD, is used so genuine spikes DO count toward the
% activity -- the guard window handles the residual instead.)
%
% Requires D.filtered (cardiac-cleaned, step 1b) and D.validMask (step 2).
%
% Adds to D:
%   D.envelope   1 x nNeural struct, per channel:
%       .t            bin-centre times (s)
%       .rms_uv       per-bin RMS (uV)
%       .validFrac    fraction of each bin that contributed
%       .meanRMS_uv   mean RMS over well-sampled bins  <- scalar activity measure
%       .guardMs
%
% For cross-condition comparison, meanRMS_uv (and the rms_uv trace) is the
% headline "how much is the nerve firing" measure. Remember it is sensitive
% to the noise floor, so carry sigma as a covariate / normalise within
% electrode, as planned.
%
% plotMode (single-file): per channel, the RMS activity trace over the whole
% recording (invalid shaded), plus a representative window of the cleaned
% signal with its continuous envelope overlaid.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end
    for f = {'filtered','validMask','sigma'}
        if ~isfield(D, f{1})
            error('step3b_envelope:missing', 'Run steps 1-2 first (missing D.%s).', f{1});
        end
    end

    fs   = D.fs;
    N    = size(D.filtered, 1);
    ch   = D.neuralChannels;
    nCh  = numel(ch);
    lab  = D.channelLabels;

    guard  = round(P.envCardiacGuardMs * 1e-3 * fs);
    binN   = max(1, round(P.envBinSec * fs));
    smN    = max(1, round(P.envSmoothSec * fs));

    % R-peak guard mask (shared across channels)
    guardMask = false(N, 1);
    if isfield(D, 'rpeakSamples') && ~isempty(D.rpeakSamples)
        rs = round(D.rpeakSamples(:));
        rs = rs(rs >= 1 & rs <= N);
        for i = 1:numel(rs)
            guardMask(max(1, rs(i)-guard) : min(N, rs(i)+guard)) = true;
        end
    end

    D.envelope = repmat(struct('t', [], 'rms_uv', [], 'sigmaFloor_uv', [], ...
        'excess_uv', [], 'validFrac', [], 'meanRMS_uv', NaN, ...
        'meanExcess_uv', NaN, 'guardMs', P.envCardiacGuardMs), 1, nCh);
    envValidAll = cell(1, nCh);

    fprintf('[step3b] Spike-band RMS envelope: %.1f s bins, R-peak guard +/- %g ms.\n', ...
        P.envBinSec, P.envCardiacGuardMs);

    artPad = round(0.005 * fs);   % +/- 5 ms around each large excursion

    for k = 1:nCh
        xf  = D.filtered(:, k);
        sig = D.sigma(:, k);

        % amplitude guard: exclude excursions above envArtifactSigma * sigma
        % (unguarded cardiac, missed beats, motion) with a small pad
        artifact = abs(xf) > P.envArtifactSigma * sig;
        if any(artifact)
            artifact = movmax(double(artifact), [artPad artPad]) > 0;
        end

        envValid = D.validMask(:, k) & ~guardMask & ~artifact & isfinite(xf);
        envValidAll{k} = envValid;

        nBins    = ceil(N / binN);
        t_c      = nan(nBins, 1);
        rms      = nan(nBins, 1);
        sigFloor = nan(nBins, 1);
        excess   = nan(nBins, 1);
        vfrac    = zeros(nBins, 1);
        for b = 1:nBins
            i0 = (b-1)*binN + 1;
            i1 = min(N, b*binN);
            sel = envValid(i0:i1);
            vfrac(b) = nnz(sel) / (i1 - i0 + 1);
            if any(sel)
                seg  = xf(i0:i1);
                rms(b) = sqrt(mean(seg(sel).^2)) * 1e6;
                ssig = sig(i0:i1);
                sigFloor(b) = median(ssig(sel), 'omitnan') * 1e6;
                % neural power above the noise floor (variances subtract)
                excess(b) = sqrt(max(0, rms(b)^2 - sigFloor(b)^2));
            end
            t_c(b) = ((i0 + i1) / 2 - 1) / fs;
        end

        wellSampled = vfrac >= 0.5 & isfinite(rms);
        meanRMS    = mean(rms(wellSampled), 'omitnan');
        meanExcess = mean(excess(wellSampled), 'omitnan');

        D.envelope(k).t            = t_c;
        D.envelope(k).rms_uv       = rms;
        D.envelope(k).sigmaFloor_uv = sigFloor;
        D.envelope(k).excess_uv    = excess;
        D.envelope(k).validFrac    = vfrac;
        D.envelope(k).meanRMS_uv   = meanRMS;
        D.envelope(k).meanExcess_uv = meanExcess;

        fprintf(['        ch %d (%s): RMS = %.2f uV | noise floor = %.2f uV | ' ...
                 'EXCESS (neural) = %.2f uV  [%d bins]\n'], ...
            ch(k), lab{k}, meanRMS, mean(sigFloor(wellSampled),'omitnan'), ...
            meanExcess, nnz(wellSampled));
    end

    if plotMode
        for k = 1:nCh
            plot_envelope(D, k, envValidAll{k}, smN, P.plotMaxPoints);
        end
    end
end

% ========================================================================
function plot_envelope(D, k, envValid, smN, maxPts)
    fs  = D.fs;
    t   = D.t;
    xf  = D.filtered(:, k);
    env = D.envelope(k);
    valid = D.validMask(:, k);
    lab = D.channelLabels{k};

    fig = figure('Color', 'w', 'Name', sprintf('Step 3b — %s', lab), ...
        'Position', [170 120 1150 760]);
    tl = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Step 3b activity envelope  |  %s  |  excess (neural) RMS = %.2f \\muV', ...
        lab, env.meanExcess_uv), 'Interpreter', 'none');

    % ---- RMS, noise floor, and noise-corrected excess over the recording ----
    ax1 = nexttile;
    plot(env.t, env.rms_uv, 'Color', [0.55 0.55 0.9], 'LineWidth', 0.9, ...
        'DisplayName', 'RMS (total)'); hold on;
    plot(env.t, env.sigmaFloor_uv, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
        'DisplayName', 'noise floor \sigma');
    plot(env.t, env.excess_uv, 'b', 'LineWidth', 1.5, 'DisplayName', 'excess (neural)');
    yline(env.meanExcess_uv, 'k--', 'LineWidth', 1, 'DisplayName', 'mean excess');
    shade_runs(ax1, t, ~valid);
    grid on; xlabel('Time (s)'); ylabel('\muV');
    title('RMS, noise floor, and noise-corrected EXCESS (neural power above noise)');
    legend('show', 'Location', 'northeast');
    xlim([t(1) t(end)]);

    % ---- representative window: cleaned signal + continuous envelope ----
    ax2 = nexttile;
    % pick the window with the highest RMS among WELL-SAMPLED bins (so an
    % artifact-contaminated edge bin is not chosen)
    rmsPick = env.excess_uv;
    rmsPick(env.validFrac < 0.5) = NaN;
    [~, bi] = max(rmsPick);
    if isempty(bi) || ~isfinite(env.t(bi)); c = t(round(numel(t)/2)); else; c = env.t(bi); end
    lo = max(0, c - 2.5); hi = min(t(end), lo + 5);
    i0 = max(1, floor(lo*fs)+1); i1 = min(numel(t), floor(hi*fs)+1);

    xmask = xf;
    xmask(~envValid) = NaN;
    cont = movmean(abs(xmask), smN, 'omitnan') * 1e6;

    plot(t(i0:i1), xf(i0:i1)*1e6, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.4, ...
        'DisplayName', 'cleaned'); hold on;
    plot(t(i0:i1), cont(i0:i1), 'r', 'LineWidth', 1.5, 'DisplayName', 'envelope');
    grid on; xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    title('Most-active window: cleaned signal with continuous |x| envelope');
    legend('show', 'Location', 'northeast');
    xlim([t(i0) t(i1)]);
end

% ========================================================================
function shade_runs(ax, t, mask)
    if ~any(mask); return; end
    yl = ylim(ax);
    d = diff([0; mask(:); 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    for i = 1:min(numel(starts), 300)
        x0 = t(starts(i)); x1 = t(min(ends(i), numel(t)));
        patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.4 0.4], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end
