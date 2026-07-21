close all force
% RUN_PIPELINE_SINGLE  Interactive single-dataset driver for the VENG pipeline.
%
% Run this section by section (Ctrl+Enter on each %% block) so you can inspect
% the result of every step before moving on. Each step is its own function and
% takes a plotMode flag; in single-file mode plotMode = true so every step
% draws its diagnostic plots. In bulk mode (later) the same functions are
% called with plotMode = false.
%
% Assumes stim/recovery splitting is already done upstream.

%% Step 0 — parameters + load one dataset (GUI)
P = pipeline_params();
plotMode = true;                 % single-file mode: show all diagnostics
figuresAtEnd = false;             % create figures hidden, reveal + save them all at the end
runExploratory = false;          % step5/5b (t-SNE/DBSCAN, over-cluster) — off by default:
                                 % they over-segment a continuum and are NOT defensible; the
                                 % modality test (step 5c) is the verdict. Set true to explore.
runSweep       = true;           % step5d threshold sweep + noise surrogate (quick check):
                                 % "are there sub-/supra-threshold populations?" — needs only
                                 % through step2. Set false to skip.
if figuresAtEnd; set(0, 'DefaultFigureVisible', 'off'); end  % (restored in the final section)

D = step0_load_data(P);
if isempty(D)
    fprintf('Cancelled at load.\n');
    return;
end

%% Step 1a — blank cardiac (QRS) windows in the RAW signal (BEFORE filtering)
% NaN-blanks +/- P.cardiacRemoveWinMs around each R-peak so the sharp QRS is
% gone before the bandpass (no filter ringing) and the heartbeat is excluded
% everywhere downstream. Replaces the old template subtraction (1b) + removal (9).
D = step1a_blank_cardiac(D, P, plotMode);

%% Step 1 — bandpass filter (100-5000 Hz), with before/after FFT + trace
% Interpolates across the cardiac blanks to filter, then restores the NaNs.
D = step1_bandpass(D, P, plotMode);

% Inspect:
%   D.filtered   -> bandpassed neural channels (samples x nNeural)
%   D.bandInfo   -> corners actually used
% One figure per neural channel shows FFT and full trace before vs after.

%% Step 2 — validity mask + time-varying noise sigma (5-s sliding)
% Flags the zero/flat recording-error region (and NaN blanks, removed
% segments, edges) so they are excluded from sigma and from later steps.
D = step2_noise_sigma(D, P, plotMode);

% Inspect:
%   D.validMask  -> true where data is usable (samples x nNeural)
%   D.sigma      -> per-sample noise sigma track (volts)
%   D.noiseInfo  -> median sigma per channel + params
% Per channel: top plot = filtered signal with +/- threshold envelope;
% bottom = sigma vs time; invalid regions shaded. The dead tail after the
% recording error should be shaded and excluded from the sigma track.

%% Step 5d (quick check) — threshold sweep + phase-randomized noise surrogate
% Exploratory: sweeps detection 2..8 sigma, detects an identical noise surrogate,
% and tests modality of TIMESCALE features (trough-to-peak, spectral centroid)
% vs the surrogate floor. Answers "are there sub-/supra-threshold populations?"
% Needs only D.filtered + D.validMask (i.e. through step2). Not in the production
% flow; run it here while exploring a single file.
if runSweep
    sweep = step5d_threshold_sweep(D, P, plotMode);   %#ok<NASGU>
end

%% Step 3 — spike detection (threshold crossings on the sigma track)
D = step3_detect(D, P, plotMode);

% Inspect:
%   D.spikes(k).times / .centers   -> detected events for neural channel k
%   D.spikes(k).peakAmp_uv         -> per-spike amplitude
%   D.spikes(k).rate_hz            -> events per valid second
% Per channel: representative window with detected events, and the
% amplitude histogram. A hard wall at the threshold line = low-amplitude
% C-fibre events are being missed (then a template-recovery pass is worth
% adding). A natural taper past the line = threshold is fine.

%% Step 3b — spike-band RMS envelope (PRIMARY low-SNR activity readout)
% Threshold-free activity measure: RMS of the cleaned signal in time bins,
% excluding a guard window around R-peaks. This is the headline "how much is
% the nerve firing" metric for cross-condition comparison.
D = step3b_envelope(D, P, plotMode);

% Inspect:
%   D.envelope(k).rms_uv     -> activity trace (per bin) for channel k
%   D.envelope(k).meanRMS_uv -> scalar activity measure for this condition

%% Step 4 — waveform extraction + alignment (detectable events)
D = step4_waveforms(D, P, plotMode);

% Inspect:
%   D.spikes(k).waveforms     -> aligned waveform matrix for channel k
%   D.spikes(k).meanWaveform  -> template shape
%   D.spikes(k).Vpp_uv / .width_ms -> per-spike shape features
% Per channel: waveform overlay with mean +/- std, and amplitude-vs-width
% scatter (preview of whether distinct waveform populations exist).

%% (optional) peri-R contamination QC figure — NOT in the production flow.
% Cardiac is already handled by step1a (blank-before-filter). To regenerate the
% "what fraction of detected events is heartbeat-locked" figure for the methods
% slide, run a separate pass with step1a skipped and call:
%   step9_cardiac_remove(D, P, false);   % tag-only (P.cardiacRemove=false)

%% Step 5 / 5b — EXPLORATORY clustering (off by default; over-segments a continuum)
% These are kept for exploration only. We established the data is a continuum
% (Step 5c is the defensible verdict), and t-SNE/DBSCAN / over-clustering invent
% spurious cluster counts here. Enable by setting runExploratory = true above.
if runExploratory
    D = step5_sort_check(D, P, plotMode);
    D = step5b_verify_populations(D, P, plotMode);
end

%% Step 5c — formal modality test (continuum vs discrete) on the features
% Run on baseline AND recovery. p<0.05 on a feature = genuine multimodality
% (discrete structure); all unimodal = continuum (no discrete populations).
D = step5c_modality_test(D, P, plotMode);

%% Step 6 — single-population spike-train report (the per-file characterisation)
% Firing rate, waveform+band, Vpp over time, log-ISI + bursts, time x ISI
% heatmap, return map, Fano-vs-window, autocorrelogram, firing-rate power
% spectrum, amplitude-vs-ISI, rolling CV2, scalar summary. All scalars stored
% in D.metrics(k) for the bulk stage.
D = step6_spike_report(D, P, plotMode);

%% Step 7 — gap-aware DFA (fractal scaling) on the blanking-corrected firing rate
D = step7_dfa_report(D, P, plotMode);

%% Save this condition's summary (run once per file, naming the condition)
% Process each file (baseline, recovery, ...) through the steps above, then:
pipeline_save_summary(D, P, 'baseline');   % <- change name per file: 'baseline','recovery',...

%% Compare conditions (after you've saved 2+ summaries)
compare_conditions();      % activity: excess RMS + rate time courses + means
compare_distributions();   % fiber question: FWHM/Vpp distributions + KS test

%% Reveal + save ALL figures (PNG, FIG, SVG)
% Restores visibility and writes every figure opened during this run to a
% "figures" folder next to the neural file, in all three formats.
set(0, 'DefaultFigureVisible', 'on');
set(findall(0, 'Type', 'figure'), 'Visible', 'on');
[nfDir, nfBase] = fileparts(D.neuralFile);
save_all_figures(fullfile(nfDir, 'figures'), nfBase, {'png','fig','svg'});

%% (next steps will be appended here: physiology / comparison, ...)
