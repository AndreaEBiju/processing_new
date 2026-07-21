# Implementation Instructions for Claude — Layer 1-2 First-Pass Model
# MATLAB implementation, integrated into the AndreaEBiju/processing_new (GEMS) repo

**Read `model_fitting_plan.md` first for equations, justifications, and limitations
(L1-L11), referenced by ID throughout. This document is execution-only: repo
background, file/variable conventions, and step-by-step build + test instructions.
No re-derivation of equations without flagging back to the user.**

**Do not worry about the L1-L11 limitations for this pass — implement, run against
real data, and report results. Limitations are documented elsewhere; this pass is
about getting a first working fit end-to-end.**

---

## Section 0 — Repo background (established facts — do not re-derive or re-verify these)

This new code is being added to the existing repo `AndreaEBiju/processing_new`
("GEMS physiology batch analysis + plotting"), MATLAB, and MUST reuse its existing
functions rather than duplicating them. The following was confirmed by directly
reading the repo's source files.

### 0.1 Pipeline stage order (`process_dataset.m`)

```
step1a_blank_cardiac  -> step1_bandpass -> step2_noise_sigma -> step3_detect
  -> step3b_envelope -> step4_waveforms -> step5c_modality_test -> step6_spike_report
```
Steps 5 and 5b (clustering-based population separation) are NOT run by default —
manual/optional only. Do not assume clustering-based unit separation exists unless
explicitly re-run.

### 0.2 Data split: baseline vs. stim vs. recovery

- Baseline is stored as its OWN separate file (its own acquisition).
- A second file ("stim_recovery") contains both the stimulation period and the
  post-stim period. `splitStimRecovery.m` splits this into two SAVED files using
  an envelope-threshold detector on a separate stimulation-monitor channel (`vib`):
  - `[baseFileName '_stim.mat']` — contains `x_stim`, `stimMask`, `dig_aligned`, `fs_sig`, `detectInfo`
  - `[baseFileName '_recovery.mat']` — contains `x_recovery`, `recoveryMask`, `dig_aligned`, `fs_sig`, `detectInfo`
- `x_recovery = x(recoveryMask, :)` is a genuine ROW SUBSET (stim samples physically
  removed, not NaN-blanked in place) — it is a clean, continuous, single-segment
  post-stim signal. **This is the correct source for ALL recovery-dynamics analysis
  in this task. Never use `_stim.mat` or unsplit raw files for fitting.**
- Confirmed: `_recovery.mat` contains ONLY post-stim samples (no baseline splice).

### 0.3 Per-file struct `D` (built by `process_dataset.m` and its step-functions)

Key fields accumulated through the pipeline, in order of appearance:

```
D.fs                    % sampling rate (Hz)
D.filtered              % [N x nChan] bandpassed signal
D.validMask             % [N x nChan] logical, false = blanked/invalid sample
D.neuralChannels        % channel indices, e.g. [1 2] for RVN/LVN
D.channelLabels         % cell array of channel name strings, e.g. {'RVN','LVN'}
D.condition             % condition string, e.g. 'M50E100' (parsed elsewhere via parse_me)
D.spikes(k)             % per-channel spike struct (k indexes neuralChannels), from step4:
    .alignedTimes       % spike times (s)
    .alignedCenters     % spike sample indices into D.filtered
    .waveforms          % [nSpikes x nSamples] waveform snippets
    .Vpp_uv             % peak-to-peak amplitude per spike (uV)
    .width_ms           % FWHM per spike (ms)
    .wf_t_ms            % waveform time axis (ms), shared across spikes
D.modality(k)           % from step5c: .features, .p, .hcrit, .verdict (Silverman test)
D.metrics(k)            % from step6 — SEE 0.4 BELOW, this is the primary metrics struct
```

### 0.4 `D.metrics(k)` fields (built by `step6_spike_report.m` — READ THIS FILE if unsure)

```
.label            % channel label string, e.g. 'RVN'
.condition        % condition string
.nSpikes          % spike count
.validDur_s       % valid (non-blanked) duration, seconds
.meanRate_hz      % scalar mean firing rate, blanking-corrected
.medianVpp_uv     % scalar
.medianFWHM_ms    % scalar
.CV, .CV2, .LV    % scalar regularity metrics (gap-clean ISIs only)
.refracViolFrac   % scalar
.fano_T, .fano_F  % Fano-vs-window curve (vectors)
.fanoCanon        % scalar, Fano at P.fanoCanonSec (default 0.1s)
.fanoSlope        % scalar, log-log slope of Fano curve
.fr_t, .fr_hz     % *** TIME-RESOLVED FIRING RATE TRACE *** (vectors, P.frBinSec=1s bins,
                  %     blanking-corrected, spans full valid recovery-file duration)
.cv2_t, .cv2_roll % *** TIME-RESOLVED ROLLING CV2 TRACE *** (vectors, P.cv2WinSec=30s
                  %     windows, gap-clean ISIs only, spans full valid recovery-file duration)
.acg_lag, .acg    % autocorrelogram
.psd_f, .psd_p    % firing-rate power spectrum
.burst            % struct: .hasBursts, .nBursts, .rate_per_min, .meanDur_s, etc.
.meanWaveform, .bandLo, .bandHi, .wf_t_ms   % waveform summary
```

**`.fr_t`/`.fr_hz` and `.cv2_t`/`.cv2_roll` are the time-resolved traces to use for
all recovery-dynamics fitting and plotting in this task. Use these directly —
do not recompute rate/CV2 from `.spikes(k)` from scratch.**

### 0.5 Bulk aggregation across files/animals/conditions (`bulk_compile.m`)

```
rawRows            % struct array, one row per (file, phase, channel), fields:
    .animal        % animal ID string
    .condition     % condition string
    .label         % channel label ('RVN'/'LVN')
    .phase         % 'baseline' or 'recovery'
    .stem          % trial ID (baseline+recovery share stem; used to pair them)
    .mean.(rate|vpp|fwhm|cv2|excess)   % per-trial scalar means
    .dist.(rate|vpp|fwhm|cv2)          % per-trial "distribution" (per-sample) arrays
    .fanoSlope

normRows = bulk_compile(rawRows)   % ONE row per recovery trial x channel:
    .animal, .condition, .label
    .dist.(rate|vpp|fwhm|cv2)      % baseline-normalized: (x - mu_baseline)/mu_baseline
    .scalar.(rate|excess|vpp|fwhm|cv2|fano)   % baseline-normalized scalar effects
```

**OPEN VERIFICATION ITEM (do this first, Section 1 below):** it is not yet confirmed
whether `rawRows(i).dist.rate` / `.dist.cv2` are literally `D.metrics(k).fr_hz` /
`.cv2_roll` carried forward, or something else. Confirm this before writing any
new time-trace code — do not assume.

### 0.6 Statistics already run on this data (`bulk_mixed_models.m`)

A fully saturated (no functional form assumed) linear mixed-effects model already
exists: presence-indicators for m10/m50/m100/e10/e100/e1000 + all 9 M:E interaction
terms, no intercept, `(1|animal)` random effect, fit via `fitlme`, ML, with
likelihood-ratio test for the interaction. Produces per-cell synergy estimates +
SE + p + FDR-corrected p for every metric. **Use this as a validation target**:
the parametric model's fitted M×E surface should approximately reproduce these
per-cell estimates. Do not re-derive this — call `bulk_mixed_models` or reuse its
output structure `R.summary` / `R.cells` directly.

### 0.7 Naming/file conventions to reuse (do not invent new ones)

- Condition strings: `'M<freq1>E<freq2>'` e.g. `'M50E100'`, parsed via regex
  `M(\d+)` / `E(\d+)` (see `parse_me` in `bulk_mixed_models.m`). Reuse this parser.
- Trial pairing: by `stem` (preferred) or fallback to `(animal, condition, channel)`
  match (see `bulk_compile.m`'s `match()` helper). Reuse this logic, do not
  re-implement trial pairing from scratch.
- Channel labels: `'RVN'`, `'LVN'` exactly (case-sensitive in places — check with
  `strcmpi` as the existing code does, not `strcmp`).
- Caching: `loadAllSeriesCached` / `loadAllWindowedCached` share ONE cache file
  (`gemsplots_series_cache.mat`) across `plotWindowedMetrics`, `plotWindowedViolins`,
  `plotTimeTraces` — whichever runs first populates it for the others. If new code
  needs cached series loading, follow this SAME shared-cache convention; do not
  create a separate cache file for the same underlying data.
- Queue building: `buildFileQueue(cacheFile)` returns `state` with `.files`,
  `.animal`, `.condition`, `.phase`, `.groups`. This is the standard entry point
  used by `plotTimeTraces.m`, `bulk_mixed_models.m`'s systemic path, etc. — use it
  rather than writing new file-discovery code.
- Metric spec structs: `struct('label',...,'unitsIn',...,'unitsOut',...,'suffix',...,
  'field',...,'seriesField',...,'timeField',...,'aggregator',...,'channel',...)` —
  see `ms()` helper in `bulk_mixed_models.m` (`hrv_specs`) and `plotTimeTraces.m`
  (`defaultWindowedMetricSpecs`). **Nerve metrics are NOT currently in either
  metric-spec list** — confirmed absent. If nerve time-series need to go through
  `loadAllSeriesCached`, a new metric spec must be added following this exact
  struct shape, pointing at wherever `.fr_hz`/`.cv2_roll` are actually saved to disk
  (verify in Section 1 first).

### 0.8 Plotting style conventions (reuse — do not restyle)

- `pubfig_setup('Theme','light','BaseFontSize',22,'LineWidth',2.0,'MarkerSize',10)`
  called before rendering — call this in any new plotting code for consistent styling.
- Time-trace figures (`plotTimeTraces.m` pattern): mean ± SEM shaded band (theme
  color, `FaceAlpha` 0.20), light individual-trial traces underneath (`Color`
  `[col 0.40]`, `LineWidth` 0.8), bold mean line on top with `DisplayName` showing
  n trials/animals, saved via `exportgraphics` to `.png`/`.svg` plus `savefig` to
  `.fig`. Reuse `saveFigureAllFormats`-style saving for consistency with the rest
  of the repo.
- Report figures (`step6_spike_report.m`'s `plot_report`): `tiledlayout`, `Interpreter`
  set to `'none'` on titles/labels (condition/animal strings often contain
  characters that would otherwise be misinterpreted as LaTeX).
- Heatmaps: `me_heatmap_render.m` is the shared heatmap renderer used by
  `bulk_mixed_models.m` and others — reuse it for any new M×E grid visualization
  rather than writing a new heatmap function.

---

## Section 1 — Required verification before writing new code

- [ ] Open `bulk_compile.m`'s caller (search the repo for what constructs `rawRows`
      before it's passed in — likely `run_pipeline_bulk.m`, `bulk_save_per_group.m`,
      or similar) and confirm exactly how `.dist.rate` / `.dist.cv2` in `rawRows`
      relate to `D.metrics(k).fr_hz` / `.cv2_roll`. Report the exact assignment
      line found.
- [ ] Confirm whether `D.metrics(k).fr_t`/`.fr_hz`/`.cv2_t`/`.cv2_roll` are saved to
      disk anywhere as part of the standard per-file `.mat` output (not just held
      in memory in `D`), and if so, under what filename/field. This determines
      whether new code can load pre-computed traces or must call
      `step6_spike_report` fresh.
- [ ] Confirm current animal/condition/trial coverage: run (or ask the user to run)
      `bulk_mixed_models` on the current dataset and report the printed
      `nObs`/`nAnimals` per metric — use this actual count, not an assumed one, when
      setting fitting expectations.
- [ ] Confirm file-naming pattern for saved per-file outputs (e.g. does a processed
      recovery file save as `<stem>_recovery_processed.mat` or similar, containing
      `D`?) — needed to know how to load existing processed data rather than
      re-running the full pipeline (steps 1a-6) unnecessarily.

Report findings for all four before proceeding to Section 2. If any cannot be
determined from the repo alone, ask the user directly rather than guessing.

---

## Section 2 — New files to create (add to the existing repo, do not fork/duplicate)

```
processing_new/
  model_layer12_equations.m       % pure function library, Section 3
  fit_layer12_stage1_electrical.m % Section 4
  fit_layer12_stage2_interaction.m % Section 5
  plot_layer12_fit_diagnostics.m  % Section 6
  test_layer12_equations.m        % Section 7 (MATLAB unit tests)
  run_layer12_first_pass.m        % Section 8, top-level driver script
```
Follow existing repo naming style (lowerCamelCase or snake_case matching the
nearest analogous existing file — e.g. match `bulk_mixed_models.m`'s style since
this is the closest existing analog). All new functions must accept/return
structs consistent with existing conventions in 0.3-0.7 — do not introduce a
parallel data representation.

---

## Section 3 — `model_layer12_equations.m`: equation building blocks

Implement as a single file with multiple local functions (MATLAB convention: one
primary function + subfunctions, matching the style of `step6_spike_report.m` and
`step5c_modality_test.m`). Implement EXACTLY per `model_fitting_plan.md` Section 2.

**Two corrections from an earlier draft, both load-bearing — implement the
corrected version below, not any version seen in earlier discussion:**
1. `w_c` and `lambda_max_cap` are symbolically non-separable (they multiply together
   everywhere either appears) and are merged into a single parameter `W_c`.
2. The original `f2(u_E) = Phi_c(u_E)` does NOT vanish at `u_E = 0` (`Phi_c(0) = 1`),
   which silently violates the "no standalone mechanical term" requirement and would
   fail the zero-at-`u_E=0` unit test in Section 7. The fix reintroduces the `h(u_E)`
   entrainment gate and gives term_2 its OWN following corner `f_max_c2` (distinct
   from term_1's `f_max_c`) so the two terms stay non-proportional.

```matlab
function out = model_layer12_equations(what, varargin)
%MODEL_LAYER12_EQUATIONS Dispatcher for Layer 1-2 model building blocks.
%   Usage: model_layer12_equations('h', u_E, u_half)
%          model_layer12_equations('Phi_c', u_E, f_max_c)
%          model_layer12_equations('g1', u_M, f_lo, f_hi)
%          model_layer12_equations('g2', u_M, u_M_half)
%          model_layer12_equations('f1', u_E, u_half, f_max_c)
%          model_layer12_equations('f2', u_E, u_half, f_max_c2)
%          model_layer12_equations('r_vagus', u_M, u_E, params)  % params = struct
    switch what
        case 'h';        out = h_fn(varargin{:});
        case 'Phi_c';     out = Phi_c_fn(varargin{:});
        case 'g1';        out = g1_fn(varargin{:});
        case 'g2';        out = g2_fn(varargin{:});
        case 'f1';        out = f1_fn(varargin{:});
        case 'f2';        out = f2_fn(varargin{:});
        case 'r_vagus';   out = r_vagus_fn(varargin{:});
        otherwise; error('model_layer12_equations:unknown','Unknown request: %s', what);
    end
end

function y = h_fn(u_E, u_half)
    y = u_E ./ (u_E + u_half);
end

function y = Phi_c_fn(u_E, f_max_c)
    % f_max_c is the CALLER-SUPPLIED corner -- pass p.f_max_c for term_c/term_1's
    % corner, or p.f_max_c2 for term_2's own, distinct corner. Do not hardcode
    % which corner this uses; it is a generic function of whichever corner is passed.
    y = f_max_c ./ (f_max_c + u_E);
end

function y = g1_fn(u_M, f_lo, f_hi)
    y = (u_M ./ (u_M + f_lo)) .* (f_hi ./ (u_M + f_hi));
end

function y = g2_fn(u_M, u_M_half)
    y = u_M ./ (u_M + u_M_half);
end

function y = f1_fn(u_E, u_half, f_max_c)
    y = h_fn(u_E, u_half) .* Phi_c_fn(u_E, f_max_c);
end

function y = f2_fn(u_E, u_half, f_max_c2)
    % CORRECTED: previously Phi_c_fn(u_E, f_max_c) alone, which does NOT vanish
    % at u_E=0 (Phi_c(0)=1) -- violated the no-standalone-mechanical-term
    % requirement. Now gated by h(u_E) like f1, but with its OWN corner f_max_c2
    % (must differ from f1's f_max_c, else f1 and f2 become proportional and the
    % model loses the ability to shift its preferred u_M with u_E -- see
    % model_fitting_plan.md Section 3).
    y = h_fn(u_E, u_half) .* Phi_c_fn(u_E, f_max_c2);
end

function y = r_vagus_fn(u_M, u_E, p)
    % p is a struct with fields: W_c, u_half, f_max_c, f_lo, f_hi, u_M_half,
    %                            f_max_c2, w_1, w_2
    % NOTE: W_c replaces the old separate w_c/lambda_max_cap (merged -- see above).
    term_c  = p.W_c .* Phi_c_fn(u_E, p.f_max_c) .* h_fn(u_E, p.u_half);
    term_1  = p.w_1 .* f1_fn(u_E, p.u_half, p.f_max_c)  .* g1_fn(u_M, p.f_lo, p.f_hi);
    term_2  = p.w_2 .* f2_fn(u_E, p.u_half, p.f_max_c2) .* g2_fn(u_M, p.u_M_half);
    y = term_c + term_1 + term_2;
    % NO standalone mechanical term (L1) -- every term above now correctly
    % vanishes at u_E=0 by construction. Verify explicitly in tests, Section 7 --
    % this specific property was violated in an earlier draft and must not
    % regress silently.
end
```

Parameter struct field names above are canonical for this task — reuse these exact
names in every downstream file (fitting, plotting, tests). Do not rename. Note there
is no `lambda_cap` dispatcher case anymore — it was folded into `term_c` directly
inside `r_vagus_fn` once `W_c` absorbed `lambda_max_cap`; do not reintroduce it as a
separate named function.

---

## Section 4 — `fit_layer12_stage1_electrical.m`: E-alone sub-model fit (POOLED ACROSS CHANNELS — updated per L12)

**This function's signature and behavior changed from an earlier draft.** An
empirical finding forced this: independent per-channel fits of `(u_half,
f_max_c)` landed in different parameter-space basins with near-identical
residuals. This is now PROVEN (not just observed) to be an exact one-parameter
symmetry of the Stage 1 objective — see `model_fitting_plan.md` L12 for the full
algebraic derivation: swapping `u_half<->f_max_c` and simultaneously rescaling
`W_c -> W_c*(f_max_c/u_half)` leaves the fit EXACTLY unchanged, for any `u_E`.
Synthetic ground-truth validation confirmed this precisely (recovered
`W_c_LVN/W_c_RVN` ratio and `sqrt(u_half*f_max_c)` both nearly exact despite raw
`u_half`/`f_max_c`/`W_c` landing in the "wrong" (swapped) basin). The fix is to
fit `u_half`/`f_max_c` ONCE, POOLED across RVN and LVN simultaneously, while
letting `W_c` remain separate per channel. Do not implement the old
single-channel-only version.

**Consequence for reporting (do not skip this):** because the symmetry is exact,
raw `u_half`/`f_max_c`/`W_c` values (i.e. `thetaE.candidates(i).u_half` etc., for
either candidate `i`) are NOT individually meaningful — only
`sqrt(u_half * f_max_c)` and the ratio `W_c_LVN / W_c_RVN` are identifiable, even
in principle, even with unlimited noiseless data (and, per L13/L15 below, these
invariants are identical for both candidates — that is the whole point). Print
BOTH the raw values (for debugging/reproducing a specific fit) AND these two
invariants explicitly labeled as "the only interpretable Stage 1 quantities" in
the console output.

```matlab
function thetaE = fit_layer12_stage1_electrical(rawRows)
%FIT_LAYER12_STAGE1_ELECTRICAL Fit E-alone sub-model POOLED across RVN and LVN.
%
%   thetaE = fit_layer12_stage1_electrical(rawRows)
%
%   NOTE signature change: no longer takes a channelLabel argument -- this
%   function now ALWAYS fits both channels together. u_half and f_max_c are
%   SHARED across channels; W_c_RVN and W_c_LVN are fit separately. See
%   model_fitting_plan.md L12 for the proof that raw u_half/f_max_c/W_c values
%   are individually non-identifiable (exact swap symmetry) -- only
%   sqrt(u_half*f_max_c) and W_c_LVN/W_c_RVN are meaningful. u_half/f_max_c
%   describe fibre physiology with no expected left/right difference; W_c
%   absorbs per-electrode coupling which plausibly DOES differ.
%
%   Uses RAW (non-normalized) per-trial mean rate from rawRows(i).mean.rate,
%   phase == 'recovery', filtered to E-alone conditions, BOTH channels pooled
%   into one fit. Reuses parse_me-style condition parsing (see
%   bulk_mixed_models.m) to extract u_M, u_E from rawRows(i).condition.
%
%   Objective: nonlinear least squares via lsqnonlin, jointly over both
%   channels:
%     minimize sum_over_BOTH_channels_and_trials[
%       (observed_rate - W_c_{channel}*Phi_c(u_E,f_max_c)*h(u_E,u_half)).^2 ]
%   -- 4 free parameters (u_half, f_max_c, W_c_RVN, W_c_LVN) against
%   3 E-alone conditions x 2 channels = 6 data-generating groups.
%
%   Multiple random restarts (>= 20), report convergence spread.
%   All parameters bounded > 0.
%
%   RETURNS thetaE as a struct with a `candidates` field: a 1x2 struct array,
%   each element holding one of the two Stage-1-equivalent solutions (raw and
%   analytically-swapped -- see below). Do NOT return a flat single-solution
%   struct -- Section 8's driver requires BOTH candidates to run the L13/L15
%   basin-selection step. This is a signature change from earlier drafts.
```

Checklist:
- [ ] Pull raw per-trial rate from `rawRows`, NOT `normRows` (raw values preferred
      per model_fitting_plan.md; percent-change is noisier)
- [ ] Filter to `u_M == 0` (confirm this is how M-alone/E-alone conditions are coded
      in `condition` strings — e.g. does `'E100'` alone appear without an `M` term,
      or as `'M0E100'`? Check actual condition strings in the data before assuming)
- [ ] Include BOTH channels in the SAME fit (do not loop and fit separately —
      that is precisely the procedure that produced the basin-splitting this
      update fixes). Use `strcmpi(rawRows(i).label, 'RVN')` /
      `strcmpi(rawRows(i).label, 'LVN')` to route each trial's contribution to
      the correct `W_c_{channel}` term inside the SAME objective function.
- [ ] Fit `[u_half, f_max_c, W_c_RVN, W_c_LVN]` (4 params total, pooled) via
      `lsqnonlin` — store this raw optimizer output in a local variable (e.g.
      `fitRaw`), NOT directly as the function's return value; it still needs to
      be packaged into the two-candidate struct below.
- [ ] Print convergence diagnostics: final residual, exit flag, parameter values
      across restarts (using `fitRaw`)
- [ ] **Print the two identifiable invariants explicitly, labeled as such (see
      L12 proof) — not just the raw parameter values:**
```matlab
rho   = sqrt(fitRaw.u_half * fitRaw.f_max_c);
ratio = fitRaw.W_c_LVN / fitRaw.W_c_RVN;
fprintf(['Stage 1 raw values (NOT individually meaningful -- exact swap ' ...
         'symmetry, see model_fitting_plan.md L12):\n']);
fprintf('  u_half = %.3f, f_max_c = %.3f, W_c_RVN = %.3f, W_c_LVN = %.3f\n', ...
    fitRaw.u_half, fitRaw.f_max_c, fitRaw.W_c_RVN, fitRaw.W_c_LVN);
fprintf(['Stage 1 IDENTIFIABLE invariants (these ARE meaningful, and are the\n' ...
         'SAME for both candidates below -- that is the point of L12):\n' ...
         '  sqrt(u_half*f_max_c) = %.3f\n  W_c_LVN/W_c_RVN = %.3f\n'], rho, ratio);
```
- [ ] **Required validation, not optional (L12):** ALSO run the original
      independent-per-channel version (2 separate 3-parameter fits) and report
      its total residual alongside the pooled fit's. If pooling increases total
      residual sharply, flag this prominently — it would be evidence against the
      shared-physiology assumption, not something to bury in a log line.
- [ ] **RESOLVED (L13) — this is now a REQUIRED step, not a diagnostic:**
      compute the swapped candidate analytically alongside the raw one (closed-form,
      no re-fit — see formula above). Both candidates get passed forward to
      Section 5/8 for Stage 2 evaluation; do NOT discard the swapped candidate
      after printing it. On real data, RVN and LVN independently selected the
      SAME specific candidate by wide margins (0.139 vs 0.347; 0.306 vs 0.447
      resnorm) — confirming Stage 2's `g1(u_M)` dependence reliably breaks the
      Stage 1 tie. The actual basin SELECTION happens in Section 8's driver
      (not here) once both channels' Stage 2 fits for both candidates exist —
      see Section 8's revised logic.
- [ ] Return a struct containing BOTH candidates (not just one), e.g.:
```matlab
rescaleFactor = fitRaw.f_max_c / fitRaw.u_half;
thetaE.candidates(1) = struct('u_half', fitRaw.u_half, 'f_max_c', fitRaw.f_max_c, ...
    'W_c_RVN', fitRaw.W_c_RVN, 'W_c_LVN', fitRaw.W_c_LVN, 'label', 'raw');
thetaE.candidates(2) = struct('u_half', fitRaw.f_max_c, 'f_max_c', fitRaw.u_half, ...
    'W_c_RVN', fitRaw.W_c_RVN*rescaleFactor, 'W_c_LVN', fitRaw.W_c_LVN*rescaleFactor, 'label', 'swapped');
```
      This is a SIGNATURE CHANGE from earlier drafts — `thetaE` is no longer a
      flat struct with one `u_half`/`f_max_c`, it is a struct with a
      `candidates` array of length 2. Do not name a field `lambda_max_cap` or
      `w_c` (removed, see Section 3).

---

## Section 5 — `fit_layer12_stage2_interaction.m`: interaction fit (STILL PER-CHANNEL — not pooled)

```matlab
function [thetaFull, resnorm2] = fit_layer12_stage2_interaction(rawRows, channelLabel, candidate)
%FIT_LAYER12_STAGE2_INTERACTION Fit [f_lo, f_hi, u_M_half, f_max_c2, w_1, w_2] to
%   the 9 combined M x E conditions for ONE channel, using ONE Stage 1 CANDIDATE
%   (either the raw or swapped member of thetaE.candidates from Section 4 --
%   this function does not know or care which; it is called once per candidate
%   by the Section 8 driver, which does the actual basin selection per L15).
%
%   candidate must have fields: u_half, f_max_c, W_c_RVN, W_c_LVN (i.e. one
%   element of thetaE.candidates(:) from Section 4).
%
%   SIGNATURE CHANGE from earlier drafts: now returns resnorm2 (the fitted
%   Stage 2 residual) as a second output -- required by Section 8's basin
%   selection logic (L15), which compares resnorm2 across candidates and
%   channels. Do not silently drop this second output.
%
%   thetaFull = fit_layer12_stage2_interaction(rawRows, 'RVN', thetaE.candidates(1))
```

Checklist:
- [ ] Filter to combined M×E conditions only (both `u_M > 0` and `u_E > 0`),
      for the requested `channelLabel` only
- [ ] Pull this channel's `W_c` correctly from the candidate struct:
```matlab
switch upper(channelLabel)
    case 'RVN'; W_c_this = candidate.W_c_RVN;
    case 'LVN'; W_c_this = candidate.W_c_LVN;
    otherwise;  error('fit_layer12_stage2_interaction:badChannel', ...
                    'channelLabel must be ''RVN'' or ''LVN'', got ''%s''.', channelLabel);
end
```
      Fix `W_c_this`, `candidate.u_half`, `candidate.f_max_c`.
- [ ] **DO NOT fit `f_lo`/`f_hi` directly — CONFIRMED exact swap symmetry (L16),
      structurally unbreakable (no analog of Section 4's `term_2` tie-break exists
      here). Fit the reparametrized pair instead:**
```matlab
% Optimizer fits rho1 (peak location) and kappa1 (log-ratio, sign-constrained)
% directly, NOT f_lo/f_hi. This is required, not optional -- without it, restarts
% land in two mirror-image families at equivalent resnorm (already observed:
% f_lo [0.3557, 37.68], f_hi [0.005314, 35.11], w_1 [6.718, 4.533e+04] on real
% data -- that spread IS this symmetry, not a convergence failure).
%
% params vector for lsqnonlin: [rho1, kappa1, u_M_half, f_max_c2, w_1_transformed, w_2]
%   f_lo  = rho1 * exp(-kappa1)
%   f_hi  = rho1 * exp(+kappa1)
%   bound kappa1 >= 0 (report sign as a FIXED CONVENTION, not a fitted finding --
%   this means f_hi >= f_lo always, by construction, not by data)
%   w_1 = w_1_transformed  -- NOTE: since f_hi/f_lo = exp(2*kappa1) is now fully
%   determined by kappa1 alone (no separate f_hi/f_lo ambiguity feeding into it),
%   w_1 no longer needs its own compensating transformation -- fit it directly.

function y = g1_reparam_fn(u_M, rho1, kappa1)
    f_lo = rho1 * exp(-kappa1);
    f_hi = rho1 * exp(kappa1);
    y = model_layer12_equations('g1', u_M, f_lo, f_hi);
end
```
      After fitting, recover `f_lo`/`f_hi` for reporting via the formulas above --
      but treat them as DERIVED display quantities, not independently meaningful
      fitted parameters. Report `rho1` and `kappa1` (or equivalently `f_hi/f_lo`)
      as the actual identifiable outputs, the same way Section 4 reports
      `sqrt(u_half*f_max_c)` and the `W_c` ratio rather than raw `u_half`/`f_max_c`.
- [ ] **Apply the IDENTICAL reparametrization to `(u_half, f_max_c)` in Section 4
      for consistency** — these are the same structural problem (L12 and L16 are
      the same phenomenon on two different parameter pairs). Section 4's current
      draft still fits `u_half`/`f_max_c` directly and relies on Stage 2 to break
      the tie (L13) rather than removing the ambiguity by construction; both are
      valid fixes for L12 specifically (since L13's tie-break works there), but
      reparametrizing Section 4 the same way as Section 5 avoids relying on two
      different resolution mechanisms for structurally identical problems. Decide
      once and document which approach was used for `(u_half, f_max_c)` — the
      tie-break (already implemented per Section 4/8) is NOT wrong, just worth
      flagging as a design choice rather than defaulting silently.
- [ ] **`f_max_c2` MUST be constrained distinct from `f_max_c`** (e.g. enforce
      `f_max_c2 < f_max_c` via bounds, or at minimum verify post-fit that they did
      not converge to the same value) — if they converge to the same value, `f1`
      and `f2` become proportional and the model loses the moving-optimum property
      that motivated having two mechanical-coupling terms in the first place (see
      model_fitting_plan.md Section 3). Report this check explicitly. This check
      is UNRELATED to the L16 reparametrization above — it is about `f_max_c`
      (term_1/term_c's shared corner) vs. `f_max_c2` (term_2's own corner), not
      about `f_lo` vs. `f_hi`. Do not conflate the two checks.
- [ ] **New (L17): flag `(u_M_half, w_2)` asymptotic ridge risk — a DIFFERENT
      kind of issue from L16, do not fix the same way.** This is a data-coverage
      problem (only 3 discrete `u_M` points tested: 10, 50, 100 Hz), not an exact
      algebraic symmetry — no reparametrization removes it, only a post-fit flag:
```matlab
maxUM_tested = 100;  % update if the tested u_M range changes
if thetaFull.u_M_half > 5 * maxUM_tested
    kappa2 = thetaFull.w_2 / thetaFull.u_M_half;
    warning('fit_layer12_stage2_interaction:asymptoticRidge', ...
        ['[%s] u_M_half = %.1f is > 5x max tested u_M (%.0f) -- L17 asymptotic ' ...
         'ridge risk. u_M_half and w_2 individually unreliable in this regime; ' ...
         'report kappa2 = w_2/u_M_half = %.4g as the more robust quantity ' ...
         'alongside (not instead of) the raw values.'], channelLabel, ...
        thetaFull.u_M_half, maxUM_tested, kappa2);
    thetaFull.asymptoticRidgeFlag = true;
    thetaFull.kappa2 = kappa2;
else
    thetaFull.asymptoticRidgeFlag = false;
end
```
      **Do NOT add an upper bound on `u_M_half` in the optimizer to suppress
      this** — a genuinely large fitted `u_M_half` (i.e. `g2` really is close to
      linear across the whole tested range) is a real finding, not an error;
      constraining it away would hide rather than fix the limitation. Only
      flag/monitor. This connects to the already-documented L4 limitation (coarse
      `u_M` sweep) — the actual fix is the finer sweep L4 already recommends, not
      anything inside this function.
- [ ] Same optimizer/restart/bounds conventions as Section 4
- [ ] Return `resnorm2` as the second output — the best (lowest) residual sum of
      squares found across restarts. This is what Section 8 uses to select the
      winning candidate (L15) — it must be the actual fitted objective value, not
      an approximation.
- [ ] **Required check — print explicitly, pass/fail, do not just report parameter
      values:**
```matlab
% thetaFull must have a scalar W_c field for model_layer12_equations to use --
% build this by copying W_c_this into thetaFull.W_c before calling r_vagus:
thetaFull.W_c = W_c_this;   % single field name expected by model_layer12_equations
pred_M50E100  = model_layer12_equations('r_vagus', 50, 100, thetaFull);
pred_M100E100 = model_layer12_equations('r_vagus', 100, 100, thetaFull);
pred_M100E10  = model_layer12_equations('r_vagus', 100, 10, thetaFull);
pred_M50E10   = model_layer12_equations('r_vagus', 50, 10, thetaFull);
fprintf('[%s] Check 1 (M50xE100 > M100xE100): %s (%.3f vs %.3f)\n', channelLabel, ...
    string(pred_M50E100 > pred_M100E100), pred_M50E100, pred_M100E100);
fprintf('[%s] Check 2 (M100xE10 > M50xE10):   %s (%.3f vs %.3f)\n', channelLabel, ...
    string(pred_M100E10 > pred_M50E10), pred_M100E10, pred_M50E10);
```
      If either check is FALSE, report this prominently — do not bury it — and
      still return the fit (do not silently retry with different starting points
      to force a pass). Run and report this SEPARATELY for RVN and LVN — a pass
      on one channel does not imply a pass on the other.
- [ ] Cross-check final fitted M×E surface against `bulk_mixed_models`'s per-cell
      synergy estimates (Section 0.6) — plot predicted vs. the saturated model's
      cell means side by side (see Section 6).

---

## Section 6 — `plot_layer12_fit_diagnostics.m`: diagnostics


Reuse `pubfig_setup` and repo plotting conventions (Section 0.8). Produce, for
each channel (RVN, LVN) separately:

- [ ] Observed vs. predicted scatter (one point per trial), Stage 1 (E-alone) and
      Stage 2 (interaction) on separate panels
- [ ] Predicted M×E heatmap (3x3 grid, M in {10,50,100} x E in {10,100,1000}) using
      `me_heatmap_render.m` for visual consistency with existing heatmaps, placed
      side by side with the actual `bulk_mixed_models` per-cell synergy heatmap for
      the same metric/channel
- [ ] Residuals vs. `u_M` and residuals vs. `u_E`, separately, as scatter plots
- [ ] Print the Section 5 pass/fail check result as text overlay or console output
      alongside the figure

---

## Section 7 — `test_layer12_equations.m`: MATLAB unit tests

Use MATLAB's built-in unit test framework (`matlab.unittest.TestCase`), matching
professional repo conventions rather than ad hoc script assertions:

```matlab
classdef test_layer12_equations < matlab.unittest.TestCase
    methods (Test)
        function test_r_vagus_zero_at_uE_zero(tc)
            % REGRESSION TEST: an earlier draft of f2 used Phi_c(u_E) alone,
            % which does NOT vanish at u_E=0 (Phi_c(0)=1) and silently broke
            % this exact property. This test must pass on the CURRENT f2
            % (gated by h(u_E) with its own corner f_max_c2) -- if it starts
            % failing again, f2 has regressed to the old broken form.
            p = default_test_params();
            uM = [10 50 100];
            for i = 1:numel(uM)
                y = model_layer12_equations('r_vagus', uM(i), 0, p);
                tc.verifyEqual(y, 0, 'AbsTol', 1e-9);
            end
        end

        function test_f2_zero_at_uE_zero_specifically(tc)
            % Isolates f2 itself (not just the full r_vagus sum) to make the
            % regression this guards against unambiguous if it ever reappears.
            y = model_layer12_equations('f2', 0, 50, 100);
            tc.verifyEqual(y, 0, 'AbsTol', 1e-9);
        end

        function test_f1_f2_not_proportional(tc)
            % f_max_c and f_max_c2 must differ, else f1 and f2 collapse to
            % rescaled copies of each other and the model loses the
            % moving-optimum property (model_fitting_plan.md Section 3).
            % Checks that the ratio f1/f2 is NOT constant across u_E.
            uE = [10 50 100 500 1000];
            u_half = 50; f_max_c = 200; f_max_c2 = 60;  % deliberately distinct
            f1v = model_layer12_equations('f1', uE, u_half, f_max_c);
            f2v = model_layer12_equations('f2', uE, u_half, f_max_c2);
            ratio = f1v ./ f2v;
            tc.verifyTrue(range(ratio) > 1e-6, ...
                'f1/f2 ratio is constant -- f_max_c and f_max_c2 must differ.');
        end

        function test_h_monotonic_increasing(tc)
            uE = linspace(0.1, 2000, 500);
            y = model_layer12_equations('h', uE, 100);
            tc.verifyTrue(all(diff(y) >= 0));
            tc.verifyTrue(all(y >= 0 & y < 1));
        end

        function test_Phi_c_monotonic_decreasing(tc)
            uE = linspace(0.1, 2000, 500);
            y = model_layer12_equations('Phi_c', uE, 100);
            tc.verifyTrue(all(diff(y) <= 0));
            tc.verifyTrue(all(y > 0 & y <= 1));
            % NOTE: Phi_c(0) = 1, not 0 -- this is correct and expected for
            % Phi_c in isolation. It is f2 (gated by h) that must vanish at
            % u_E=0, not Phi_c itself. Do not "fix" this test to expect 0.
        end

        function test_g1_has_single_interior_peak(tc)
            uM = linspace(0.1, 500, 1000);
            y = model_layer12_equations('g1', uM, 20, 80);
            [~, idx] = max(y);
            tc.verifyTrue(idx > 1 && idx < numel(uM));  % interior, not monotonic
            tc.verifyTrue(all(diff(y(1:idx)) >= -1e-9));   % rising before peak
            tc.verifyTrue(all(diff(y(idx:end)) <= 1e-9));  % falling after peak
        end

        function test_g2_monotonic_increasing_bounded(tc)
            uM = linspace(0.1, 500, 500);
            y = model_layer12_equations('g2', uM, 50);
            tc.verifyTrue(all(diff(y) >= 0));
            tc.verifyTrue(all(y >= 0 & y < 1));
        end

        function test_no_division_errors_at_zero(tc)
            p = default_test_params();
            tc.verifyWarningFree(@() model_layer12_equations('r_vagus', 0, 0, p));
            tc.verifyWarningFree(@() model_layer12_equations('r_vagus', 0, 100, p));
        end
    end
end

function p = default_test_params()
    p = struct('W_c',1,'u_half',50,'f_max_c',200, ...
                'f_lo',20,'f_hi',80,'u_M_half',50,'f_max_c2',60, ...
                'w_1',1,'w_2',1);
    % f_max_c2 deliberately != f_max_c (200 vs 60) -- see
    % test_f1_f2_not_proportional above for why this must hold.
end
```

Checklist:
- [ ] Run `runtests('test_layer12_equations')` and confirm ALL PASS before proceeding
      to Section 4 fitting. Do not fit against real data with failing unit tests.

---

## Section 8 — `run_layer12_first_pass.m`: top-level driver

**MAJOR REVISION per L13/L15 — this section now implements required basin
selection, not just sequential fitting.** Earlier drafts ran Stage 1 once, then
Stage 2 once per channel, using whatever Stage 1 converged to. That is no longer
correct: Stage 1's `(u_half, f_max_c)` is provably non-unique (L12), and Stage 2
reliably picks a winner between the two candidates (L13) — but only if the driver
actually evaluates both and selects, rather than passing through Stage 1's raw
output unexamined.

```matlab
function results = run_layer12_first_pass(rawRows)
%RUN_LAYER12_FIRST_PASS Top-level driver for the Layer 1-2 first-pass model fit.
%   results = run_layer12_first_pass(rawRows)   % rawRows from existing pipeline
%
%   Runs unit tests, Stage 1 (pooled, both candidates), Stage 2 for BOTH
%   candidates x BOTH channels, SELECTS the Stage-2-winning candidate (L15),
%   then diagnostics using ONLY the winning candidate's fits, then cross-check
%   against bulk_mixed_models. Returns a struct with the selected parameters,
%   full audit trail of both candidates' fits, and diagnostic figure handles.
```

Checklist:
- [ ] Step 1: `runtests('test_layer12_equations')` — abort with clear error
      message if any test fails
- [ ] Step 2: `thetaE = fit_layer12_stage1_electrical(rawRows)` — called ONCE,
      POOLED across both channels. Per Section 4's signature change, `thetaE`
      now contains `thetaE.candidates(1)` (raw) and `thetaE.candidates(2)`
      (swapped) — NOT a single flat `u_half`/`f_max_c`.
- [ ] Step 3 — **run Stage 2 for BOTH candidates, BOTH channels (4 fits total):**
```matlab
channels = {'RVN','LVN'};
resnorm2 = nan(2,2);      % rows = candidate index, cols = channel index
thetaFullAll = cell(2,2); % same indexing
for ci = 1:2
    for k = 1:numel(channels)
        [thetaFullAll{ci,k}, resnorm2(ci,k)] = fit_layer12_stage2_interaction( ...
            rawRows, channels{k}, thetaE.candidates(ci));
    end
end
```
- [ ] Step 4 — **select the winning candidate (L15), by TOTAL resnorm summed
      across channels (joint/shared-parameter decision):**
```matlab
totalResnorm = sum(resnorm2, 2);   % [total for candidate 1; total for candidate 2]
[bestTotal, winnerIdx] = min(totalResnorm);
otherIdx = 3 - winnerIdx;          % the non-winning candidate's index
marginPct = 100 * (totalResnorm(otherIdx) - bestTotal) / bestTotal;

MIN_MARGIN_PCT = 10;   % ask user to confirm/adjust this threshold before relying
                        % on auto-selection for a new dataset
if marginPct < MIN_MARGIN_PCT
    warning('run_layer12_first_pass:ambiguousBasin', ...
        ['Stage 1 basin selection margin only %.1f%% (< %.0f%% threshold) -- ' ...
         'AMBIGUOUS, review manually before trusting auto-selected results.'], ...
        marginPct, MIN_MARGIN_PCT);
    results.basinAmbiguous = true;
else
    results.basinAmbiguous = false;
end
```
- [ ] Step 5 — **required per-channel agreement check, every run (not just this
      one-off manual comparison) — this is a NEW automated safety check, not
      previously in any draft:**
```matlab
[~, winnerPerChannel] = min(resnorm2, [], 1);   % winner index per channel
if ~all(winnerPerChannel == winnerPerChannel(1))
    warning('run_layer12_first_pass:channelDisagreement', ...
        ['Channels DISAGREE on which Stage 1 candidate is correct -- this is ' ...
         'evidence AGAINST the L12 pooling assumption (shared u_half/f_max_c ' ...
         'across channels), not just a selection nuance. Report prominently; ' ...
         'do not silently default to the total''s answer.']);
    results.channelDisagreement = true;
else
    results.channelDisagreement = false;
end
```
      On the real data seen so far, RVN and LVN agreed (both preferred the same
      candidate) — this check exists for future datasets, where agreement is not
      guaranteed just because it held once.
- [ ] Step 6 — promote the winning candidate. **Do NOT re-fit a third time** —
      the winning candidate's Stage 2 fits from Step 3 already exist:
```matlab
results.thetaE_selected = thetaE.candidates(winnerIdx);
results.thetaE_rejected = thetaE.candidates(otherIdx);
results.selectionMarginPct = marginPct;
for k = 1:numel(channels)
    results.perChannelDetail.(channels{k}).thetaFull = thetaFullAll{winnerIdx,k};
    results.perChannelDetail.(channels{k}).resnorm2  = resnorm2(winnerIdx,k);
end
```
- [ ] Step 7: for each channel, run `plot_layer12_fit_diagnostics` using ONLY
      `results.perChannelDetail.(channel).thetaFull` (the winning candidate's
      fit) — diagnostics/plots must never be built from the rejected candidate
- [ ] Step 8: run (or call) `bulk_mixed_models` on the same data if not already
      available, for cross-check plots in Section 6
- [ ] Step 9: print a final summary table (channel x parameter x fitted value) to
      console and save to a `.mat`/`.csv` file following existing repo output
      conventions. Report `u_half`/`f_max_c` ONCE (shared, from
      `thetaE_selected`), not duplicated per channel. **Also save the rejected
      candidate and both channels' resnorms for BOTH candidates as an audit
      trail** — do not discard this once selection is made.
- [ ] Step 10: report the Section 5 pass/fail pattern check result for BOTH
      channels, computed from the SELECTED candidate's fits — prominently, this
      is the single most important result of this first pass
- [ ] Step 11: report the Section 4 pooled-vs-independent Stage 1 residual
      comparison (L12 validation) prominently as well
- [ ] Step 12: report the basin selection margin (Step 4) and any warnings from
      Steps 4-5 prominently in the final printed summary — this is new,
      required reporting, not optional diagnostic output

---

## Section 9 — Final report back to user

After running against real data, report:
- [ ] Fitted parameter values, both channels, with convergence diagnostics —
      report `u_half`/`f_max_c` as SHARED values (one pair, not two) — **from
      the SELECTED candidate only (`results.thetaE_selected`)**
- [ ] **Basin selection outcome (L13/L15) — new, required:** which candidate
      was selected, the margin (`results.selectionMarginPct`), whether flagged
      ambiguous (`results.basinAmbiguous`), and whether channels agreed
      (`results.channelDisagreement`). Report the rejected candidate's values
      too, for audit purposes — do not report only the winner as if the other
      never existed.
- [ ] Stage 1 pooled-vs-independent residual comparison (L12) — state explicitly
      whether pooling was defensible for this dataset or not
- [ ] **`(u_M_half, w_2)` asymptotic ridge flag (L17), both channels** — if
      `asymptoticRidgeFlag` is true for either channel, report it prominently and
      show `kappa2` alongside the raw values; note this is a data-coverage
      limitation (L4), not an error, and does not block interpreting the rest of
      the fit — but raw `u_M_half`/`w_2` should not be reported as precise when
      flagged
- [ ] Pass/fail on the M50xE100 / M100xE10 pattern check, both channels
      (reported separately — a pass on one channel does not imply the other),
      **computed from the selected candidate's fits**
- [ ] Side-by-side comparison plots: this model's predicted M×E surface vs.
      `bulk_mixed_models`'s saturated per-cell synergy estimates
- [ ] Any condition/channel combinations with insufficient trials to fit (report
      counts, do not silently exclude)
- [ ] Any unexpected data format issues encountered (e.g. if Section 1's open
      verification items turned out different from documented above)
