# DFA Build Plan — Gap-Aware Fractal Scaling for GEMS Pipeline

For a Claude Code session with access to `github.com/AndreaEBiju/processing_new`.
This document specifies the algorithm and
exact build steps for two cases -- cardiac (RR-interval-indexed) and nerve
(1 Hz-binned firing rate) -- both discontinuous due to artifact removal. The
core engine (`dfaGapAware.m`) is delivered complete; what remains is two thin
domain-specific wrappers plus wiring into the existing pipeline outputs.

## 0. Delivered

**`dfaGapAware.m`** -- complete, generic gap-aware DFA engine. Do not
reimplement the core algorithm; call this function from both wrappers.

Signature: `out = dfaGapAware(x, validMask, P)`

- `x`: any single index-domain series (RR intervals or firing-rate bins).
- `validMask`: same-length logical, true = usable. OR pass
  `P.dfaRunsOverride` (Nx2 index list) instead when gap structure isn't
  naturally per-sample (see Cardiac spec below) -- in that case pass
  `validMask = true(size(x))` and let `runs` do the real work.
- `P`: `pipeline_params()` struct; reads `P.dfa*` fields (see §3).

Returns `out.alpha1`, `out.alpha2`, `out.alphaFull` (scaling exponents),
`out.R2_1`/`out.R2_2`/`out.R2_full` (fit quality), `out.nCross` (crossover
scale), `out.nWindows` (pooled window count per retained scale -- QC, not a
physiological metric), `out.localAlpha` (scale-dependent slope curve),
`out.excludedScales`, `out.runs`.

Method: pooled valid-window F²(n) accumulation across contiguous clean
runs, forward + reverse tiling within each run, no interpolation across
gaps, scales dropped when too few pooled windows survive
(`P.dfaMinWindowsPerScale`).

## 1. Cardiac case -- build spec

**Inputs available:** `RR_intervals` (vector, seconds), `RR_times` (vector,
seconds, onset time of each RR interval), both already computed by
`HR_BR_HRVAnalysis_new.m` with gap-spanning intervals pre-removed by
`computeValidRRIntervals`.

**Key fact driving the algorithm:** `RR_intervals`/`RR_times` already
exclude any interval that itself spans a gap, but consecutive *surviving*
entries can still have a real-time jump between them wherever a beat was
dropped. A `validMask` over `RR_intervals` cannot represent this -- every
individual entry is legitimately valid; it's the *transition* between two
valid entries that can be broken. Therefore: build a **runs list**, not a
`validMask`.

**Build `dfaRR_gapAware.m`:**

```matlab
function out = dfaRR_gapAware(RR_intervals, RR_times, P)
% Build a runs list: consecutive RR_intervals entries i, i+1 belong to the
% same run iff no beat was dropped between them.
if nargin < 3 || isempty(P); P = pipeline_params(); end
tol = P.dfaBeatGapTolSec; % ~1e-3 s; floating-point slop only, NOT physiological

n = numel(RR_intervals);
gapAfter = [RR_times(2:end) - RR_times(1:end-1) - RR_intervals(1:end-1); 0];
breakAfter = abs(gapAfter) > tol; % true at index i means run breaks between i, i+1

starts = [1; find(breakAfter(1:end-1)) + 1];
ends   = [find(breakAfter(1:end-1)); n];
runs = [starts, ends];

P.dfaRunsOverride = runs;
out = dfaGapAware(RR_intervals, true(n,1), P);
end
```

Un-comment the DFA call in `HR_BR_HRVAnalysis_new.m`:

```matlab
% replace:
% DFA intentionally omitted -- not appropriate for interrupted/non-stationary recordings
% [alpha1, alpha2, nVals, F] = dfaRR(RR_intervals);
% with:
dfaOut = dfaRR_gapAware(RR_intervals, RR_times, pipeline_params());
```

Add `dfaOut.alpha1`, `dfaOut.alpha2`, `dfaOut.alphaFull`, `dfaOut.R2_1`,
`dfaOut.R2_2`, `dfaOut.nCross`, `dfaOut.nWindows`, `dfaOut.excludedScales`
to the `out` struct and to `save(hrvFile, ...)`'s variable list, grouped
near `hrv`/`rmssd`/`sd1`/`sd2`. Suggested saved names: `dfa_alpha1`,
`dfa_alpha2`, `dfa_alphaFull`, `dfa_R2_1`, `dfa_R2_2`, `dfa_nCross`,
`dfa_nWindows`, `dfa_excludedScales`.

## 2. Nerve case -- build spec

**Inputs available:** `D.validMask` (samples x nCh, from
`step2_noise_sigma.m`), `D.removedSegmentIdx` (Nx2 sample indices, shared
across channels), firing rate via `[t, fr] = firing_rate(cen, valid, N, fs,
binSec)` (confirmed signature; source file and full per-bin output set NOT
yet confirmed in this repo -- resolve via step 2.0 below before writing
anything else).

**Key fact driving the algorithm:** heartbeat blanking (30 ms per beat,
already dead-time-corrected in `firing_rate()`'s denominator) is invisible
to DFA and needs no handling. Motion-artifact blanking is the only
DFA-relevant gap source for this channel, and unlike the cardiac case, the
firing-rate series IS naturally bin-indexed with one value per bin --
so a per-bin `validMask` (not a runs override) is the right representation
here, straightforwardly built from whichever bins got excluded by the
already-implemented valid-duration floor.

### Step 2.0 -- confirm before writing code (required first step)

```
grep -n "function.*firing_rate" step*.m
grep -rn "fr_hz\|fr_t\|fr_valid\|validFrac" step6*.m
```

Resolve:
1. Exact file containing `firing_rate()` (assumed `step6_spike_report.m`;
   confirm, don't assume).
2. Whether `firing_rate()` returns a third output (per-bin valid-duration
   or valid-fraction), or whether the caller separately computes the
   floor-exclusion using something like `validSecBin = ... * binSec`
   (this exact pattern -- `env.validFrac(:) * binSec` -- already exists
   for the Step 3b envelope in `pipeline_save_summary.m`; check whether
   Step 6 mirrors it).
3. Exact field name/struct (`M.fr_t`/`M.fr_hz` or similar) these live on,
   and whether that struct is returned from Step 6 or written onto `D`.

### Step 2.1 -- build the per-bin validMask

Once §2.0 gives you the real field names, build:

```matlab
binValidMask = <per-bin valid-fraction or keep-vector> >= <existing floor threshold>;
```

using whatever threshold the existing floor-exclusion logic already uses
(do not invent a new threshold -- reuse the one already in place, per the
user's confirmation that this exclusion is already implemented).

### Step 2.2 -- build `step7_dfa_report.m`

Follow the `D = stepN(D, P, plotMode)` convention from steps 1-6.

```matlab
function D = step7_dfa_report(D, P, plotMode)
if nargin < 2 || isempty(P); P = pipeline_params(); end
if nargin < 3 || isempty(plotMode); plotMode = true; end

% D.fr_t, D.fr_hz, D.fr_validMask (or equivalent) must exist from Step 6.
D.dfa = dfaGapAware(D.fr_hz, D.fr_validMask, P);

if plotMode
    plot_dfa(D.fr_t, D.dfa); % see step 2.3
end
end
```

### Step 2.3 -- diagnostic plot

Add a `plot_dfa` subfunction inside `step7_dfa_report.m` matching
`step2_noise_sigma.m`'s `plot_sigma`/`shade_runs` conventions: log-log F(n)
vs n, excluded scales marked distinctly (e.g. hollow markers or red x),
crossover scale annotated if finite, `tiledlayout`, `'Interpreter','none'`,
`figure('Color','w', ...)`.

### Step 2.4 -- wire into summary output

Extend `pipeline_save_summary.m`'s per-channel loop to also write
`S.dfa_alpha1(k)`, `S.dfa_alpha2(k)`, `S.dfa_alphaFull(k)`, `S.dfa_nCross(k)`,
`S.dfa_nWindows{k}` alongside the existing `S.meanExcess`, `S.meanRate`,
following the same per-channel array/cell pattern already in that
function.

## 3. Shared -- `pipeline_params.m` additions

Add this block (values match `dfaGapAware.m`'s built-in defaults, so this
step only makes them visible/tunable in the one place the rest of the repo
expects tunables to live):

```matlab
% ---- Step 7: gap-aware DFA (fractal scaling) --------------------------
P.dfaMinScale           = 4;    % smallest box size n (beats/bins)
P.dfaMaxScaleFrac       = 0.25; % largest box size as frac of longest valid run
P.dfaNScales            = 20;   % number of log-spaced box sizes
P.dfaShortRange         = [4 16];  % alpha1 fit range
P.dfaLongRangeMinScale  = 16;   % alpha2 fit uses scales above this
P.dfaMinWindowsPerScale = 8;    % scale excluded from fit below this window count
P.dfaDetrendOrder       = 1;    % linear detrending
P.dfaMinRunLen          = 4;    % skip valid runs shorter than this
P.dfaBeatGapTolSec      = 1e-3; % cardiac-only; floating-point slop tolerance
                                 % for detecting a dropped beat between two
                                 % surviving RR_intervals entries -- NOT a
                                 % physiological threshold, keep small
```

## 4. Files to create / edit

| File | Action |
|---|---|
| `dfaGapAware.m` | Create (delivered alongside this plan -- copy in as-is) |
| `dfaRR_gapAware.m` | Create (§1) |
| `step7_dfa_report.m` | Create (§2.2-2.3) |
| `HR_BR_HRVAnalysis_new.m` | Edit: un-comment + replace DFA call (§1), extend `out`/`save` |
| `pipeline_save_summary.m` | Edit: add per-channel DFA fields (§2.4) |
| `pipeline_params.m` | Edit: add Step 7 block (§3) |
| `defineMetricsUI.m` (or wherever `defaultWindowedMetricSpecs()` lives) | Optional: register `dfa_alpha1`/`dfa_alpha2` as plottable metrics like `HRV`/`RMSSD` |
| `test_dfaGapAware.m` | Create (§5) -- includes `plot_dfaTest_diagnostic` helper (§5.0) |
| `test_figs/` | Created automatically by the test file; every test's saved `.png`/`.fig` diagnostic lands here for manual review |

## 5. Test suite

All tests below should be written as an actual runnable MATLAB test file
(e.g. `test_dfaGapAware.m`, using plain `assert` or the `matlab.unittest`
framework if preferred) and run before either wrapper touches real data.

**Every test must save a diagnostic figure, not just print assert
pass/fail.** Claude Code cannot always execute MATLAB in this environment,
so a test suite that only reports results as console text may end up never
actually having been run by anything -- the figures are the fallback
verification path: if Claude Code cannot execute the tests, run
`test_dfaGapAware.m` yourself in MATLAB and review the saved `.png` files
directly. Even when Claude Code can execute them, the figures let you
sanity-check that a passing assertion actually corresponds to a
sensible-looking fit rather than, say, a degenerate all-NaN case that
happens to satisfy a loosely-written assertion.

### 5.0 Shared diagnostic plot helper

Build this once, call it from every test below (`savePath` argument makes
each test's figure land in a dedicated `test_figs/` folder with a
per-test filename, e.g. `test_figs/T1_H0.9.png`).

```matlab
function fig = plot_dfaTest_diagnostic(x, validMask, out, trueAlpha, titleStr, savePath)
% PLOT_DFATEST_DIAGNOSTIC Standardized 3-panel diagnostic figure for
% visually verifying dfaGapAware output. Used by the test suite AND by
% manual review when Claude Code cannot execute MATLAB directly.
%
% x         - input series passed to dfaGapAware
% validMask - logical mask used (for shading gaps on panel 1); pass
%             true(size(x)) if the call used P.dfaRunsOverride instead
% out       - dfaGapAware() output struct
% trueAlpha - (optional) known ground-truth alpha/H for synthetic tests;
%             NaN or omit for real-data / mock-struct calls with no
%             ground truth
% titleStr  - figure title / test ID, e.g. 'T4: long/rare gaps, H=0.9'
% savePath  - full path with NO extension; saves both .fig and .png

if nargin < 4 || isempty(trueAlpha); trueAlpha = NaN; end
if nargin < 5; titleStr = ''; end

fig = figure('Color', 'w', 'Name', titleStr, 'Position', [100 100 900 800]);
tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, titleStr, 'Interpreter', 'none');

% Panel 1: input series with gaps shaded (same shading idiom as
% step2_noise_sigma.m's shade_runs)
ax1 = nexttile(tl);
plot(ax1, x, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold(ax1, 'on');
shade_invalid(ax1, ~validMask);
if ~isempty(out.runs)
    yl = ylim(ax1);
    for r = 1:size(out.runs, 1)
        plot(ax1, [out.runs(r,1) out.runs(r,1)], yl, 'g:', 'LineWidth', 0.8);
        plot(ax1, [out.runs(r,2) out.runs(r,2)], yl, 'g:', 'LineWidth', 0.8);
    end
end
title(ax1, sprintf('Input series (%d pts, %.0f%% valid, %d run(s))', ...
    numel(x), 100*mean(validMask), size(out.runs,1)), 'Interpreter', 'none');
ylabel(ax1, 'x'); xlabel(ax1, 'index');

% Panel 2: log-log F(n) with fitted ranges, excluded scales, crossover
ax2 = nexttile(tl);
loglog(ax2, out.nVals, out.F, 'ko-', 'MarkerFaceColor', 'k'); hold(ax2, 'on');
if isfinite(out.alpha1)
    idxS = out.nVals >= out.P.dfaShortRange(1) & out.nVals <= out.P.dfaShortRange(2);
    plot(ax2, out.nVals(idxS), out.F(idxS), 'r-', 'LineWidth', 2);
end
if isfinite(out.alpha2)
    idxL = out.nVals > out.P.dfaLongRangeMinScale;
    plot(ax2, out.nVals(idxL), out.F(idxL), 'b-', 'LineWidth', 2);
end
if isfinite(out.nCross)
    xline(ax2, out.nCross, '--', sprintf('crossover n=%.1f', out.nCross));
end
if ~isempty(out.excludedScales)
    for es = out.excludedScales(:)'
        xline(ax2, es, 'r:', 'HandleVisibility', 'off');
    end
end
grid(ax2, 'on'); xlabel(ax2, 'n (log scale)'); ylabel(ax2, 'F(n) (log scale)');
legTxt = sprintf('alpha1=%.3f (R2=%.2f)  alpha2=%.3f (R2=%.2f)  alphaFull=%.3f (R2=%.2f)', ...
    out.alpha1, out.R2_1, out.alpha2, out.R2_2, out.alphaFull, out.R2_full);
if isfinite(trueAlpha)
    legTxt = sprintf('%s  |  TRUE=%.3f', legTxt, trueAlpha);
end
title(ax2, legTxt, 'Interpreter', 'none', 'FontSize', 9);

% Panel 3: pooled window count per scale (QC floor) + local scale-dependent slope
ax3 = nexttile(tl);
yyaxis(ax3, 'left');
bar(ax3, out.nVals, out.nWindows, 'FaceColor', [0.7 0.7 0.9]); hold(ax3, 'on');
yline(ax3, out.P.dfaMinWindowsPerScale, 'r--');
ylabel(ax3, 'pooled windows (QC)');
yyaxis(ax3, 'right');
plot(ax3, out.localAlpha.n, out.localAlpha.slope, 'g.-');
ylabel(ax3, 'local slope alpha(n)');
xlabel(ax3, 'n');
title(ax3, 'Window-count floor (bars, red line = P.dfaMinWindowsPerScale) and local exponent (line)', ...
    'Interpreter', 'none', 'FontSize', 8);

if nargin >= 6 && ~isempty(savePath)
    [d,~,~] = fileparts(savePath);
    if ~isempty(d) && ~isfolder(d); mkdir(d); end
    savefig(fig, [savePath '.fig']);
    exportgraphics(fig, [savePath '.png'], 'Resolution', 150);
end
end

function shade_invalid(ax, mask)
if ~any(mask); return; end
yl = ylim(ax);
d = diff([0; mask(:); 0]);
starts = find(d == 1); ends = find(d == -1) - 1;
for i = 1:numel(starts)
    patch(ax, [starts(i) ends(i) ends(i) starts(i)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.85 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
ylim(ax, yl);
end
```

### 5.1 Core engine correctness (no gaps)

- **T1 -- Known-H recovery, ungapped.** Generate fractional Gaussian noise
  with a known Hurst exponent H (e.g. via `wfbm(H, N)` if Wavelet Toolbox
  is available, else a documented cumulative-sum construction), N >= 5000.
  Run `dfaGapAware(x, true(N,1), P)`. Assert `abs(out.alphaFull - H) <
  0.05` (loosen if the specific generator used has known bias -- state the
  tolerance basis in a comment). Repeat for H = 0.3 (anti-persistent),
  H = 0.5 (white noise), H = 0.9 (strongly persistent) -- three separate
  assertions, not just one H value, since Ma et al. (2010) found gap
  sensitivity differs sharply between persistent and anti-persistent
  signals and later tests depend on this baseline being solid for both
  regimes.
  **Figure:** call `plot_dfaTest_diagnostic` for each of the three H
  values with `trueAlpha` set to that H, save to
  `test_figs/T1_H0.3.png`, `test_figs/T1_H0.5.png`, `test_figs/T1_H0.9.png`.
  Pass looks like: panel 2's `alphaFull` annotation lands close to the
  `TRUE=` value; panel 3's bars are all comfortably above the red QC
  floor line at every scale (no gaps means no scale should be excluded).
- **T2 -- Agreement with `dfaRR.m` on ungapped data.** Same synthetic
  series as T1 (H=0.7 case), run both `dfaRR.m` and
  `dfaGapAware(x, true(N,1), P)`. Assert `abs(out.alpha1 - alpha1_old) <
  0.02` and same for `alpha2` (tight tolerance -- these should be
  near-identical since gap-handling is a no-op when there are no gaps).
  **Figure:** save `test_figs/T2_H0.7.png` via `plot_dfaTest_diagnostic`
  with `trueAlpha = 0.7`, and separately print/annotate `dfaRR.m`'s
  `alpha1`/`alpha2` as a text overlay or console-printed comparison line
  next to the figure filename, so the two numbers can be eyeballed
  side by side even without re-running the assertion.

### 5.2 Core engine correctness (synthetic gaps)

- **T3 -- Short/sparse gaps, persistent signal.** Take the H=0.9 series
  from T1, punch out 5-10 gaps of 1-5 seconds each (short, matching the
  low end of your actual motion-artifact range), scattered randomly.
  Build `validMask` accordingly. Assert `abs(out.alphaFull - 0.9) < 0.08`
  (looser than T1 since gaps are present, but should still be close --
  per Ma et al., persistent signals are robust even to large data loss).
  **Figure:** `test_figs/T3_H0.9_shortGaps.png`, `trueAlpha = 0.9`. Pass
  looks like: panel 1 shows the gap shading clearly separating several
  green-dotted run boundaries; panel 2's fit still tracks close to
  `TRUE=0.9` despite the gaps.
- **T4 -- Long/rare gaps, persistent signal.** Same H=0.9 series, punch out
  1-3 gaps of 60-180 seconds within a ~600 s series (matching your
  "minutes-long motion artifact in a 10-min baseline" case). Assert the
  function still returns finite `alpha1`/`alpha2` (not NaN) as long as at
  least one clean run exceeds `P.dfaMinRunLen * P.dfaMinScale`-ish length,
  and that `out.excludedScales` is non-empty at the large-n end (this is
  the EXPECTED, correct behavior -- the test should assert truncation
  happened, not treat it as a failure).
  **Figure:** `test_figs/T4_H0.9_longGaps.png`, `trueAlpha = 0.9`. Pass
  looks like: panel 1 shows 1-3 wide shaded gap regions with only 2-4
  runs total; panel 2 shows red dotted excluded-scale lines clustered at
  the large-n end only; panel 3's bars drop below the red QC floor line
  exactly at those same large-n scales, visually confirming the exclusion
  reason.
- **T5 -- Long/rare gaps, anti-persistent signal.** Same as T4 but with
  H=0.3. Per Ma et al. (2010), anti-persistent signals are NOT robust to
  gaps -- assert that `out.alphaFull` measurably drifts toward 0.5 relative
  to the ungapped H=0.3 baseline from T1, and add a comment noting this
  drift is an expected literature-documented effect, not a bug, so nobody
  "fixes" it later by mistake.
  **Figure:** `test_figs/T5_H0.3_longGaps.png`, `trueAlpha = 0.3`. Pass
  looks like: panel 2's `alphaFull` visibly sits BETWEEN `TRUE=0.3` and
  0.5, closer to 0.5 than T1's ungapped H=0.3 figure was -- put the T1
  H=0.3 figure and this one side by side when reviewing, since the drift
  is only obvious in comparison, not from this figure alone.
- **T6 -- Gap severe enough to break the whole run.** Construct a series
  where valid runs are all shorter than `P.dfaMinRunLen`. Assert
  `dfaGapAware` throws the documented `dfaGapAware:noValidRuns` error
  rather than returning a silently-wrong result.
  **Figure:** none expected (the function should error before producing
  `out`) -- instead, wrap the call in try/catch and save a trivial text
  figure or console log confirming the error identifier matched
  `dfaGapAware:noValidRuns` exactly, so a human reviewer can confirm it
  failed for the RIGHT reason and not some unrelated bug.
- **T7 -- Too few windows at large scales only.** Construct a series where
  small scales have plenty of pooled windows but large scales fall below
  `P.dfaMinWindowsPerScale`. Assert `alpha1`/`R2_1` are finite while
  `alpha2`/`R2_2` are NaN, and that the excluded large scales appear in
  `out.excludedScales`.
  **Figure:** `test_figs/T7_largeScaleExclusion.png`. Pass looks like:
  panel 2 shows a red fitted line (alpha1) at small n and NO blue line
  (alpha2 is NaN) at large n; panel 3's bars visibly cross below the red
  QC floor line partway through the scale range, exactly where the blue
  fit stops appearing in panel 2 -- the two panels should tell the same
  story.

### 5.3 `P.dfaRunsOverride` path

- **T8 -- Runs override matches equivalent validMask.** Construct a
  validMask with a known gap pattern, run `dfaGapAware` normally. Then
  construct the equivalent runs list by hand and run again via
  `P.dfaRunsOverride` with `validMask = true(N,1)`. Assert both calls
  return identical `out.F`, `out.nWindows`, `out.alpha1`, `out.alpha2`
  (this is a pure refactor-equivalence test -- the two code paths must
  agree exactly, not approximately).
  **Figure:** save both as `test_figs/T8_viaValidMask.png` and
  `test_figs/T8_viaRunsOverride.png`. Pass looks like: the two figures
  are visually indistinguishable in panels 2 and 3 (panel 1 will differ
  cosmetically since the second call's `validMask` argument is
  `true(N,1)`, so no red shading appears there -- that's expected, not a
  discrepancy).

### 5.4 Cardiac wrapper (`dfaRR_gapAware.m`)

- **T9 -- No dropped beats.** Construct `RR_intervals`/`RR_times` with
  exact consecutive timing (`RR_times(i+1) = RR_times(i) +
  RR_intervals(i)` for all i). Assert the built `runs` list is a single
  run spanning the whole series.
  **Figure:** plot `RR_times` on the x-axis vs. cumulative beat index on
  the y-axis (a simple `plot(RR_times, 1:n)` is enough here, not the full
  3-panel helper), with the single detected run boundary marked at each
  end. Save `test_figs/T9_noDroppedBeats.png`. Pass looks like: a single
  smooth line, no vertical break markers in the middle.
- **T10 -- One dropped beat.** Same construction, but insert one artificial
  jump in `RR_times` (e.g. add 2 seconds at index 50, simulating a beat
  that was dropped) without changing `RR_intervals` values around it.
  Assert `runs` correctly splits into two runs at index 50, and that
  `dfaRR_gapAware`'s output differs from what you'd get by (incorrectly)
  treating the series as one unbroken run -- confirm this by also running
  `dfaGapAware(RR_intervals, true(N,1), P)` (ignoring the break) and
  asserting the two `alpha1` values are NOT the same, demonstrating the
  break-detection actually changes the computation rather than being a
  no-op.
  **Figure:** run `plot_dfaTest_diagnostic` on BOTH the correct
  (break-aware) and incorrect (break-ignored) calls, saving
  `test_figs/T10_correct_withBreak.png` and
  `test_figs/T10_incorrect_ignoredBreak.png`. Pass looks like: the
  break-aware figure's panel 1 shows a green dotted boundary exactly at
  index 50 with two separate runs; the ignored-break figure shows one
  continuous run (no boundary at 50) and a visibly different `alpha1` in
  its panel 2 title -- reviewing the two side by side should make the
  practical impact of break-detection obvious, not just a number in a
  console log.
- **T11 -- Floating-point slop does not falsely break a run.** Construct
  exact consecutive timing as in T9, then perturb every `RR_times` entry
  by up to 1e-6 seconds (simulating float rounding, not a real dropped
  beat). Assert `runs` is still a single run (i.e. `P.dfaBeatGapTolSec =
  1e-3` correctly absorbs this and does not fragment the series).
  **Figure:** same style as T9's figure, save
  `test_figs/T11_floatSlopNoFalseBreak.png`. Pass looks like: still a
  single smooth line with no break markers, despite the added jitter --
  visually near-identical to T9's figure, which is the point.

### 5.5 Nerve wrapper (`step7_dfa_report.m`)

- **T12 -- Mock `D` struct, no artifacts.** Build a minimal mock `D` with
  synthetic `D.fr_hz` (known H, as in T1) and `D.fr_validMask =
  true(size(D.fr_hz))` (field name pending §2.0 confirmation -- update
  this test once real field names are known). Run `step7_dfa_report(D, P,
  false)`. Assert `D.dfa.alphaFull` is close to the known H.
  **Figure:** even though `plotMode = false` in the assertion call above,
  separately call `plot_dfaTest_diagnostic` directly on `D.dfa` and save
  `test_figs/T12_mockD_noArtifacts.png` (do this as an explicit step in
  the test, not by relying on `step7_dfa_report`'s own `plotMode` branch,
  so the test's figure output doesn't depend on that branch being
  correctly wired yet).
- **T13 -- Mock `D` struct, with artifacts.** Same as T12 but with a
  `validMask` containing both short and long false-runs, matching the
  "ms to minutes" real artifact distribution. Assert the function runs
  without error and produces a `D.dfa` struct with the expected fields
  populated (finite where scales have enough windows, NaN where they
  don't, per T7's pattern).
  **Figure:** `test_figs/T13_mockD_withArtifacts.png`, same approach as
  T12. Pass looks like: panel 1 shows a mix of short and long shaded
  regions; panel 3 shows the QC floor being crossed only at large n, same
  pattern as T4/T7.
- **T14 -- Real-data smoke test.** Once §2.0 is resolved and real field
  names are wired in, run `step7_dfa_report` on 2-3 real baseline files
  (expected: low artifact burden, most scales retained) and 2-3 real
  recovery files (expected: heavier artifact burden, more scales
  excluded, possibly NaN `alpha2`). This is a smoke test, not an
  assertion-based test -- the goal is confirming it runs end-to-end on
  real `.mat` files and manually eyeballing that `nWindows`/
  `excludedScales` look sane (e.g. compare a heavily-blanked recovery file
  against `browseMotionArtifacts.m`'s view of the same file's annotated
  artifacts, to catch any unit mismatch between the artifact annotation
  format and what feeds into `validMask`).
  **Figure:** this is the one test where `step7_dfa_report`'s own
  `plotMode = true` diagnostic plot (§2.3) IS the figure -- save one per
  file as `test_figs/T14_<condition>_<animal>_<phase>.png`, following the
  existing `saveFigureAllFormats`/`saveFigure` naming pattern used
  elsewhere in the repo. Since there's no ground truth here, "pass" is
  necessarily a human judgment call: does the baseline file's figure look
  like T3 (light gaps, most scales retained) and does the recovery file's
  figure look like T4 (heavier gaps, large-n scales excluded)? If a
  baseline file's figure instead looks like T4, that's a red flag worth
  investigating before trusting its `alpha` values.

### 5.6 Regression / integration

- **T15 -- `HR_BR_HRVAnalysis_new.m` end-to-end.** After un-commenting the
  DFA call, run the full function on one real condition file. Assert it
  completes without error, `dfa_alpha1`/`dfa_alpha2` land in a
  physiologically plausible range for rat HRV (roughly 0.5-1.3; flag, do
  not silently accept, anything far outside this as a likely masking bug
  rather than a real finding), and the new fields appear correctly in the
  saved `_HRVMeasures.mat`.
  **Figure:** call `plot_dfaTest_diagnostic` on `dfaOut` from this real
  run, save `test_figs/T15_<conditionLabel>.png`. Since
  `HR_BR_HRVAnalysis_new.m` already produces its own Figure 1/Figure 2
  diagnostics (with `invalidMask`/`edgeMask` shading) when `figtrue=true`,
  open those alongside this new DFA figure and confirm the shaded invalid
  regions in this figure's panel 1 line up with the same regions shown in
  the existing Figure 1 panel 1 -- this cross-check catches any mismatch
  between the RR-level `runs` list this plan introduces and the
  sample-level `invalidMask` the rest of the file already uses.
