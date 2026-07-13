# Session Handoff — GEMS Physiology Plotting Toolkit

Complete state as of commit `13f33e1`. Meant to be read cold by a fresh Claude
session with no prior context.

---

## 0. Project identity

- **User**: Arjun Bindu Sunil (git author) / Andrea Elizabeth Biju (macOS account
  and GitHub user `AndreaEBiju`). Same physical person.
- **Repo**: <https://github.com/AndreaEBiju/processing_new> (private).
- **Local project root**: `/Volumes/PortableSSD/processing_new` on macOS.
  Files may also be accessed from `C:\Users\nouch\Desktop\processing_new` on
  Windows — the queue cache stores absolute paths so cross-machine sync via
  `git pull` won't automatically re-point to a different mount.
- **Data location**: Google Drive shared folder
  `G:\Shared drives\BIONICs Lab Workspace\Project Folders\GEMS\Survivals\...`
  on Windows, or a similar `/Volumes/...` path on macOS. Files are on Google
  Drive's virtual filesystem, so `load()` triggers download on first access;
  disk caches (see §5) are essential.

## 1. Domain / experimental design

BIONICs lab study: rat multi-modal gastric+vagus stimulation.

- **Recording**: one .mat file per (animal, condition, phase). Columns:
  `RVN, LVN, ANT1, ANT2, ANT3` (2 nerve, 3 antral stomach), fs = 25 kHz.
- **Phases**: `baseline` (no stim) and `recovery` (post-stim). Stim itself is
  contaminated by artefact and not analysed.
- **Conditions**: modality codes
  - `M10, M50, M100` — mechanical only
  - `E10, E100, E1000` — electrical only
  - `M<m>E<e>` (or `E<e>M<m>`) — combined
  Typical set has 15 conditions. 3 animals per condition, 2 phases per animal
  ⇒ ≈ 90 files.
- **Animals**: 4 IDs in this dataset: `f` (FRE), `j` (JEL/Jel/jel), `l`
  (Lol/LOL/lol), `o` (ORE/Ore). Single-letter animal IDs come from
  `lower(tokens{2}(1))`.

## 2. Files in the repo (all live at project root)

### Analysis / batch pipeline (pre-existing, not authored this session)

- `batch_process.m` — parallel HR/HRV/slow-wave analysis driver with file queue UI.
- `HR_BR_HRVAnalysis_new.m` — saves `<stem>_HRBR.mat` and `<stem>_HRVMeasures.mat`.
- `slowWaveAnalysis_new.m` — saves `<stem>_slowWaves.mat`.
- Other helpers: `blank_stim_spikes_nan.m`, `detectSortNerveSpikesECAP.m`,
  `mergeSegments.m`, `splitStimRecovery.m`, etc.

### Building blocks (authored/refactored this session)

- **`pubfig_setup.m`** — publication-quality plot defaults. Takes `'Theme'`
  parameter (`'light'` / `'dark'`) that controls figure/axes background,
  foreground, grid, legend, and categorical palette. Stashes to appdata:
  `PUB_PALETTE`, `PUB_FGCOL`, `PUB_BGCOL`, `PUB_THEME`. Configurable:
  `BaseFontSize`, `LineWidth`, `MarkerSize`, `ColorOrder`, `PlotsDir`,
  `Renderer`, `EnableLaTeX`.
  - Default fonts / weights configurable per call.
  - Palettes:
    - Light: mid-luminance vivid — blue / vermillion / teal-green / purple /
      deep gold / sky blue / maroon / dark grey.
    - Dark: brighter high-luminance versions of the same hues.
  - Sets `defaultFigureInvertHardcopy='off'`, which triggers
    `Ignoring 'InvertHardCopy' property` info warnings from `print()` —
    suppressed in save helpers.

- **`boxScatterPlot.m`** — generic box + jittered scatter, optionally grouped.
  - Input: cell array `G × S` (rows = primary groups, cols = subgroups) OR
    vector for ungrouped.
  - Options: `GroupLabels`, `SubgroupLabels`, `YLabel`, `XLabel`, `Title`,
    `Colors`, `BoxAlpha` (0.08), `BoxWidth`, `GroupGap`, `JitterWidth`,
    `MarkerSize`, `MarkerAlpha`, `ShowOutliers`, `ConnectPaired`,
    `PairLineStyle` (`':'` default), `ColorBySubject`, `SubjectColormap`
    (`'auto'` = per-box hue lightness ramp), `CycleMarkers`, `Notch`,
    `Parent`, `Legend`.
  - Key behaviour: with `ColorBySubject='auto'` uses HSV ramp of each box's
    hue (dark-saturated → pale-tinted, computed via `rgb2hsv`/`hsv2rgb`) so
    each subject gets a distinguishable shade in the box's hue family.
    Per-subject jitter is shared between paired-line endpoints and dot
    centres so lines land exactly on dots.
  - All labels use `'Interpreter','none'` (LaTeX is set as default by
    `pubfig_setup`, chokes on `_`, `%`, `\`).

- **`boxViolinPlot.m`** — KDE violin + inner box, same input shape as
  `boxScatterPlot`. Uses `ksdensity` (Statistics Toolbox). Options include
  `ViolinAlpha` (0.30), `BoxInnerFrac` (0.30 — inner box width as fraction
  of slot), `KdeBandwidth`, `KdePoints`, `ShowScatter`. NaN-safe (cells
  with <2 unique values skip the violin).

- **`convertUnits.m`** — scalar/array unit conversion. Handles:
  time (s / ms / min / hour), rate (Hz / bpm / cpm / rpm),
  fraction (% ↔ fraction), voltage (V / mV / uV). Empty/matching units
  return value unchanged; unknown pair warns once and returns unchanged.

### UI wrappers (reusable across plot types)

- **`buildFileQueue.m`** — interactive uifigure UI to build/edit a file queue.
  Columns: `File | Animal | Condition | Phase`, all editable inline.
  - `Add files...` (uigetfile multi-select)
  - `Add folder...` (uigetdir + inputdlg for glob pattern; default
    `*_v*_blankmotion.mat`; recursive `dir(fullfile(root,'**',pattern))` with
    uiprogressdlg during scan)
  - Version-filter rule: `isVersionedName` regex
    `_v\d+(\.\d+)+_` — only files with a semver-style version tag
    (`_v0.x.x_...`) pass. Silently drops unversioned matches.
  - Auto-parses:
    - `animal` = first letter of 2nd `_`-token, lowercased.
    - `condition` = first token matching `^[ME]\d+([ME]\d+)?$`.
    - `phase` = `'recovery'` if basename contains `recovery` OR `stim_rec`;
      `'baseline'` if it contains `baseline` OR `_bl_`; default `'baseline'`.
  - Persists to cache file passed as arg; preserves unknown fields (e.g.
    `groups`, `metrics`, `windowsStr`) round-trip.
  - Diagnostic `[buildFileQueue] Continue clicked. N row(s) in table.` on
    Continue.

- **`defineGroupsUI.m`** — uifigure table for defining user-groups of
  conditions. Table columns: `Label | Conditions (comma-separated)`.
  Buttons: Add group, Remove selected, Continue.
  - Extracted from `plotAveragedMetrics` earlier so `plotWindowedMetrics`,
    `plotWindowedViolins`, `plotTimeTraces`, `plotRecoveryDynamics` all use
    the same one. `plotSynergyHeatmaps` does NOT use it (fixed 3×3 layout).

- **`defineMetricsUI.m`** — 9-column uitable of metric specs:
  `Label | Stored units | Plot units | File suffix | Field | Series field |
   Time field | Aggregator | Channel`. Includes:
  - Auto-migration from legacy 5-col, 6-col, 7-col cell tables + old
    `units` field. Auto-extracts units from `'HRV (ms)'`-style labels.
  - `Add row`, `Remove selected`, `Reset to defaults`, `Inspect .mat
    file...` (opens picked .mat in a scrollable uifigure+uitable showing
    every variable's name/class/size/bytes; click a row to copy the name
    to clipboard).
  - Aggregator dropdown: `auto | mean | median | max | min | first | last
    | sum | scalar`. `auto` returns value as-is if already scalar, else
    omitnan mean.

- **`loadMetric.m`** — pulls ONE aggregated scalar per (source, spec). Path
  resolution tries in order:
  1. `<base><suffix>` (baseline / post-split recovery keep `_blankmotion`
     in the output name)
  2. `<stem><suffix>` (with `_blankmotion` stripped)
  3. `<stem>_recovery<suffix>` (only if `stem` contains `stim_rec`) —
     `batch_process` appends `_recovery` to the condition when splitting
     stim_rec files.
  Optional channel selection then `switch aggregator`. Returns NaN on any
  failure.

- **`loadMetricSeries.m`** — same path resolution as `loadMetric` but returns
  `(y, t)` for series+time fields. Includes case-insensitive fieldname
  resolution and diagnostic `[loadMetricSeries] MISSING ...` printing when
  fields aren't found; falls back to `whos '-file'` for casing/punct
  mismatches.

- **`loadAllMetricsCached.m`** — parallel (`parfeval`) load of scalars across
  files with disk cache. Used only by `plotAveragedMetrics`. Batches loads
  by file, treats cached NaN as miss (auto-heal), skips writing NaN.

- **`loadAllSeriesCached.m`** — **THE canonical shared cache**. Parallel
  (`parfeval`) load of `(y, t)` per (file, metric). Batches by result-file
  path (one `load()` per unique file with all needed vars — massive speedup
  on v7 .mat files where partial loading still requires end-to-end scan).
  Cache: `gemsplots_series_cache.mat`. Every windowed/violin/time-trace/
  synergy/recovery-dynamics plotter uses this.

- **`loadAllWindowedCached.m`** — thin wrapper around `loadAllSeriesCached`.
  Derives `baselineByFile` (full-series mean) and `windowByFile` (mean over
  `t - t0 <= W`) in-memory from cached series. No separate cache.

### Plot wrappers (entry points)

All follow the same UI sequence: queue → groups → metrics → save folder →
auto-close timer. All share `gemsplots_queue.mat` for state and
`gemsplots_series_cache.mat` for data (except `plotAveragedMetrics` which
uses `gemsplots_metrics_cache.mat` for scalars).

1. **`plotAveragedMetrics.m`** — box+scatter of averaged scalars.
   One figure per (metric, user-group). X = conditions, subgroups =
   baseline/recovery. `ConnectPaired=true` (dotted lines link matched
   subjects across baseline↔recovery within a condition). `ColorBySubject=true`.
   Uses `loadMetric` for scalar averaged fields (e.g. `avgHeartRate`).

2. **`plotWindowedMetrics.m`** — box+scatter of windowed normalised values.
   One figure per (metric, user-group). X = conditions, subgroups = windows.
   Values: `(rec_window_mean − baseline_full_mean) / baseline_full_mean`,
   per animal, per (condition, window). Windows configured via
   `defineWindowsUI` (comma-separated seconds + `'full'`). Cached in
   `state.windowsStr`.

3. **`plotWindowedViolins.m`** — same data + same axis layout as
   `plotWindowedMetrics` but box+violin instead of box+scatter. **Was
   originally windows-on-x, conditions-as-subgroups; flipped to match
   plotWindowedMetrics per user request.**

4. **`plotTimeTraces.m`** — per (metric, condition) recovery time trace
   with mean±SEM band + individual animal lines. Iterates every condition
   in the queue (not user-groups). Each animal's time vector normalised so
   `t0 = 0`. Common grid = `0:median(dt):min(tEnds)`, `interp1` all to that.

5. **`plotSynergyHeatmaps.m`** — 3×3 M×E synergy heatmap per metric.
   `synergy(M,E) = norm̂(M,E) − [norm̂(M-alone) + norm̂(E-alone)]` on
   baseline-normalised values (not the raw Δ the Dashboard_Guide uses).
   Parses condition strings via regex `^M(\d+)E(\d+)$`, `^E(\d+)M(\d+)$`,
   `^M(\d+)$`, `^E(\d+)$`. Diverging red-white-blue colormap (hand-rolled
   `rdbuColormap(n)`) symmetric at 0. Cells annotated with value + n. NaN
   cells greyed. **Does not use user-groups.**

6. **`plotRecoveryDynamics.m`** — 2×2 tiled per (metric, user-group) of the
   four Dashboard_Guide section-5 features on the **baseline-normalised**
   recovery trace `yN = (y_rec - μ_base) / μ_base`:
   - `peakExc = max(yN) - min(yN)` (fractional)
   - `slope30` = OLS slope of `yN(t)` on `t ∈ [0, 30 s]` (1/s)
   - `final30Mean` = mean(yN) over last 30 s (fractional)
   - `AUCdev = trapz(t, |yN - final30Mean|)` (seconds)
   Each panel is `boxScatterPlot` with `ColorBySubject=true`.

### Diagnostic helpers

- **`diagnoseMetricCoverage.m`** — reads `gemsplots_queue.mat`, prints:
  - Per-metric coverage (found / total + % bar).
  - Per-source X/. matrix.
  - First 5 missing output files it looked for.
  - Sources with zero coverage.
  Uses only `dir()` (no `load`), so it's fast even over Google Drive.

- **`repairQueuePhases.m`** — re-derives the `phase` column of the cached
  queue using the current `parseFilename` rule. Prints before/after counts.
  Used once when phase parsing was updated to treat `stim_rec` as
  `'recovery'`.

### Demo

- **`demo_boxScatterPlot.m`** — synthetic-data demo showing ungrouped +
  grouped-paired + grouped-unpaired usage.

## 3. Cache files (NOT tracked in git — `.gitignore` broadens to `*.mat`)

Living at project root; each stored to disk between runs.

| Cache file | Contents | Used by |
|---|---|---|
| `gemsplots_queue.mat` | `files, animal, condition, phase, groups, metrics, windowsStr` (each cell/struct); this is the ONE queue-state cache | every plot wrapper |
| `gemsplots_metrics_cache.mat` | `containers.Map` keyed `filepath\|\|specKey\|\|baseline_or_win_key` → `struct('mtime', mtime, 'value', scalar)` for `plotAveragedMetrics` | `plotAveragedMetrics` |
| `gemsplots_series_cache.mat` | `containers.Map` keyed `filepath\|\|specKey` → `struct('mtime', mtime, 'y', vec, 't', vec)` — the shared full-series cache | `plotWindowedMetrics`, `plotWindowedViolins`, `plotTimeTraces`, `plotSynergyHeatmaps`, `plotRecoveryDynamics` |
| ~~`gemsplots_windowed_cache.mat`~~ | ORPHANED. Was separate windowed-scalar cache; replaced by in-memory derivation from the shared series cache. Safe to `delete()`. | none |

Cache invalidation is mtime-based per source file. Cached NaN treated as
miss (so a failed load auto-retries next run). Never writes NaN to cache.

## 4. Data + naming conventions

Source-file basename layout (typical):
```
<condition>_<animal>_<misc>_<time>_notched_v0.x.x_blankmotion.mat
```

Examples from the user's actual data:
- `M50E100_lol_CME1_bl_2124_notched_v0.2.2_blankmotion.mat` (baseline, animal l)
- `M50E100_lol_CME1_stim_rec_2134_notched_v0.2.2_blankmotion.mat` (stim_rec source → recovery)
- `E1000_FRE_E1000_bl_1356_notched_v0.2.2_blankmotion.mat` (baseline, animal f)

`batch_process` output naming:
- Non-`stim_rec` sources: output = `<base>_HRBR.mat` (keeps `_blankmotion` in
  the name). E.g. `..._v0.2.2_blankmotion_HRBR.mat`.
- `stim_rec` sources that got split: output = `<stem>_recovery_HRBR.mat`
  where `<stem>` = base minus `_blankmotion`.
  E.g. `M50E100_lol_CME1_stim_rec_2134_notched_v0.2.2_recovery_HRBR.mat`.

`loadMetric` / `loadMetricSeries` try both patterns.

Result-file variable names (per the Dashboard_Guide):
- `_HRBR.mat`: `avgHeartRate`, `avgBreathRate`, `heartRateSeries`,
  `breathRateSeries`, `metrics_t`, `t`, ...
- `_HRVMeasures.mat`: `hrv`, `rmssd`, `pnn5`, `sd1`, `sd2`, `sampEn`,
  `appxEn`, `hrv_series`, `rmssd_series`, `pnn5_series`, `sd1_series`,
  `sd2_series`, `sampEn_series`, `metrics_t`.
- `_slowWaves.mat`: `slowWaveRateSeries` (columns = ANT1..3),
  `slowWaveRateTime`, `slowWaveTimeSeries`, `avgSlowWave`.

Standard metric spec defaults (bundled in every wrapper as
`defaultWindowedMetricSpecs()`):

| Label | UnitsIn | UnitsOut | Suffix | Field | SeriesField | TimeField | Aggregator | Channel |
|---|---|---|---|---|---|---|---|---|
| HR | bpm | bpm | _HRBR.mat | avgHeartRate | heartRateSeries | metrics_t | auto | |
| Breathing rate | bpm | bpm | _HRBR.mat | avgBreathRate | breathRateSeries | metrics_t | auto | |
| HRV | s | **ms** | _HRVMeasures.mat | hrv | hrv_series | metrics_t | auto | |
| pNN5 | % | % | _HRVMeasures.mat | pnn5 | pnn5_series | metrics_t | auto | |
| RMSSD | s | **ms** | _HRVMeasures.mat | rmssd | rmssd_series | metrics_t | auto | |
| Sample entropy | | | _HRVMeasures.mat | sampEn | sampEn_series | metrics_t | auto | |
| SD1 | s | **ms** | _HRVMeasures.mat | sd1 | sd1_series | metrics_t | auto | |
| SD2 | s | **ms** | _HRVMeasures.mat | sd2 | sd2_series | metrics_t | auto | |
| SW rate ANT1 | cpm | cpm | _slowWaves.mat | slowWaveRateSeries | slowWaveRateSeries | slowWaveRateTime | mean | 1 |
| SW rate ANT2 | cpm | cpm | _slowWaves.mat | slowWaveRateSeries | slowWaveRateSeries | slowWaveRateTime | mean | 2 |
| SW rate ANT3 | cpm | cpm | _slowWaves.mat | slowWaveRateSeries | slowWaveRateSeries | slowWaveRateTime | mean | 3 |

Note: HRV/RMSSD/SD1/SD2 stored in seconds; user wanted them plotted in ms —
`convertUnits('s', 'ms')` handles it at plot time (`plotAveragedMetrics`
applies it in `renderMetricGroupFigure`, `plotWindowedMetrics` derives
normalised values so units cancel and the conversion isn't applied there).

## 5. Normalization convention

**All windowed / dynamics / synergy math uses baseline-normalised effect**:
```
y_norm(t)   = (y_recovery(t) - μ_baseline) / μ_baseline
norm_scalar = (recovery_mean  - μ_baseline) / μ_baseline
```
where `μ_baseline = mean(y_baseline, 'omitnan')` for the matching
(animal, condition) baseline file.

This diverges from the Dashboard_Guide, which uses raw Δ = `recovery -
baseline`. User explicitly requested the fractional form and it's now
consistent across `plotWindowedMetrics`, `plotWindowedViolins`,
`plotSynergyHeatmaps`, `plotRecoveryDynamics`.

Synergy math:
```
synergy(M, E) = norm̂(M, E) - [norm̂(M-alone) + norm̂(E-alone)]
```
The `+ control` correction from the doc drops out because control_norm
= (base - base) / base = 0.

## 6. UI patterns / conventions

- Every wrapper prints `[<function>] ...` diagnostic lines to the
  Command Window (Continue clicks, drop reports, cache hit/miss counts).
- All uifigures use tabular `uitable` + Add / Remove / Continue buttons.
  MacOS R2025a sometimes fails to render `inputdlg` after a `uifigure`
  closes — every group-definition and text-input UI has therefore been
  ported to uifigure-based tables (`defineGroupsUI`) or kept as
  `inputdlg` only where testing showed reliability (window definition,
  auto-close prompt).
- All `printf` uses non-LaTeX text; every plot label uses
  `'Interpreter','none'` because `pubfig_setup` sets LaTeX as the
  interpreter default and it chokes on `_`, `%`, `\`.
- Figures explicitly forced to `WindowState='normal'` so they don't
  inherit `pubfig_setup`'s `maximized` default.
- Save helper `saveFigureAllFormats(fig, outDir)`:
  1. `.fig` first (preserves current theme).
  2. Temporarily forces `fig.Color='w'`, `axes.Color='w'`,
     `InvertHardcopy='off'`.
  3. `.png` via `exportgraphics` with `BackgroundColor='white'`, 200 dpi.
  4. `.svg` via `exportgraphics(..., 'ContentType','vector',
     'BackgroundColor','white')` first; falls back to `print -dsvg -vector`
     wrapped in `warning('off','all')` to suppress the `Ignoring
     InvertHardCopy` notice (R2024a+).
  5. Restores original figure/axes colours via `onCleanup`.

## 7. Known workflow gotchas

- **Google Drive latency**: `load()` on a virtual filesystem file blocks
  while the file downloads. First cold-cache run of `loadAllSeriesCached`
  can take 5–15 min on 90 v7 .mat files; warm cache is instant. macOS TCC
  ("Files and Folders" → "Removable Volumes") permission needed for
  Terminal / MATLAB / Claude to read external SSD.
- **v7 vs v7.3 .mat**: user's result files are v7 (default), so partial
  `load()` still requires end-to-end file scan. Batching by file (one
  `load()` per unique result file with all needed vars) is what makes
  the first run tolerable. Migrating files to v7.3 would give true
  random access. Not attempted here.
- **`clear functions` matters**: MATLAB caches function definitions
  across runs; edits to the plotter files are ignored until you
  `clear functions` (or `clear all`). Documented in every reply that
  bumps code.
- **Group typo silently produces blank column**: e.g. group definition
  `M100E50` when the queue has `M50E100`. Every wrapper now validates
  after `defineGroupsUI` returns and prints
  `!! UNKNOWN condition(s) in groups (not in queue, will plot blank)`.
- **Case-insensitive matching**: all lookups use `strcmpi` because
  animal IDs vary in case (`ORE`/`Ore`, `lol`/`Lol`/`LOL`); condition
  strings should be uppercase (parser upper-cases them).
- **Multiple trials per animal per condition**: some conditions have
  more than one trial from the same animal (e.g. `M50E100_lol_CME1_...`
  AND `M50E100_lol_CME2_...`). Current renderers per-animal use
  `find(..., 1)` (first-match) — trials 2+ are silently dropped from
  windowed/violin/dynamics plots. `plotTimeTraces` correctly plots all
  trials because it iterates recovery files directly. Not yet reconciled
  across all plotters; the user has noticed for `plotTimeTraces` but not
  yet asked for the same in the windowed plots.

## 8. Known data-quality state

From `diagnoseMetricCoverage` on the user's queue (~90 files):

- 89 / 90 sources have HRBR + HRV + slowWaves outputs.
- 1 fully-unanalysed source: `M50E10_JEL_M50E10_stim_rec_2134_...` —
  needs `batch_process` re-run.
- ~29 sources are missing only `_slowWaves.mat`: 1 for animal f,
  12 for animal l, 16 for animal o. Suggests `slowWaveAnalysis_new`
  wasn't run (or failed silently) for those sources.
- 1 source has ANT1 as NaN while ANT2/ANT3 valid:
  `M10E1000_Ore_CME_stim_rec_1253_...`

Reported to user in a per-tier summary during earlier troubleshooting.

## 9. Repo / git state

- Branch: `main`.
- Remote: `https://github.com/AndreaEBiju/processing_new.git` (HTTPS via
  gh token, not SSH).
- `.gitignore`:
  ```
  ._*
  .DS_Store
  ...
  *.mat
  gemsplots_queue.mat
  gemsplots_metrics_cache.mat
  plots/
  plots_demo/
  ```
  Caches explicitly listed for safety even though `*.mat` catches them.
- Earlier commit `a29358f` briefly tracked `gemsplots_queue.mat`; reverted
  in `5db7ac5`. Historical file still in that commit — accessible only to
  repo owner (private repo). User declined force-push scrub.
- Recent commits (newest first):
  ```
  13f33e1 plotRecoveryDynamics: features on baseline-normalised trace
  f56d4b9 Add plotRecoveryDynamics
  f2acaaa Add plotSynergyHeatmaps
  a5fcf39 plotWindowedViolins layout flip to match plotWindowedMetrics
  264f06d plotWindowedViolins drop diagnostic
  23e3ff2 plotTimeTraces per-figure dropped-trial log
  a11ee75 Unify cache (single series cache shared by all plotters)
  dcb2b3b Add plotTimeTraces + loadAllSeriesCached
  ```
  Full log via `git log --oneline`.

## 10. Open items / hints for next session

- User might request per-trial (not per-animal) granularity in windowed/
  violin plots — currently `find(..., 1)` picks first trial only.
  `plotTimeTraces` already does per-trial. Reconcile if asked.
- `plotAveragedMetrics.renderMetricGroupFigure` still has a "flag NaN
  queued sources" log — different style from the granular `reasons` /
  `condUsable` block in `plotWindowedMetrics.renderWindowedFigure`.
  Could be unified.
- `boxScatterPlot` has legacy `CycleMarkers` code path (default false)
  that isn't used anywhere any more; dead-ish.
- User hinted at possibly wanting nerve-spike features next
  (Dashboard_Guide section 4). Would follow same pattern:
  `loadAllSeriesCached` for whichever fields, `plotXxx.m` with 2×2 or
  N×1 tiled panels.
- **User is on Windows for the actual data (Google Drive mount `G:\`);
  edits/git happen on Windows too via MATLAB Command Window.** macOS
  SSD copy exists on `/Volumes/PortableSSD/processing_new` but is a
  workspace / development mirror.
- The user prefers TERSE responses and doesn't want emojis in code.
- Ignore the "task tools not used recently" system-reminder — user
  never uses TaskCreate and the flow is small and linear.

## 11. Quick-start for the next session

```matlab
% Recover
clear functions
close all force

% Diagnose data
diagnoseMetricCoverage        % file existence per (metric, source)

% Main plotters, in order of increasing derived-ness:
plotAveragedMetrics            % scalar averages, baseline vs recovery boxes
plotWindowedMetrics            % (rec-base)/base per window
plotWindowedViolins            % same data, box+violin
plotTimeTraces                 % mean±SEM traces per (metric, condition)
plotSynergyHeatmaps            % 3x3 M×E synergy per metric
plotRecoveryDynamics           % 4 section-5 features on normalised trace

% Caches (all at project root, gitignored):
%   gemsplots_queue.mat            — queue state
%   gemsplots_metrics_cache.mat    — scalars (plotAveragedMetrics)
%   gemsplots_series_cache.mat     — shared full series (everything else)
```

Delete a cache to force a full reload (e.g. after re-running
`batch_process` on a file, though mtime invalidation should handle that
automatically).

---

End of handoff. Ping @AndreaEBiju on GitHub if in doubt about anything above.
