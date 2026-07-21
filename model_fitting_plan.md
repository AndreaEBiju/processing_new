# Layer 1–2 First-Pass Model: Final Equations and Fitting Plan

**Status: first-pass, deliberately incomplete.** This document consolidates everything established through iterative critique. It is not a finished model — it is the most defensible starting point given current data, with every known gap flagged rather than hidden. Read the "Known limitations" section before building anything from this.

---

## 1. Scope of this model

- Predicts: mean firing rate, CV2, FWHM — all three, from one shared parameter set, at the vagal cuff (RVN/LVN).
- Inputs: `u_M` (mechanical frequency: 10, 50, 100 Hz), `u_E` (electrical frequency: 10, 100, 1000 Hz).
- Explicitly out of scope: Layer 3 (NTS), Layer 4 (cardiac/respiratory output), amplitude dependence (`A_M`, `A_E` fixed in current data), any time-domain/nonstationary behavior within the 2-minute stimulation train.
- This is Layer 1+2 only.

---

## 2. Final equations

### 2.1 Shared building blocks

Electrical entrainment (saturating, no threshold — established from direct E-alone vagal neurogram data, NOT from HRV-based literature thresholds):

```
h(u_E) = u_E / (u_E + u_half)
```

C-fibre following, term-c/term-1 corner (degrades at high u_E — explains why E1000 doesn't dominate despite highest raw rate):

```
Phi_c(u_E) = f_max_c / (f_max_c + u_E)
```

C-fibre following, term-2's OWN corner (distinct from f_max_c; f_max_c2 < f_max_c so that
term_2 is relatively larger at low u_E and relatively smaller at high u_E than term_1 —
this is what actually produces the moving-optimum property, NOT f2 alone; see 2.2):

```
Phi_c2(u_E) = f_max_c2 / (f_max_c2 + u_E)
```

**NOTE on symbolic normalization:** `w_c` and `lambda_max_cap` are NEVER separable —
they multiply together in every place either appears, so no amount of data could ever
distinguish them. They are collapsed into a single parameter `W_c = w_c * lambda_max_cap`
throughout this document and the implementation. Do not reintroduce them as separate
parameters.

### 2.2 Mechanical-coupling functions (two required — see Section 3 for why)

Bandpass, tuned toward the M50xE100 hotspot:

```
g1(u_M) = [u_M / (u_M + f_lo)] * [f_hi / (u_M + f_hi)]
```

Monotonic saturating, tuned toward the M100xE10 hotspot:

```
g2(u_M) = u_M / (u_M + u_M_half)
```

Electrical-side partners for each (reuse existing blocks — no new functions):

```
f1(u_E) = h(u_E) * Phi_c(u_E)       # peaks at moderate u_E, pairs with g1
f2(u_E) = h(u_E) * Phi_c2(u_E)      # ALSO gated by h(u_E) -- see fix below.
                                      # Relatively favors low (nonzero) u_E vs f1,
                                      # via its own corner f_max_c2, pairs with g2
```

**CORRECTED from an earlier draft: `f2(u_E) = Phi_c(u_E)` alone was WRONG.**
`Phi_c(0) = 1`, not 0, so a term built from `Phi_c` alone does not vanish at `u_E = 0`
— it directly contradicts the "no standalone mechanical term" requirement (L1) and
would fail the `r_vagus(u_M, 0) == 0` unit test. The fix reintroduces the `h(u_E)`
entrainment gate (physically required: with zero electrical drive there is nothing
to entrain) while keeping a SEPARATE following corner `f_max_c2` so that `f1` and `f2`
remain non-proportional to each other (required for the moving-optimum property from
Section 3 — if they shared the same corner, term_1 and term_2 would just be rescaled
copies of one another and the model would collapse back to the rank-1 case that
cannot reproduce the hotspot shift).

### 2.3 Compound afferent signal (the core object — everything else derives from this)

```
r_vagus = W_c  * Phi_c(u_E) * h(u_E)
        + w_1  * f1(u_E) * g1(u_M)
        + w_2  * f2(u_E) * g2(u_M)
```

**No standalone mechanical-only term.** Every term above now correctly vanishes when
`u_E = 0`: term_c and term_1 via the shared `h(u_E)` factor, term_2 via its own `h(u_E)`
factor introduced in the 2.2 fix. This automatically matches the observed flat
mechanical-alone data — a *parsimony choice under identifiability constraints*, not a
claim that the true standalone mechanical effect is proven zero (see Limitation L1).

### 2.4 Derived observables (the three fit targets)

**Mean firing rate** — direct output of the compound signal:

```
rate_hat(u_M, u_E) = r_vagus(u_M, u_E)
```

**CV2 and FWHM** — NOT separate equations. Both derive from the *relative proportions* of the three additive terms in `r_vagus`, since different terms represent contributions with different underlying temporal/waveform character (cap-channel-dominated vs. g1-dominated vs. g2-dominated regimes).

Define the three term magnitudes and their normalized shares:

```
T_cap = W_c  * Phi_c(u_E) * h(u_E)
T_1   = w_1  * f1(u_E) * g1(u_M)
T_2   = w_2  * f2(u_E) * g2(u_M)

total = T_cap + T_1 + T_2
p_cap = T_cap / total
p_1   = T_1   / total
p_2   = T_2   / total
```

CV2 and FWHM as weighted combinations of per-channel reference values (`CV2_cap`, `CV2_1`, `CV2_2`, `FWHM_cap`, `FWHM_1`, `FWHM_2` — six ADDITIONAL free parameters representing the "pure" CV2/FWHM value if each term totally dominated):

```
CV2_hat(u_M, u_E)  = p_cap * CV2_cap  + p_1 * CV2_1  + p_2 * CV2_2
FWHM_hat(u_M, u_E) = p_cap * FWHM_cap + p_1 * FWHM_1 + p_2 * FWHM_2
```

**This is a first-pass linear-mixture approximation, not a derivation from spike-train statistics** (see Limitation L2).

### 2.5 No random effect for this first pass

All animals share identical parameters (fully pooled model):

```
r_vagus_observed(u_M, u_E, animal=a) = r_vagus(u_M, u_E)     # same for every animal a
```

**A per-animal random effect was considered and deliberately dropped for this pass** (see Limitation L3 — this replaces an earlier draft that included one). The original justification — absorbing motor-to-tissue coupling gain into a per-animal threshold shift — no longer applies, since amplitude dependence (and the threshold it would have shifted) was already collapsed out in Section 2.1–2.3. A generic catch-all random scalar could be reintroduced, but given the parameter budget is already tight relative to available conditions (see L9), any animal-to-animal variability will show up as unexplained residual noise for now. Revisit only if residuals show a clear animal-clustered pattern once real fit results exist — and if so, check first whether that variability looks common-mode across the cap and mechanical terms (consistent with a shared scalar) or channel-specific (e.g., driven separately by electrode contact quality vs. motor contact quality), since a single shared multiplicative term is only the right shape for the former.

---

## 3. Why this specific structure (one-paragraph justification per decision, for your own reference when explaining this to Claude or your committee)

- **Two mechanical-coupling terms, not one:** the M50xE100 and M100xE10 hotspots require the model's preferred `u_M` to shift depending on `u_E`. A single separable term `f(u_E)*g(u_M)` cannot do this — its optimal `u_M` is fixed regardless of `u_E` (provable by taking d/du_M = 0; `f(u_E)` cancels out of that equation). Two terms, each paired with a *different* `f(u_E)`, allow the ratio of their electrical weightings to shift with `u_E`, which is what actually moves the optimum.
- **No standalone mechanical term:** would predict a `u_M`-dependent effect at `u_E=0`, and while the mechanical-alone data is *not statistically conclusive* about a zero effect (small n, noisy baseline-ratio metric), the model is already at/near its parameter budget relative to available conditions (9 params vs. 9 M×E combinations) — this exclusion is about affordability under current data, not a claim the true effect is zero.
- **CV2/FWHM from shared term-shares, not independent equations:** avoids fitting three unrelated models to three metrics that plausibly share the same underlying channel-mixture mechanism; keeps total parameter count from tripling.
- **Electrical entrainment has no threshold:** direct vagal cuff data shows real effects at E10, ruling out the ~19 Hz HRV-derived threshold from prior literature for this preparation.
- **C-fibre attribution (not A-fibre):** capsaicin abolition data (Peles 2003; Liu & Chen 2004) — establishes TRPV1+ C-fibre mediation for the electrically-driven channel specifically in this serosal-electrode/short-pulse preparation.

---

## 4. Known limitations — carry ALL of these into any code/writeup, do not silently drop

| ID | Limitation | Why it's there | What would resolve it |
|----|-----------|-----------------|------------------------|
| L1 | No standalone mechanical term | Parsimony under n=3-4, NOT a proven zero effect | More animals; fasted/fed protocol to reduce baseline noise |
| L2 | CV2/FWHM via linear term-share mixture | Placeholder; not derived from actual spike-train/ISI statistics | Explicit point-process/spike-train model of the three channels |
| L3 | No per-animal random effect (fully pooled model) | Original coupling-gain justification no longer applies once amplitude was collapsed out (Sec 2.5); parameter budget already tight (L9) | Check residuals for animal-clustered patterns once real fit exists; if present, determine common-mode vs. channel-specific before reintroducing |
| L4 | `g1`, `g2` mechanical shapes are NOT independently validated at receptor level | Justified only by the interaction pattern (hotspot shift), not by controlled single-modality mechanical dose-response data | Bench-top motor characterization; fasted/fed intragastric pressure manipulation; finer u_M sweep (e.g., 25/50/75/100 Hz). **See L17 for a concrete manifestation of this limitation observed in real fitting: `(u_M_half, w_2)` becomes poorly identified specifically because the tested u_M range (10-100 Hz) is too narrow to constrain `g2`'s saturation point when the optimizer pushes `u_M_half` well beyond it.** |
| L5 | IGLE/IMA/mucosal receptor identity NOT assigned to any term | No experimental tool currently separates them | Piezo2 knockout (mouse only); topical mucosal anesthesia; fasted/fed pressure manipulation |
| L6 | No nonstationarity/adaptation within the 2-min stimulation train | Explicitly deferred; C-fibre activity-dependent slowing literature suggests this is non-negligible at 2 min | Add time-dependent state variable once basic structure is validated |
| L7 | Electrical-to-mechanoreceptor cross-leak (`rho_E`) not included | Deferred; assumed 0 | Gadolinium + electrical-alone experiment (imperfect tool, but partial info) |
| L8 | Amplitude dependence (`A_M`, `A_E`) collapsed out entirely | Both fixed in current dataset — not identifiable | Future cohort with amplitude variation |
| L9 | Parameter count (~9) vs. distinct M×E conditions (9) is very thin | Near-zero residual degrees of freedom | Expand cohort (n and/or conditions); or fix some parameters from literature/independent data before this fit |
| L10 | `g1`/`g2` split is a rank-2 mathematical convenience representing ONE lumped, unresolved mechanoreceptor population — NOT two distinct receptor types | Avoid re-introducing an unsupported IGLE/IMA claim | N/A — this is a permanent framing point, not something to "resolve" |
| L11 | **Static (steady-state) equations are used PROVISIONALLY, not because dynamics were ruled out** | Time-resolved CV2/FWHM traces exist in principle (raw vagal recordings are available) but have not yet been computed/inspected | **REQUIRED CHECK before treating any fit from this plan as final:** once CV2/FWHM time traces are computed within the 2-min stimulation window, check whether they (a) plateau within the window -> static equation's steady-state is a reasonable target, safe to proceed as-is; or (b) show ongoing drift without settling (as the C-fibre activity-dependent slowing literature, L6, suggests is plausible) -> static equations are structurally wrong regardless of fit quality, and Layer 2 needs an accumulating-state dynamic (analogous to the FDD variable `D(t)` in Layer 3, but peripheral) instead of a fixed-point relaxation. Also check whether the CV2 sliding-window width used to compute the trace is itself wide enough to be manufacturing an apparent plateau — decouple windowing-smoothness from genuine physiological settling before concluding (a). |
| L12 | **RVN and LVN Stage 1 fits (`u_half`, `f_max_c`) share `(u_half, f_max_c)` across channels; only `W_c` is fit per-channel. UPGRADED: this is a PROVEN exact one-parameter symmetry of the Stage 1 objective, not merely an empirical ridge.** | Confirmed both empirically (independent per-channel fits landed in different basins with near-identical residuals) AND algebraically: `term_c(u_E) = W_c * f_max_c/(f_max_c+u_E) * u_E/(u_E+u_half)`. Swapping `u_half <-> f_max_c` leaves the denominator `(f_max_c+u_E)(u_E+u_half)` EXACTLY unchanged (same two factors, either order) and only rescales the numerator -- so `term_c` is EXACTLY invariant under the swap if simultaneously `W_c -> W_c * (f_max_c/u_half)`. This is an exact symmetry for every `u_E`, not an approximation. Synthetic ground-truth validation (Python, prior to MATLAB) confirmed all three predicted consequences precisely: (1) pooled-vs-independent resnorm changed by only +0.7% (the true difference is exactly zero; the residual gap is optimizer/numerical, not a real fit-quality difference between basins); (2) `W_c_LVN/W_c_RVN` recovered the true ratio (1.74 vs true 1.75) despite absolute `W_c` values being off by ~22x -- because BOTH channels share the SAME `u_half`/`f_max_c` and therefore get multiplied by the IDENTICAL compensating factor, which cancels exactly out of the ratio; (3) `sqrt(u_half*f_max_c)` recovered near-exactly (397 vs true 394) -- trivially invariant since it is symmetric in its two arguments by construction. The ~22x factor itself matched the algebraic prediction `f_max_c_true/u_half_true` (~1867/84.5 ≈ 22.1) almost exactly. **FURTHER CORROBORATED by the L13 resolution below: on real fit output, RVN and LVN independently (each using only its own interaction data) preferred the exact same specific candidate pair `(u_half≈84.9, f_max_c≈1799)` over the same alternative, by wide margins (RVN resnorm 0.139 vs 0.347; LVN resnorm 0.306 vs 0.447). This is NOT guaranteed by the L12 algebra alone -- the exact symmetry only says Stage 1 cannot distinguish the two candidates; it says nothing about whether two independent channels' interaction data would agree on which one is physically correct. That they did agree is real empirical support for the pooling assumption itself (shared fibre-level physiology across the bilateral vagus), not just a restatement of the symmetry.** | **Never report or interpret raw `u_half`/`f_max_c`/`W_c` individually from Stage 1 -- only `sqrt(u_half*f_max_c)` (location/shape invariant) and cross-channel `W_c` ratios are identifiable even with infinite noiseless data.** Concrete fix beyond post-hoc checking: reparametrize Stage 1 to fit `rho = sqrt(u_half*f_max_c)` and `kappa = f_max_c/u_half` (or `log(kappa)`, which flips sign rather than value under the swap -- may be a cleaner optimizer target) directly, plus one reference `W_c` and the cross-channel ratio, instead of fitting `u_half`/`f_max_c` separately and relying on a post-hoc invariant check. Not blocking for the current smoke-test pass, but should be done before raw Stage 1 values are ever reported as if meaningful on their own. |
| L13 | **RESOLVED: Stage 2 reliably and decisively breaks the Stage 1 `(u_half, f_max_c)` tie -- both channels independently select the same specific candidate.** | On real fit output: RVN and LVN, each using ONLY its own M×E interaction data, independently preferred the identical `(u_half≈84.9, f_max_c≈1799)` pair over the swapped alternative `(u_half≈1799, f_max_c≈84.9)` -- RVN resnorm 0.139 (selected) vs 0.347; LVN resnorm 0.306 (selected) vs 0.447; pooled/total-channel driver showed the rejected candidate was 141.8% worse. This is NOT the same finding as L12 -- L12 is an exact algebraic fact about Stage 1 alone (the two candidates are indistinguishable there); L13 is an empirical finding about Stage 2 (once `g1(u_M)`'s dependence on `u_M` is brought in via `term_1`, the tie breaks cleanly). Mechanistically: `term_1 = w_1*f1(u_E)*g1(u_M)` shares `u_half`/`f_max_c` with `term_c` via `f1`, but ALSO multiplies against `g1(u_M)`, which has no compensating degree of freedom under the swap -- this asymmetry is what allows Stage 2's interaction data to discriminate what Stage 1's E-alone data cannot. | **Architectural consequence, now required (see L15): the pipeline must not report/plot/save whichever candidate Stage 1's optimizer happened to converge to first -- it must run Stage 2 for both candidates and propagate the WINNING one to every downstream output.** This was previously only printed as a diagnostic; L15 specifies the required restructuring. |
| L15 | **Stage 1 → Stage 2 handoff must select and propagate the Stage-2-preferred `(u_half, f_max_c)` candidate, not whatever Stage 1's optimizer converged to.** | Direct consequence of L12 (Stage 1 cannot distinguish the two candidates) + L13 (Stage 2 can, decisively). Without this, summary tables/plots/pattern-check margins reflect an arbitrary, optimizer-initialization-dependent choice rather than the data-supported one -- despite the correct answer already being computable from fits the pipeline is already running. | Required pipeline change (see `matlab_implementation_instructions.md` Sections 4/5/8 for the exact restructuring): (1) compute the swapped candidate analytically from Stage 1's raw output (closed-form, no re-fit); (2) run Stage 2 for BOTH candidates, BOTH channels (this already happens for the L13 check -- no new fitting required, only using the results differently); (3) select the winner by TOTAL resnorm summed across channels (a joint/shared-parameter decision, not per-channel, since `u_half`/`f_max_c` are pooled); (4) promote the winning candidate's already-computed Stage 2 fits to be the official `thetaFull` per channel -- do not re-fit a third time; (5) ALSO verify per-channel agreement on the winner as an automated check every run, not just this one-off manual comparison -- if per-channel winners ever disagree, flag prominently as "channel disagreement -- pooling assumption in question" rather than silently defaulting to the total's answer; (6) add a minimum-margin threshold (e.g. 10-15%) below which the pipeline flags "AMBIGUOUS -- review manually" rather than auto-selecting, since current margins are decisive but future datasets may not be. |
| L16 | **`(f_lo, f_hi, w_1)` has the SAME exact swap symmetry as `(u_half, f_max_c, W_c)` (L12) -- CONFIRMED algebraically, and structurally UNBREAKABLE (no Stage-2-style tie-break exists for this one).** | Algebra: `g1(u_M) = u_M*f_hi / [(u_M+f_lo)(u_M+f_hi)]`. The denominator is the same two factors regardless of order, so swapping `f_lo<->f_hi` only rescales the numerator -- exactly compensated by `w_1 -> w_1*(f_hi/f_lo)`. Unlike L12/L13, there is NO analogous tie-break available: `f_lo`/`f_hi` appear ONLY inside `g1`, nowhere else in the model (no term_2-style analog with a differently-paired dependence to discriminate against), and Stage 2 is not pooled across channels the way Stage 1 is -- there is no second data source, pooled or otherwise, that could break this tie the way `term_2`'s independent `f_max_c2` broke L12's. This was ALREADY VISIBLE, unlabeled, in prior restart-spread output: the swapped-candidate RVN fit showed `f_lo [0.3557, 37.68]`, `f_hi [0.005314, 35.11]`, `w_1 [6.718, 4.533e+04]` -- wide spread across restarts at equivalent resnorm, while `u_M_half`, `f_max_c2`, `w_2` stayed tight. That was this exact symmetry manifesting as restart instability, not a convergence failure. **Interpretive consequence (distinct from the fitting-mechanics consequence): the model's original motivating narrative -- from the earlier LaTeX specification (`Model_Revised_Full_Specification_v3.pdf`), predating this fitting plan -- assigns `f_lo` to a mechanical recruitment threshold and `f_hi` SPECIFICALLY to Piezo2 channel inactivation kinetics (Moroni et al. 2018). That narrative is not restated in this document's Section 2.2, which only gives the bare equation, but it is the reason this specific bandpass functional form was chosen and it is worth revisiting given this finding. The fit cannot confirm or refute that specific assignment -- it can only confirm a bandpass shape exists with peak `sqrt(f_lo*f_hi)`. Which fitted number is "the recruitment corner" vs "the inactivation corner" is not data-driven, exactly as `u_half` vs `f_max_c`'s physical assignments (entrainment vs following-limit) were not data-driven in L12.** No equations, predictions, resnorm, or pattern-check outcomes are affected -- `r_vagus` is identical under the swap by construction. | **Fix by construction (reparametrize), not by search (no data-driven tie-break exists to find).** Fit `rho_1 = sqrt(f_lo*f_hi)` (peak location, already exactly invariant -- make this the PRIMARY fit variable) and `kappa_1 = 0.5*ln(f_hi/f_lo)` directly, instead of `f_lo`/`f_hi` separately. The swap corresponds exactly to `kappa_1 -> -kappa_1`; constrain `kappa_1 >= 0` (report sign as a fixed convention, not a finding) so the fit lands in one place regardless of restart initialization. Apply the identical fix to `(u_half, f_max_c)` in L12/Section 4 for consistency -- these are structurally the same problem solved twice with two different mechanisms otherwise. **Never claim the fit confirms which specific parameter corresponds to which named physical mechanism (recruitment vs. Piezo2 inactivation, or entrainment vs. following-limit) -- only claim the fit confirms a bandpass/tuning shape with a given peak location and width exists.** Testing the specific mechanistic attribution would require an independent external estimate (e.g. a literature-derived predicted Piezo2 inactivation timescale) to check against, not just the r_vagus fit alone. |
| L17 | **`(u_M_half, w_2)` shows a near-degenerate ridge specifically when the optimizer's `u_M_half` lands well beyond the observed `u_M` range (10-100 Hz) -- CONFIRMED numerically, but this is an ASYMPTOTIC ridge, NOT an exact symmetry like L12/L16, and is a direct manifestation of the already-documented L4 limitation (coarse 3-point u_M sweep, no independent mechanical validation).** | Algebra: `g2(u_M) = u_M/(u_M+u_M_half)`. When `u_M_half >> max(u_M tested) = 100`, a first-order expansion gives `g2(u_M) ~ u_M/u_M_half`, so `term_2 ~ [w_2/u_M_half] * f2(u_E) * u_M` -- only the RATIO `kappa2 = w_2/u_M_half` is identified in this regime, not the two parameters individually. This is fundamentally DIFFERENT IN KIND from L12/L16: it is a statement about where the DATA stops constraining the curve (only 3 discrete u_M points, none near or beyond `u_M_half`'s fitted value), not an algebraic identity holding for all parameter values -- at `u_M_half` = 500 (5x the max tested u_M) the linear approximation already has ~4% error; by 1000 (10x) it is under ~1%; but for `u_M_half` within or near the tested 10-100 Hz range, `w_2` and `u_M_half` ARE well-separated by the curvature of `g2`. Did NOT break this run's pattern checks or basin selection (L13/L15) -- lower severity than L12/L16, no urgent structural fix needed. | **Numerical safeguard, not reparametrization or selection -- no exact substitution removes an asymptotic ridge the way it removes an exact symmetry.** Add a post-fit check: if fitted `u_M_half` exceeds ~5x the max tested `u_M` (i.e. > 500 Hz given the current 10-100 Hz range), flag "asymptotic ridge risk -- `u_M_half`/`w_2` individually unreliable in this regime" and report `kappa2 = w_2/u_M_half` alongside the raw values as the more robust quantity. Do NOT hard-bound `u_M_half` in the optimizer to suppress this -- a genuinely large true `u_M_half` (i.e. `g2` really is close to linear across the whole tested range) is a real, meaningful finding, not an error to be constrained away; only flag/monitor, don't silence. Monitor every run, not just this one -- a different restart or dataset could land here even when this run didn't. The actual resolution remains what L4 already specifies: the finer u_M sweep (25/50/75/100 Hz, or wider) would directly constrain `g2`'s saturation point and shrink or eliminate this ridge, unlike L12/L16 which no amount of `u_M`/`u_E` data could ever resolve. |

---

## 5. Fitting plan

### 5.1 Data preparation

```
INPUT: per-animal, per-condition trial-level data for all 15 M×E combinations
       (9 combined M×E, 3 M-alone, 3 E-alone; M-alone/E-alone needed for
       partial validation even though not in the core r_vagus fit target)

FOR each animal a, each condition (u_M, u_E):
    compute trial-level: mean_rate, CV2(RVN), CV2(LVN), FWHM(RVN), FWHM(LVN)
    (decide RVN/LVN: fit separately first, pool only if parameter estimates
     are consistent between sides — do NOT pool by default)

Structure as long-format table:
    columns = [animal_id, u_M, u_E, side(RVN/LVN), metric_type, value]
```

### 5.2 Parameter vector (first pass, RATE ONLY — see 5.4 for full joint fit)

**UPDATED per L12: `u_half`/`f_max_c` are now SHARED across RVN/LVN; everything
else is fit per-channel.** The parameter count below is no longer "one theta
vector for the model" — it is split by what's shared vs. per-channel.

```
theta_shared    = [u_half, f_max_c]                  # Stage 1, SHARED across both channels

theta_perchan   = [W_c,                                # Stage 1, per-channel
                    f_lo, f_hi,                         # Stage 2, per-channel (g1 shape)
                    u_M_half,                            # Stage 2, per-channel (g2 shape)
                    f_max_c2,                            # Stage 2, per-channel
                    w_1, w_2]                            # Stage 2, per-channel

n_params_shared           = 2
n_params_perchan           = 7
n_params_total_2channels   = 2 + 2*7 = 16   # 2 shared + 7 x 2 channels

n_conditions_stage1_pooled = 6    # 3 E-alone conditions x 2 channels
n_conditions_stage2_perchan = 9   # combined M×E, PER channel (not pooled)
```

**Note on parameter count:** the 9-parameter figure from earlier drafts of this plan
described a SINGLE channel's fit in isolation. With two channels and the L12
pooling decision, the honest total across both channels is 16 (not 18, since
`u_half`/`f_max_c` are no longer duplicated per channel) — but the more relevant
comparison is per-stage: Stage 1 now fits 4 params (`u_half`, `f_max_c`, `W_c_RVN`,
`W_c_LVN`) against 6 data-generating conditions (better-conditioned than the old
3-against-3 that produced the basin-splitting in the first place), while Stage 2
still fits 7 params against 9 conditions PER CHANNEL (unchanged — no pooling
argument yet exists for these, see Section 5.3).

Separately, and NOT to be confused with this: merging the old separate `w_c` and
`lambda_max_cap` into `W_c` (Section 2.1) is a genuine, provable symbolic
degeneracy — `W_c` was never separable from `lambda_max_cap` under ANY amount of
data. The `(u_half, f_max_c)` ridge (L12) is a DIFFERENT kind of problem — a
statistical weak-identifiability that showed up empirically with current data
and is addressed here by pooling, not a permanent structural degeneracy. Do not
conflate the two when writing this up.

**Flag immediately: Stage 2's 7 params vs. 9 conditions, per channel, still leaves
very little residual degrees of freedom.** Two options, present both to Claude
when building code:

- **Option A (recommended first):** Stage 1 pooled as above, THEN Stage 2 run
  independently per channel using that channel's own `W_c` but the shared
  `u_half`/`f_max_c`. This is now the DEFAULT procedure per L12, not just one of
  two equally-weighted options — Option B below is retained only as a validation
  check, not a live alternative.
- **Option B (validation only — see L12):** also fit `u_half`/`f_max_c`
  independently per channel (the ORIGINAL, now-superseded procedure) and compare
  total residual against the pooled version. Report both. If pooling increases
  residual sharply, that is evidence against the pooling assumption and should be
  flagged prominently, not buried.

**Decide and document which option was used — do not let Claude silently pick one.**

### 5.3 Optimization procedure (rate-only first pass)

```
1. Fit E-alone sub-model, POOLED ACROSS BOTH CHANNELS (see L12 -- confirmed
   empirically necessary, not just a theoretical option):
   theta_E_shared = [u_half, f_max_c]              # SHARED across RVN and LVN
   theta_E_perchan = [W_c_RVN, W_c_LVN]             # SEPARATE per channel

   objective: minimize sum_over_animals_and_BOTH_channels[
       (observed_rate - W_c_{channel}*Phi_c(u_E; f_max_c)*h(u_E; u_half))^2 ]
   over the 3 E-alone conditions x 2 channels x n_animals trials
   -> fix u_half, f_max_c (shared) and W_c_RVN, W_c_LVN (per-channel) for step 2

   RATIONALE (L12): independent per-channel fits of (u_half, f_max_c) landed in
   different basins with near-identical residuals -- direct empirical evidence
   of the ridge already anticipated in Section 5.5 item 3. u_half/f_max_c
   describe fibre-level physiology with no expected left/right difference;
   W_c absorbs per-electrode coupling, which plausibly DOES differ by channel.

   VALIDATION (required, not optional): also fit the two channels independently
   (as in the original single-channel procedure) and compare total residual
   against the pooled fit. Report both. If pooling increases residual sharply,
   flag this prominently -- do not silently prefer the pooled result.

2. Fit interaction parameters, PER CHANNEL (u_M x u_E structure may genuinely
   differ by channel -- no pooling assumption made here, unlike step 1):
   theta_interaction = [f_lo, f_hi, u_M_half, f_max_c2, w_1, w_2]
   (u_half, f_max_c, W_c_{channel} held FIXED from step 1, using that channel's
    own W_c but the SHARED u_half/f_max_c)
   objective: minimize sum_over_animals[
       (observed_rate - rate_hat(u_M, u_E))^2 ] over 9 combined conditions,
       for ONE channel at a time
```

**On both channels separately:** step 2 is deliberately NOT pooled the way step
1 is — there is no equivalent physiological argument yet for why `f_lo, f_hi,
u_M_half, f_max_c2, w_1, w_2` should be shared between RVN and LVN, so each
channel gets its own full interaction fit. Revisit this if the same kind of
basin-splitting shows up here too.

Model is fully pooled across ANIMALS (Section 2.5) — no per-animal parameters
to fit anywhere in this procedure. Standard nonlinear least squares
(e.g., `lsqnonlin`) is sufficient; no mixed-effects machinery needed for this
pass. Animal identity is only used post-hoc, in residual diagnostics
(Section 5.5), to check whether pooling ACROSS ANIMALS was reasonable — this is
a separate question from the channel-pooling in step 1 above.

### 5.4 CV2/FWHM extension (second pass — do NOT attempt until rate fit is validated)


```
ADDITIONAL parameters:
    CV2_cap, CV2_1, CV2_2        (3 params)
    FWHM_cap, FWHM_1, FWHM_2     (3 params)
n_params_total = 9 + 6 = 15   # against 9 conditions x 3 metrics = 27 data "points"
                                 # (still thin once you account for correlation
                                 #  between the 3 metrics at each condition —
                                 #  they are NOT independent observations)

Fit jointly: minimize weighted sum of squared residuals across all three
metrics simultaneously (rate, CV2, FWHM), holding theta_rate parameters
shared with step 5.3, only adding the 6 new linear-mixture coefficients.

Weight residuals by inverse variance per metric (rate, CV2, FWHM are on
different scales and likely different noise levels) — do NOT use raw
unweighted squared error across metrics of different units.
```

### 5.5 Model evaluation

```
1. Residual diagnostics (rate fit):
   - Plot observed vs. predicted rate per condition, per animal
   - Check residuals vs. u_M, u_E separately for structure
     (structured residuals -> missing term; random scatter -> ok for this pass)
   - Specifically check: does the model correctly reproduce
     M50xE100 > M100xE100 AND M100xE10 > M50xE10?
     (this is the exact pattern that justified the two-term structure —
      if the fit doesn't reproduce it, the structure has failed at its
      one required job)

2. Leave-one-animal-out (LOAO) cross-validation:
   FOR each animal a:
       refit theta on all animals EXCEPT a
       predict rate/CV2/FWHM for animal a's conditions
       record prediction error
   Report mean +/- range of LOAO error across animals
   (expect this to be poor given n=3-4 and near-zero residual df —
    report honestly, do not oversell)

3. Parameter identifiability check (DO THIS BEFORE trusting any fitted value):
   - `(u_half, f_max_c)`: NOT merely a ridge -- PROVEN exact one-parameter
     symmetry of the Stage 1 objective (see L12 for the full algebraic proof,
     confirmed by synthetic ground-truth validation). Swapping `u_half<->f_max_c`
     and rescaling `W_c -> W_c*(f_max_c/u_half)` leaves the Stage 1 fit EXACTLY
     unchanged. DO NOT report or interpret raw `u_half`/`f_max_c`/`W_c` values
     individually -- report ONLY `sqrt(u_half*f_max_c)` and cross-channel `W_c`
     ratios, which are exactly invariant under this symmetry and are the only
     Stage 1 quantities that are identifiable even with infinite noiseless data.
   - **CONFIRMED (L16), and unlike `(u_half, f_max_c)`, UNBREAKABLE:** `(f_lo, f_hi, w_1)`
     has the identical exact swap symmetry, but no analog of `term_2`/`f_max_c2`
     exists to discriminate against it, and Stage 2 has no pooling across channels
     to lean on either. Reparametrize as `rho_1 = sqrt(f_lo*f_hi)`,
     `kappa_1 = 0.5*ln(f_hi/f_lo)` with `kappa_1` sign-constrained -- do not expect
     a data-driven tie-break the way L13 found for Stage 1.
   - **RESOLVED (L13):** the `(u_half, f_max_c)` symmetry does NOT survive into
     Stage 2 -- `term_1`'s additional `g1(u_M)` dependence breaks the tie
     decisively. Confirmed on real fit output: RVN and LVN independently selected
     the same specific candidate pair by wide margins (RVN 0.139 vs 0.347; LVN
     0.306 vs 0.447 resnorm). **Required consequence (L15): the pipeline must
     select and propagate the Stage-2-winning candidate to every downstream
     output** -- see L15 for the exact required restructuring. Do not report
     results from whichever candidate Stage 1's optimizer happened to converge
     to first.

4. Compare against the SPECIFIC failure mode expected if this model is
   too simple: since CV2/FWHM use a linear-mixture approximation (L2),
   check whether residuals for CV2/FWHM are systematically worse at
   conditions where the term-shares (p_cap, p_1, p_2) are similar to
   each other (i.e., no term dominates) -- this is where a linear
   mixture is most likely to break down.

5. Explicitly report all items from Section 4 (Known Limitations) as
   "not addressed by this fit" in any writeup of results -- do not let
   a good-looking fit imply these are resolved.

6. BEFORE finalizing any conclusions from this static-equation fit,
   compute time-resolved CV2/FWHM within the 2-min stimulation window
   from the raw vagal recordings (these exist but have not yet been
   examined this way) and check plateau-vs-drift per L11. A good
   steady-state fit does NOT substitute for this check -- a model can
   fit condition-level means well while still being structurally wrong
   about the dynamics that produced them.
```

---

## 6. What to tell Claude when building this

- Implement Section 2 equations exactly as written; do not "improve" or add terms without flagging it back to me first.
- Implement Section 5.3 as TWO separate fitting stages (Option A), clearly labeled, with the fixed-then-interaction structure explicit in the code, not hidden inside a single black-box optimizer call.
- Do NOT implement CV2/FWHM (Section 5.4) until I confirm the rate-only fit (5.3) has been reviewed.
- Build the identifiability checks (5.5, item 3) as a required step, not optional — both `(u_half, f_max_c)` and `(f_lo, f_hi)` are CONFIRMED exact swap symmetries (L12, L16), not just anticipated ridges. `(u_half, f_max_c)` gets resolved by Stage 2 (L13/L15 — select the winning candidate). `(f_lo, f_hi, w_1)` has no such resolution available and must be fixed by reparametrization instead (L16 — fit `sqrt(f_lo*f_hi)` and a sign-constrained log-ratio directly, not `f_lo`/`f_hi` separately).
- Surface every one of the Section 4 limitations in code comments at the relevant point of implementation (e.g., comment above the "no standalone mechanical term" line pointing to L1).
- The static equations in Section 2 are used provisionally (L11). As a SEPARATE deliverable from the fit itself, write a short script that computes CV2 and FWHM in sliding time windows across the 2-minute stimulation period from the raw vagal recordings, and plots these traces per condition. This is a diagnostic, not part of the fitting pipeline — its only job is to let me visually check plateau vs. drift before I decide whether the static equations in this plan remain valid. Do not fold this into or make it a prerequisite for the Section 5.3 fit; run them in parallel and bring both back for review.
