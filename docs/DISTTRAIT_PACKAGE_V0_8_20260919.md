# disttrait v0.8 package status — 2026-09-19

Status: validated candidate with permutation-calibrated species-specific signed-slope meta-analysis.

Package:

`packages/disttrait/`

Version:

`0.8.0`

## 1. New in v0.8

v0.8 extends the v0.7 direction-heterogeneity benchmark with a model-based
alternative that does not force every species to share one signed slope.

Public API:

- `species_slope_estimate`
- `random_effects_slope_summary`
- `random_effects_from_groups`
- `random_effects_slope_permutation_test`

For each species, an ordinary least-squares signed slope is estimated together
with its sampling variance. The species slopes are then summarized by a
DerSimonian-Laird random-effects model.

The key v0.8 addition is **matched within-species permutation calibration**.
Coordinates/predictors and the complete observed trait-value multiset are fixed
within every species; trait values are permuted only within species. Every null
world re-estimates:

- all species-specific signed slopes;
- all within-species slope variances;
- the random-effects mean;
- Cochran Q;
- tau-squared.

The calibrated mean and heterogeneity probabilities are combined with a
Bonferroni omnibus.

## 2. Why calibration was required

An exploratory preflight first used the usual large-sample random-effects
normal/chi-square calibration.

That preflight failed under small per-species sample sizes with 50% MCAR
missingness:

- analytic mean-effect maximum null rejection = **0.175**;
- analytic Q heterogeneity maximum null rejection = **0.30**;
- analytic combined omnibus maximum null rejection = **0.30**.

The failed preflight is retained in:

`docs/DISTTRAIT_V0_8_ANALYTIC_META_PREFLIGHT_DIAGNOSIS_20260919.md`

and is not presented as a validated comparator.

After matched within-species permutation calibration:

- calibrated mean-component maximum null rejection = **0.075**;
- calibrated heterogeneity-component maximum null rejection = **0.075**;
- calibrated combined omnibus maximum null rejection = **0.05**.

The component values are finite 40-world simulation estimates and are not a
universal calibration theorem. The final combined omnibus stays within the
nominal 0.05 ceiling across the six null cells in this benchmark.

## 3. Benchmark design

Canonical benchmark:

`packages/disttrait/benchmarks/random_slope_meta_surface.py`

Frozen outputs:

- `results/disttrait_random_slope_meta_v0_8_20260919/cells.csv`
- `results/disttrait_random_slope_meta_v0_8_20260919/result.json`

The benchmark uses the same direction-heterogeneity design as v0.7:

- effect size: **0, 0.4, 0.8, 1.2**;
- reversal fraction: **0, 0.25, 0.5**;
- MCAR missingness: **0%, 50%**;
- 40 replicate worlds per cell;
- 20 species per world;
- 20 observations per species before missingness;
- 39 spatial matched-null permutations;
- 99 signed-slope/meta matched permutations.

## 4. Shared-direction result

At effect size 0.4 and reversal fraction 0:

- common signed-slope model detection = **1.00**;
- calibrated slope-meta omnibus = **1.00**;
- direction-invariant distance-dissimilarity matched-null detection =
  **0.50–0.825**.

When all species share direction, signed-slope models use the information more
efficiently.

## 5. Direction-heterogeneity result

At effect size 0.4 and 50% reversal:

- common signed-slope detection = **0.025–0.05**;
- distance-dissimilarity matched-null detection = **0.40–0.90**;
- calibrated slope-heterogeneity detection = **0.75–1.00**;
- calibrated slope-meta omnibus detection = **0.525–1.00**.

At effect size 0.8 and 50% reversal:

- distance-dissimilarity matched-null detection = **1.00**;
- calibrated slope-meta omnibus detection = **1.00**.

Thus a flexible species-slope analysis can recover directional heterogeneity
that cancels out of a common signed slope.

## 6. Interpretation

v0.8 does not identify a single universally best model.

It distinguishes three estimands.

### Common signed slope

Asks whether species share one directional response.

It is efficient when this assumption is scientifically correct.

### Species-specific signed-slope meta-analysis

Separates:

- average signed response;
- among-species directional heterogeneity.

This is useful when direction itself matters and may vary among species.

### Distance-dissimilarity matched-null organization

Asks whether individual traits become more dissimilar with geographic
separation within species, regardless of which signed direction generates that
organization.

It does not require a shared directional coordinate system across species.

## 7. Methodological lesson from the failed preflight

The failure of the uncalibrated random-effects test is itself informative.

With small species-level samples, plugging estimated slope variances into
standard meta-analytic normal/chi-square reference distributions can be
anti-conservative.

The supported v0.8 statement is:

> Species-specific slope models can represent directional heterogeneity, but
> their final mean/heterogeneity tests should not automatically inherit
> large-sample calibration. When the scientific null is exchangeability of
> traits among observed positions within species, matched within-species
> randomization can calibrate the slope/meta statistics to the actual sampling
> geometry.

## 8. Validation

Dedicated workflow:

`.github/workflows/disttrait-random-slope-meta.yml`

The validated exploratory run before receipt freeze was:

- PR #48;
- head: `a8b8a6f3eb52809da1f1a42c752ed7432c188a40`;
- run: `35391979183`;
- job: `105752183293`;
- artifact: `10565558831`;
- conclusion: success.

The final workflow additionally regenerates all 24 cells and checks them against
the frozen v0.8 receipt.

## 9. Methods-paper state after v0.8

The comparison stack now includes:

1. naive pooled inference;
2. equal-species matched-null spatial organization;
3. pair-weighted matched-null inference;
4. species-rho summary tests;
5. correctly specified fixed-intercept logistic models;
6. continuous common-slope models;
7. **species-specific signed slopes with random-effects mean and heterogeneity**;
8. **matched permutation calibration of those slope/meta statistics**;
9. an external non-flower continuous-trait empirical transport.

The most important remaining extensions are broader empirical transport,
nonlinear/multivariate dissimilarities and observation processes that are
missing-not-at-random.
