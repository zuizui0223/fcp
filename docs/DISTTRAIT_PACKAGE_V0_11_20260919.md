# disttrait v0.11 package status — 2026-09-19

Status: validation candidate with nonlinear within-species estimand benchmark.

Package: packages/disttrait/

Candidate version: 0.11.0

## 1. New in v0.11

v0.11 adds a nonlinear continuous-trait benchmark using a centered quadratic local-position effect within species.

Canonical benchmark:

packages/disttrait/benchmarks/nonlinear_curvature_surface.py

Frozen outputs:

- results/disttrait_nonlinear_curvature_v0_11_20260919/cells.csv
- results/disttrait_nonlinear_curvature_v0_11_20260919/result.json

Receipt checker:

packages/disttrait/scripts/check_nonlinear_curvature_receipt.py

Dedicated workflow:

.github/workflows/disttrait-nonlinear-curvature.yml

## 2. Benchmark design

Every simulated species has:

- a species-specific geographic centre;
- a species-specific baseline continuous trait level linked to that centre;
- 24 local observations before missingness;
- a centered quadratic local-position response;
- Gaussian residual noise.

The benchmark crosses:

- effect size: 0, 0.4, 0.8, 1.2;
- curvature reversal fraction: 0, 0.25, 0.5;
- MCAR observation loss: 0%, 50%;
- 40 replicate worlds per cell;
- 20 species per world;
- 39 matched spatial-null permutations.

There are 24 cells.

## 3. Compared estimands

### Naive pooled pairwise analysis

All species are pooled before relating geographic distance to absolute trait difference.

### Equal-species matched-null distance-dissimilarity

Within species, the statistic is the Spearman association between geographic distance and absolute pairwise trait difference. Species contribute equally.

This estimand measures strength of spatial organization without requiring the sign of a response to be shared.

### Common linear signed slope

A species fixed-intercept linear model estimates one shared local-position slope.

For the symmetric quadratic data-generating process, this model is intentionally misspecified.

### Common quadratic curvature

A species fixed-intercept model includes local position and squared local position and tests one shared quadratic coefficient.

This model is correctly aligned when species share curvature direction.

## 4. Null behavior

Across the six effect-zero cells:

- minimum naive pooled rejection = 1.00;
- maximum matched-null rejection = 0.075;
- maximum common-linear rejection = 0.10;
- maximum common-quadratic rejection = 0.075.

These are observed 40-world cell frequencies, not universal type-I error guarantees. The 0.075 values correspond to three rejections among 40 worlds in the worst cells.

## 5. Shared-curvature result

At effect size 0.4 with no curvature reversal:

- common quadratic detection = 1.00 under both missingness levels;
- matched-null detection = 0.10–0.375;
- common linear detection = 0–0.025.

The correctly specified nonlinear model is therefore much more efficient when the shared curvature estimand is scientifically correct.

## 6. Curvature-heterogeneity result

At effect size 0.4 with 50% curvature reversal:

- common quadratic detection = 0.025–0.05;
- matched-null detection = 0.15–0.35.

At effect size 0.8 with 50% curvature reversal:

- common quadratic detection = 0–0.075;
- matched-null detection = 0.85–1.00.

The fitted common quadratic coefficient approaches zero because convex and concave species cancel, even though individual species retain strong geographic trait organization.

## 7. Interpretation

The benchmark does not show that nonlinear models are inferior.

It separates scientific targets:

- a common quadratic coefficient asks whether species share one curvature direction;
- a distance-dissimilarity statistic asks whether traits become more different with geographic separation, regardless of the sign of the underlying curve.

The supported statement is:

> correctly specified nonlinear response models are highly efficient when species share the modeled shape, while direction-invariant dissimilarity estimands retain information when nonlinear response signs differ among species.

This extends the v0.7 direction-heterogeneity result from linear signed slopes to a nonlinear process.

## 8. Relation to the methods paper

v0.11 closes the previously listed nonlinear within-species benchmark gap.

The package now has explicit evidence for:

- between-species confounding;
- observation imbalance and MCAR loss;
- model-based versus matched-null efficiency;
- signed direction heterogeneity;
- calibrated species-slope heterogeneity;
- joint MNAR observation-process limits;
- nonlinear curvature and curvature-direction heterogeneity;
- two external empirical transports.

## 9. Remaining strengthening gaps

The highest-value remaining additions are now:

1. a multivariate continuous-trait benchmark;
2. standalone package repository/licence/release metadata.

Arbitrary nonlinear functions remain outside the benchmark claim; v0.11 tests one interpretable quadratic family rather than every possible nonlinear process.
