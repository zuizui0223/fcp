# disttrait v0.7 package status — 2026-09-19

Status: validated candidate with continuous-trait direction-heterogeneity benchmark.

Package:

`packages/disttrait/`

Version:

`0.7.0`

## 1. What v0.7 adds

v0.7 retains the v0.6 package, synthetic validation stack and external street-tree transport, and adds an explicit benchmark for **direction heterogeneity among species**.

The benchmark asks a methods question that the earlier binary-logistic comparator did not address:

> What happens when a trait is spatially organized within each species, but the signed direction of change is not shared across species?

This situation is important because two species can show equally strong geographic trait differentiation while one increases and the other decreases along its local spatial axis.

Canonical script:

`packages/disttrait/benchmarks/direction_heterogeneity_surface.py`

Frozen outputs:

- `results/disttrait_direction_heterogeneity_v0_7_20260919/cells.csv`
- `results/disttrait_direction_heterogeneity_v0_7_20260919/result.json`

Dedicated benchmark workflow:

`.github/workflows/disttrait-direction-heterogeneity.yml`

## 2. Data-generating process

Every simulated species has:

- a species-specific geographic centre;
- a species-specific baseline continuous trait level linked to that centre, creating strong between-species geographic confounding;
- repeated local observations around the species centre;
- a local continuous trait effect of specified magnitude;
- Gaussian residual variation.

The signed local effect is reversed in a prespecified fraction of species.

The benchmark crosses:

- effect size: **0.0, 0.4, 0.8, 1.2**;
- reversal fraction: **0, 0.25, 0.5**;
- MCAR observation loss: **0%, 50%**;
- **40** replicate worlds per cell;
- **20** species per world;
- **20** observations/species before missingness;
- **39** matched spatial-null permutations.

There are 24 benchmark cells.

## 3. Compared estimands/methods

### Naive pooled pairwise analysis

All species are pooled before relating geographic distance to absolute trait difference.

This deliberately ignores species identity and is included as a confounded baseline.

### Equal-species matched-null spatial omnibus

Within each species:

`rho_i = Spearman(pairwise geographic distance, absolute pairwise trait difference)`.

Species contribute equally, and significance is calibrated against matched vertex-permutation null worlds.

Because the response is an absolute pairwise difference, the statistic measures **strength of geographic organization without requiring a shared signed direction**.

### Species-rho one-sample test

A one-sided one-sample t-test is applied to species-specific `rho_i` values.

### Species fixed-intercept common-slope model

A linear model with:

- species-specific intercepts;
- one common signed within-species slope.

The model is correctly aligned when all species share a direction. As reversal fraction increases, the common-slope assumption is deliberately violated.

## 4. Null calibration

Across the six null cells:

- naive pooled minimum false-positive fraction = **1.00**;
- equal-species matched-null maximum FPR = **0.05**;
- species-rho t-test maximum FPR = **0.10**;
- common-slope model maximum FPR = **0.10**.

Thus the matched-null omnibus remains within the nominal 0.05 ceiling across these null cells, whereas the two parametric summaries reach 0.10 in at least one 40-world cell.

This is a finite simulation result, not a universal calibration theorem.

## 5. Shared-direction result

At effect size 0.4 with reversal fraction 0:

- common-slope model detection = **1.00** under both missingness levels;
- matched-null detection = **0.45–0.80**.

This is expected. The common-slope model is correctly aligned with the data-generating process and uses signed observation-level information efficiently.

## 6. Direction-heterogeneity result

At effect size 0.4 with 50% reversal:

- common-slope detection = **0.05–0.10**;
- matched-null detection = **0.475–0.90**.

At effect size 0.8 with 50% reversal:

- common-slope detection = **0.05–0.125**;
- matched-null detection = **1.00** under both missingness levels.

The median fitted common slope approaches zero as opposing signed effects cancel, despite strong species-specific geographic organization.

## 7. Interpretation

The result should **not** be phrased as “matched-null beats regression”.

The supported interpretation is:

> Method performance depends on the estimand. A common signed slope is highly efficient when species share a directional response. A direction-invariant distance–dissimilarity estimand retains information when species are spatially organized but differ in the sign of that organization.

This clarifies the conceptual role of the disttrait/RGFCA-style species-conditioned statistic. It asks whether individual traits become more dissimilar with geographic separation **within species**, not whether all species share one signed geographic response.

## 8. Relation to FCP/RGFCA

This benchmark formalizes a distinction already present in the FCP programme.

The original RGFCA shared-geography question weakened because different species need not differentiate in the same places or directions. The current flower-colour paper instead retains:

- species-specific spatial organization;
- a recurrent phenotype-space axis;
- cross-species association between diversity amount and spatial organization.

v0.7 demonstrates in a generic continuous-trait setting why direction-specific model summaries and direction-invariant spatial-organization summaries should not be conflated.

## 9. What v0.7 does not show

The benchmark does not establish that:

- common-slope models are generally invalid;
- absolute trait difference is the optimal dissimilarity for every phenotype;
- matched permutation inference is always more powerful;
- species-specific directions are biologically random;
- all hierarchical models require a shared slope;
- a random-slope hierarchical model would behave like the common-slope comparator;
- the benchmark covers nonlinear, multivariate or MNAR processes.

A random-slope hierarchical model is an obvious future comparator and would target the direction-heterogeneity setting more appropriately than the deliberately constrained common-slope model.

## 10. Methods-paper state after v0.7

The validation stack now includes:

1. exact FCP/RGFCA algorithm-equivalence fixtures;
2. optional raw-artifact replay;
3. synthetic non-flower transport;
4. between-species confounding benchmark;
5. effect/imbalance/missingness performance surface;
6. nonparametric/species-summary comparator surface;
7. correctly specified binary model comparator;
8. external non-flower continuous-trait empirical transport;
9. **continuous direction-heterogeneity / estimand-alignment benchmark**.

The largest remaining model-comparison gap is now a flexible random-slope/hierarchical alternative rather than a simple common-slope model. Broader external empirical replication and multivariate/nonlinear trait transport remain desirable for a standalone methods publication.
