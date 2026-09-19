# disttrait v0.10 package status — 2026-09-19

Status: validation candidate with explicit observation-process / MNAR stress testing.

Package: packages/disttrait/

Candidate version: 0.10.0

## 1. New in v0.10

v0.10 adds a benchmark that separates two different claims:

1. calibration conditional on the observed sample;
2. recovery of the latent biological process under outcome-dependent observation.

The package had already shown calibration under MCAR observation loss. v0.10 asks what happens when the probability that a row is observed depends on the trait, local position, or their interaction.

Canonical benchmark:

packages/disttrait/benchmarks/mnar_observation_surface.py

Frozen outputs:

- results/disttrait_mnar_observation_v0_10_20260919/cells.csv
- results/disttrait_mnar_observation_v0_10_20260919/result.json

Receipt checker:

packages/disttrait/scripts/check_mnar_observation_receipt.py

Dedicated workflow:

.github/workflows/disttrait-mnar-observation.yml

## 2. Latent biological null

Every simulated species has:

- a species-specific geographic centre;
- a species-specific baseline binary-trait frequency linked to that centre;
- 100 latent individual observations;
- local positions around the species centre;
- no within-species effect of local position on trait state.

Thus the benchmark deliberately retains strong between-species geographic turnover while the latent within-species biological association is null.

## 3. Observation mechanisms

Before analysis, observations are retained under one of four mechanisms.

### MCAR

Observation probability is independent of trait and local position.

### Trait-only selection

Observation probability depends on trait state but not position.

### Position-only selection

Observation probability depends on local position but not trait state.

### Joint trait-by-position selection

Observation probability depends on the interaction between trait state and local position.

The joint mechanism can make trait and position associated among observed rows even though they were independent in the latent biological process.

Two selection strengths are crossed with the four mechanisms: 0.8 and 1.6.

Each of the 8 cells contains 40 replicate worlds and all cells retain 40/40 evaluable worlds.

## 4. Compared inference

The benchmark reports:

- naive pooled distance–trait analysis;
- equal-species matched-null spatial organization;
- species fixed-intercept logistic inference.

The matched-null statistic retains the frozen disttrait construction:

rho_i = Spearman(pairwise geographic distance, pairwise trait dissimilarity).

The randomization fixes the observed rows, positions and trait multiset and permutes trait assignments within species.

## 5. Result under non-joint observation

Across MCAR, trait-only and position-only selection:

- maximum equal-species matched-null rejection = 0.05;
- median equal-species rho remains close to zero;
- mean absolute observed within-species state-position rho remains about 0.10–0.11.

The species fixed-effect logistic model reaches a maximum rejection fraction of 0.10 in one 40-world trait-only cell, while being at or below 0.075 in the other non-joint cells.

These are finite simulation frequencies, not universal calibration guarantees.

## 6. Result under joint trait-by-position selection

Under the joint observation process:

- equal-species matched-null rejection = 0.80–1.00;
- species fixed-effect logistic rejection = 1.00;
- median equal-species spatial rho rises from 0.043 to 0.120 as selection strength increases;
- mean absolute observed within-species state-position rho rises to 0.21–0.32.

Thus both nonparametric and model-based species-conditioned analyses detect strong spatial structure in the observed sample.

## 7. Why this is not a conventional type-I error result

The latent biological process is null, but the joint observation process has changed the conditional distribution of the retained data.

After selection, trait and position are no longer exchangeable within species.

Therefore a test that rejects exchangeability in the observed rows is not miscalibrated merely because the pre-selection biological process was null.

The methodological distinction is:

> matched randomization can calibrate a statistic to the observed sampling geometry, but cannot reconstruct a counterfactual unobserved dataset when the observation process itself depends jointly on trait and location.

This applies to the fixed-effect logistic comparator as well: adding a parametric response model does not identify whether the observed trait-position association is biological or selection-induced.

## 8. Implication for opportunistic data

v0.10 sharpens the package claim.

Supported:

> species-conditioning protects against between-species geographic turnover, and matched randomization provides assumption-light conditional calibration under an exchangeability null for the observed rows.

Not supported:

> arbitrary outcome-dependent observation processes are solved by species-conditioning or permutation alone.

Applications should therefore use observation-process metadata, matched controls, observer restrictions, sampling design, sensitivity analyses or external validation when trait-dependent spatial observation is plausible.

## 9. Relation to FCP

This limitation is directly relevant to the original FCP design choices.

The flower-colour application used:

- outcome-blind acquisition rules;
- observer caps;
- fixed photo budgets;
- geographic maximin selection;
- location-blind image measurement;
- flower/background matched diagnostics;
- explicit missingness states.

Those controls do not prove absence of every MNAR mechanism. They reduce the scope for the most obvious trait-dependent acquisition and measurement feedbacks.

v0.10 therefore strengthens rather than weakens the paper-level claim boundary: conditional calibration and observation-process validity are separate problems.

## 10. Methods-paper state after v0.10

The validation stack now includes:

1. FCP/RGFCA exact algorithm-equivalence fixtures;
2. pooled-confounding benchmarks;
3. effect/imbalance/MCAR performance surfaces;
4. nonparametric and model-based comparators;
5. direction-heterogeneity benchmarks;
6. permutation-calibrated species-slope meta-analysis;
7. two external non-flower empirical transports;
8. an explicit MNAR observation-process stress test.

The highest-value remaining strengthening steps are now:

1. nonlinear within-species processes;
2. multivariate continuous traits;
3. standalone package repository/licence/release metadata.

MNAR is no longer an untested gap, but arbitrary MNAR robustness remains a hard nonclaim.
