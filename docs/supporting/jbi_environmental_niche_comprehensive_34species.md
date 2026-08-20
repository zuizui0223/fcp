# Comprehensive environmental-niche analysis for the frozen 34-species baseline

This note records the symmetric five-metric analysis used by the current Journal of Biogeography manuscript and prevents selective reporting of only the strongest climatic metric.

## Analysis scope

The production input is the checksum-locked file `data/frozen/frozen_34species_five_metric_dataset.csv`: 34 species, 25 families, 20 within-population and 14 among-population cases, each with at least 20 occupied climate cells. The response is `among` (1 = geographically structured flower-colour variation; 0 = within-population coexistence).

Each climatic metric is fitted separately as:

`among ~ metric_z + effort_z`

where `metric_z` is the standardized focal metric and `effort_z = z(log1p(n_climate_cells))`. Five metrics are analysed symmetrically: temperature breadth, moisture breadth, climatic heterogeneity, PCA dispersion and PCA hull area.

The active production workflow is `.github/workflows/34species-paper.yml`. It reads the committed frozen dataset directly and runs family-clustered GLMs, 9,999 label permutations, leave-one-family-out refits, collinearity diagnostics, OpenTree + Grafen `phyloglm`, V.PhyloMaker2 time-scaled `phyloglm`, CR2/Satterthwaite inference and design-based power/precision diagnostics.

## Family-clustered GLM, permutation and family deletion

Canonical species ordering is fixed before finite Monte Carlo permutation so results are invariant to caller row order.

| Metric | OR | 95% CI | Clustered Wald p | Permutation p | Holm Wald p | Holm permutation p | LOFO OR range |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 0.817 | 0.384–1.739 | 0.6000 | 0.6131 | 1.0000 | 1.0000 | 0.636–0.973 |
| Moisture breadth | 0.412 | 0.180–0.947 | 0.0368 | 0.0423 | 0.1840 | 0.2115 | 0.306–0.465 |
| Climatic heterogeneity | 0.681 | 0.294–1.577 | 0.3700 | 0.3567 | 1.0000 | 1.0000 | 0.492–0.782 |
| PCA dispersion | 0.712 | 0.306–1.660 | 0.4317 | 0.3859 | 1.0000 | 1.0000 | 0.525–0.813 |
| PCA hull area | 0.577 | 0.312–1.067 | 0.0797 | 0.2372 | 0.3189 | 0.9488 | 0.489–0.671 |

All five point estimates are below one and every leave-one-family-out estimate remains below one. Moisture breadth has the largest negative effect estimate and is the only metric whose unadjusted family-clustered confidence interval excludes one. Neither its clustered nor permutation result remains below 0.05 after Holm adjustment across all five metrics.

Earlier workflow logs contain slightly different finite-draw permutation p-values because species row order was not explicitly canonicalized. That implementation detail has now been removed: the durable freeze and model helper define a single species order. Coefficients, clustered intervals, LOFO results and the biological interpretation did not change.

## Collinearity diagnostics

The five climatic metrics are fitted in separate models. Collinearity diagnostics therefore concern each focal standardized metric and standardized occurrence-effort covariate.

| Metric | Predictor correlation | Max VIF | Condition number |
|---|---:|---:|---:|
| Temperature breadth | 0.329 | 1.121 | 1.407 |
| Moisture breadth | 0.298 | 1.097 | 1.360 |
| Climatic heterogeneity | 0.412 | 1.204 | 1.549 |
| PCA dispersion | 0.295 | 1.095 | 1.355 |
| PCA hull area | 0.541 | 1.415 | 1.833 |

These values do not indicate problematic predictor collinearity.

## OpenTree + Grafen phylogenetic logistic sensitivity

OpenTree matching retained 30 species. One hundred polytomy-resolution fits completed for every metric and all fitted metric coefficients were negative.

| Metric | Species | Median OR | Median 95% CI | Median p | Holm median p | Negative fraction |
|---|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 30 | 0.908 | 0.436–1.893 | 0.7974 | 1.0000 | 1.00 |
| Moisture breadth | 30 | 0.573 | 0.234–1.403 | 0.2227 | 1.0000 | 1.00 |
| Climatic heterogeneity | 30 | 0.726 | 0.322–1.636 | 0.4398 | 1.0000 | 1.00 |
| PCA dispersion | 30 | 0.819 | 0.387–1.731 | 0.6006 | 1.0000 | 1.00 |
| PCA hull area | 30 | 0.654 | 0.261–1.641 | 0.3653 | 1.0000 | 1.00 |

## Time-scaled V.PhyloMaker2 sensitivity

All 34 species are retained under placement scenarios S1–S3. Point estimates remain below one for all five metrics in every scenario. S1 and S3 are numerically identical; S2 differs only slightly.

| Metric | S1 OR (95% CI; p) | S2 OR (95% CI; p) | S3 OR (95% CI; p) |
|---|---|---|---|
| Temperature breadth | 0.827 (0.385–1.777; 0.6258) | 0.838 (0.393–1.789; 0.6484) | 0.827 (0.385–1.777; 0.6258) |
| Moisture breadth | 0.448 (0.167–1.206; 0.1121) | 0.454 (0.171–1.211; 0.1147) | 0.448 (0.167–1.206; 0.1121) |
| Climatic heterogeneity | 0.670 (0.290–1.550; 0.3491) | 0.679 (0.295–1.565; 0.3635) | 0.670 (0.290–1.550; 0.3491) |
| PCA dispersion | 0.698 (0.319–1.527; 0.3681) | 0.707 (0.324–1.542; 0.3838) | 0.698 (0.319–1.527; 0.3681) |
| PCA hull area | 0.599 (0.232–1.545; 0.2892) | 0.612 (0.242–1.549; 0.3002) | 0.599 (0.232–1.545; 0.2892) |

For moisture breadth, Holm-adjusted p-values are 0.5603 for S1/S3 and 0.5736 for S2. Phylogenetic correction preserves the direction and approximate magnitude of the moisture estimate while increasing inferential uncertainty; all phylogenetic confidence intervals include one.

## Interpretation

The data do not establish moisture breadth as a unique driver. Rather, all five point estimates indicate lower odds of among-population organization with broader sampled occupied climatic niches, with the largest observed contrast along moisture. Direction is stable to family deletion and phylogenetic treatment, while multiplicity and phylogenetic uncertainty argue for effect-size and directional-consistency language rather than a binary significance claim. Occupied climate is not physiological tolerance, the analysis is species-level rather than morph-level, and the association is not causal.
