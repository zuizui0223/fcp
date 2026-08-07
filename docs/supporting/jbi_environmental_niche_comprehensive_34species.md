# Comprehensive environmental-niche analysis for the frozen 34-species baseline

This document records the symmetric five-metric analysis requested for the frozen baseline-unambiguous dataset. It is intended to support the JBI manuscript revision and to prevent selective reporting of only the strongest climatic metric.

## Analysis scope

The analysis used the frozen 34-species baseline only (25 families; 20 within-population and 14 among-population cases) at the minimum 20 occupied-climate-cell threshold. The response was `among` (1 = geographically structured flower-colour variation; 0 = within-population polymorphism). For each climatic metric, the model was

`among ~ metric_z + effort_z`

where `metric_z` is the standardized focal metric and `effort_z` is standardized `log1p(n_climate_cells)`. Five metrics were analysed symmetrically: temperature breadth, moisture breadth, climatic heterogeneity, PCA dispersion, and PCA hull area. Family-clustered sandwich uncertainty, 9,999 label permutations, leave-one-family-out refits, collinearity diagnostics, OpenTree + Grafen `phyloglm`, and V.PhyloMaker2 time-scaled `phyloglm` sensitivities were evaluated. Holm-adjusted p-values across the five climatic metrics are reported as multiplicity context rather than as a replacement for effect-size interpretation.

Source workflow run: `31142541223`; artifact: `34species-environmental-niche-comprehensive` (artifact ID `8980386463`); head commit: `b82e7cd71db77b6f184aaa0bd0847d13c14858bd`.

## Family-clustered GLM, permutation, and family deletion

| Metric | OR | 95% CI | Clustered Wald p | Permutation p | Holm Wald p | Holm permutation p | LOFO OR range |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 0.817 | 0.384–1.739 | 0.6000 | 0.6029 | 1.0000 | 1.0000 | 0.636–0.973 |
| Moisture breadth | 0.412 | 0.180–0.947 | 0.0368 | 0.0475 | 0.1840 | 0.2375 | 0.306–0.465 |
| Climatic heterogeneity | 0.681 | 0.294–1.577 | 0.3700 | 0.3574 | 1.0000 | 1.0000 | 0.492–0.782 |
| PCA dispersion | 0.712 | 0.306–1.660 | 0.4317 | 0.3820 | 1.0000 | 1.0000 | 0.525–0.813 |
| PCA hull area | 0.577 | 0.312–1.067 | 0.0797 | 0.2382 | 0.3189 | 0.9528 | 0.489–0.671 |

All five point estimates were below one, and every leave-one-family-out estimate remained below one. Moisture breadth had the largest negative effect estimate and was the only metric whose unadjusted family-clustered confidence interval excluded one. Its unadjusted clustered and permutation p-values were 0.0368 and 0.0475, respectively; neither remained below 0.05 after Holm adjustment across the five metrics.

## Collinearity diagnostics

Because the five climatic metrics were fitted in separate models, collinearity among the climatic metrics themselves does not enter any single model. Diagnostics therefore concern the focal standardized metric and standardized occurrence-effort covariate in each model.

| Metric | Predictor correlation | Max VIF | Condition number |
|---|---:|---:|---:|
| Temperature breadth | 0.329 | 1.121 | 1.407 |
| Moisture breadth | 0.298 | 1.097 | 1.360 |
| Climatic heterogeneity | 0.412 | 1.204 | 1.549 |
| PCA dispersion | 0.295 | 1.095 | 1.355 |
| PCA hull area | 0.541 | 1.415 | 1.833 |

These diagnostics do not indicate problematic predictor collinearity.

## OpenTree + Grafen phylogenetic logistic sensitivity

OpenTree matching retained 30 species. One hundred polytomy-resolution fits completed for every metric. The point estimate was negative in every completed fit.

| Metric | Species | Median OR | Median 95% CI | Median p | Holm median p | Negative fraction |
|---|---:|---:|---:|---:|---:|---:|
| Temperature breadth | 30 | 0.908 | 0.436–1.893 | 0.7974 | 1.0000 | 1.00 |
| Moisture breadth | 30 | 0.573 | 0.234–1.403 | 0.2227 | 1.0000 | 1.00 |
| Climatic heterogeneity | 30 | 0.726 | 0.322–1.636 | 0.4398 | 1.0000 | 1.00 |
| PCA dispersion | 30 | 0.819 | 0.387–1.731 | 0.6006 | 1.0000 | 1.00 |
| PCA hull area | 30 | 0.654 | 0.261–1.641 | 0.3653 | 1.0000 | 1.00 |

## Time-scaled V.PhyloMaker2 sensitivity

All 34 species were retained under placement scenarios S1–S3. Point estimates remained below one for all five metrics in all scenarios. S1 and S3 were numerically identical; S2 differed only slightly.

| Metric | S1 OR (95% CI; p) | S2 OR (95% CI; p) | S3 OR (95% CI; p) |
|---|---|---|---|
| Temperature breadth | 0.827 (0.385–1.777; 0.6258) | 0.838 (0.393–1.789; 0.6484) | 0.827 (0.385–1.777; 0.6258) |
| Moisture breadth | 0.448 (0.167–1.206; 0.1121) | 0.454 (0.171–1.211; 0.1147) | 0.448 (0.167–1.206; 0.1121) |
| Climatic heterogeneity | 0.670 (0.290–1.550; 0.3491) | 0.679 (0.295–1.565; 0.3635) | 0.670 (0.290–1.550; 0.3491) |
| PCA dispersion | 0.698 (0.319–1.527; 0.3681) | 0.707 (0.324–1.542; 0.3838) | 0.698 (0.319–1.527; 0.3681) |
| PCA hull area | 0.599 (0.232–1.545; 0.2892) | 0.612 (0.242–1.549; 0.3002) | 0.599 (0.232–1.545; 0.2892) |

For moisture breadth, Holm-adjusted p-values across the five metrics were 0.5603 for S1/S3 and 0.5736 for S2. The phylogenetic analyses therefore preserve the effect direction and a comparatively large moisture-breadth effect size, while widening uncertainty so that all confidence intervals include one.

## Interpretation for the manuscript

The most defensible biological summary is not that moisture breadth is a uniquely established driver. Rather, geographically structured flower-colour variation is consistently associated with narrower sampled occupied climatic niches across all five metrics, with the strongest contrast for moisture breadth. This pattern is stable to leave-one-family-out deletion and remains negative under both phylogenetic treatments. However, only the unadjusted non-phylogenetic moisture model excludes the null, and the result does not retain conventional significance after correction across the five exploratory metrics or after phylogenetic correction. Occupied climatic breadth is not physiological tolerance, the analysis is species-level rather than morph-level, and the association is not causal.
