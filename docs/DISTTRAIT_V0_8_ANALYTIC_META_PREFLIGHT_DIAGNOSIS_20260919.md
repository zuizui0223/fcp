# disttrait v0.8 preflight — analytic random-effects calibration diagnosis

Date: 2026-09-19 JST

Status: diagnostic only; not a promoted benchmark result.

## Purpose

Evaluate whether a conventional species-slope DerSimonian-Laird random-effects
summary can serve as the flexible model-based comparator requested after v0.7.

The exploratory 24-cell benchmark reused the v0.7 continuous direction-
heterogeneity design and added:

- one signed OLS slope per species;
- estimated within-species slope variance;
- DerSimonian-Laird tau^2;
- random-effects mean;
- Cochran Q heterogeneity test;
- Bonferroni combination of mean and heterogeneity probabilities.

## Exploratory result

GitHub Actions:

- PR: #48
- exploratory run: 35385955045
- artifact: disttrait-random-slope-meta-exploratory

The method behaved as intended under true directional heterogeneity, but its
analytic calibration failed in null worlds when per-species sample size was
reduced by 50% MCAR missingness.

Across the six null cells:

- matched distance-dissimilarity maximum FPR: 0.05;
- common-slope maximum FPR: 0.10;
- analytic random-effects mean maximum FPR: **0.175**;
- analytic Q heterogeneity maximum FPR: **0.30**;
- analytic combined omnibus maximum FPR: **0.30**.

The worst inflation occurred under 50% missingness.

## Diagnosis

The problem is not evidence that species-specific slopes are unusable.

The exploratory implementation treated estimated within-species slope variances
as if their large-sample meta-analytic calibration were adequate with roughly
8–20 observations per species. Under missingness, sampling variability in the
variance estimates propagates into the random-effects mean and Q statistic,
making the chi-square/normal approximations anti-conservative in this design.

## Decision

The analytic random-effects p-values are **not promoted** as a validated v0.8
comparator.

The species-slope estimands are retained, but final inference is being
recalibrated with matched within-species trait permutations:

1. preserve each species' predictor/coordinate values;
2. preserve each species' complete observed trait-value multiset;
3. permute trait assignments only within species;
4. refit every species slope and variance in every null world;
5. recompute the random-effects mean statistic and Q statistic;
6. calibrate both against their empirical matched nulls;
7. combine the two calibrated probabilities conservatively.

This preserves the model-based signed-slope interpretation while conditioning
on the actual species-specific sample sizes and trait distributions.

## Claim boundary

The failed analytic preflight must not be described as a validated hierarchical
comparator. If the matched-permutation calibration succeeds, the methodological
lesson is stronger:

> flexible species-specific slope summaries can represent directional
> heterogeneity, but small-sample meta-analytic test calibration should not be
> assumed; matched species-conditioned randomization can calibrate the final
> slope/meta statistics to the observed sampling geometry.
