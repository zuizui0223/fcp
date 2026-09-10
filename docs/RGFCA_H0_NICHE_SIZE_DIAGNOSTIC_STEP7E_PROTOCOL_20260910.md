# RGFCA Step 7E — H0 niche-size diagnostic

Date: 2026-09-10 JST
Status: prospectively frozen after Step 7D and before opening Step 7E diagnostics.

## Aim

Diagnose the unexpected two-tranche Step-7D signal in the negative-control total climatic hypervolume. This is a post-result diagnostic, not a rescue or a new confirmatory hypothesis.

## Fixed inputs

- Step 7D `species_environment_metrics.csv` and tranche-specific WorldClim cell tables.
- Frozen Step-7A photo organization states.
- Geographic extent is computed from the same colour-blind raw candidate coordinates used to construct Step-7D environment predictors.

## Geographic extent

Reproduce the preserved Chapter-1 definition from `scripts/data/compute_geographic_extent_metrics.py`:

1. spherical centroid of all valid raw candidate coordinates for each species;
2. great-circle distance of every coordinate to that centroid;
3. primary extent covariate = 95th percentile radius in km;
4. model covariate = `log1p(geographic_radius_95_km)`.

No colour labels enter this calculation.

## Analysis population

Within discovery and reserve separately, retain Step-7D species with `metric_status=complete` and pure photo states only:

- C* only = `local_cooccurrence_only`;
- S* only = `spatial_segregation_only`.

Do not pool tranches.

## Diagnostics

### E1 — state/extent association

Statistic: median log1p(radius95) in S* minus C*.
Two-sided 20,000-label permutation p-value.

Purpose: ask whether the photo-derived S* state preferentially occurs in geographically broader sampled species.

### E2 — hypervolume/extent association

Use `y = log1p(hv3d_rarefied20)` and `x = log1p(radius95)`.
Report Pearson and Spearman correlations in the full pure-state sample of each tranche.

### E3 — extent-adjusted S* coefficient

Fit OLS `y ~ 1 + S_star + z(x)` within each tranche. Report the S* coefficient and HC3 95% CI.

Inference uses a Freedman-Lane-style residual permutation under the reduced model `y ~ 1 + z(x)`: permute reduced-model residuals 20,000 times with fixed seed, add them back to fitted values, refit the full model, and compare the absolute S* coefficient. Two-sided p-value.

### E4 — radius-stratified label permutation

Create four quantile strata from `x` within each tranche using rank-based quartiles, so ties cannot collapse bins. Shuffle C*/S* labels only within those fixed strata, preserving the observed S* count in each stratum. Statistic: S* minus C* median difference in `y`. 20,000 permutations, two-sided.

## Diagnostic interpretation

- `extent_explains_signal`: Step-7D H0 is significant but E3 and E4 are both non-significant in both tranches, with E1 and/or E2 showing strong extent coupling.
- `extent_insufficient`: E3 or E4 retains p<=0.05 in both tranches with the same positive S* direction.
- otherwise `mixed_or_unresolved`.

No Step-7D directional hypothesis is reopened regardless of Step-7E outcome.

## Boundary

A residual S* association after geographic adjustment is still observational and may arise from photo-state construction, sampling geometry, omitted environmental dimensions, or biological niche breadth. Step 7E cannot identify selection or causation.
