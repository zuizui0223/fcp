# Dynamic TerraClimate protocol for the FCP space-time test

Status: **analysis design frozen before dynamic TerraClimate outcomes are inspected (2026-08-27)**. This is not an externally preregistered study. It is a prospective extension of the literature-driven v1 space-time hypothesis screen.

## Biological question

The core comparison is not whether a species has a large climatic niche. It is whether environmental variation is allocated differently through **space and time** in species whose flower-colour variation is documented as local coexistence (C), spatial segregation (S), or both (C+S).

Working predictions:

- stronger seasonal and interannual environmental movement/turnover -> C;
- stronger persistent among-location environmental variance -> S;
- a higher spatial:temporal environmental-variance ratio -> S;
- hypervolume size alone is not treated as the mechanism.

## Why occurrence year is not used

The existing GBIF file was built for spatial climatic sampling. It pages the occurrence-search API to a per-species cap and therefore does not constitute an unbiased temporal sample. A year audit found a large concentration of recent records and missing historical coverage for some informative species. Consequently, occurrence year is not used to define the climate time series.

Instead, the analysis fixes spatial locations first and follows **the same locations through climate time**.

## Fixed spatial support

Inputs are the QC-filtered GBIF occurrence sample from the fixed full-climate rebuild and the display-core-v6 focal-consistent species universe.

For each of the 66 climate-eligible core species:

1. calculate spherical distance from the occurrence-sample spherical centroid;
2. remove the outer 5% of distances as a geographic-outlier guard;
3. map retained records to the 1/24-degree TerraClimate grid and deduplicate cells;
4. select exactly 20 cells using deterministic spherical maximin sampling, starting from the retained cell nearest the spherical centroid.

All 66 species retain at least 20 unique TerraClimate cells after the radial guard. The resulting 1,320-point key is frozen by SHA-256 in `docs/supporting/jbi_dynamic_terraclimate_v1_plan_manifest.json`.

This does not turn the capped GBIF search sample into a complete species-range census. Inference is explicitly about the sampled occupied distribution.

## Dynamic climate source

Use only the current TerraClimate v1.1/current aggregated service for 1958-2025. Monthly variables requested are:

- `tmin` and `tmax`, combined to `tmean`;
- `ppt`;
- climate water deficit `def`;
- vapor-pressure deficit `vpd`.

The 1991-2020 period is the scaling/climatology baseline. The analysis does not mix legacy TerraClimate v1.0 with v1.1 and does not interpret long-term trends as independent evidence beyond the parent reanalysis. TerraClimate is used here for **environmental geometry and variability**, not trend attribution.

## Space-time geometry

Monthly climate is reduced to four dimensions: `tmean`, `ppt`, `def`, and `vpd`.

### Seasonal movement

For each species, calculate the 1991-2020 monthly environmental centroids at the 20 fixed points in globally standardized climate coordinates. The mean cyclic distance between adjacent calendar-month centroids is `seasonal_centroid_cycle_mean`.

Prediction: larger values favour C over persistent S.

### Annual hypervolumes

For each fixed point and year, summarize annual mean temperature, annual precipitation, annual climate-water deficit, and annual mean VPD. Standardize the four dimensions using the 1991-2020 baseline across all fixed point-years.

For each species-year, estimate a 95% shrinkage-Gaussian environmental ellipsoid using Ledoit-Wolf covariance. This avoids treating a high-dimensional KDE volume estimated from only 20 points as exact.

Prospective temporal metrics are:

- `interannual_centroid_drift_mean`;
- `interannual_overlap_loss_mean`, using Bhattacharyya affinity between consecutive annual Gaussian hypervolumes;
- `annual_hypervolume_log_sd`.

Prediction: greater temporal turnover/volatility favours C.

### Direct variance decomposition

In the same standardized annual climate coordinate system, decompose variability into:

- `spatial_variance_component`: variance among the long-term means of the 20 fixed locations;
- `temporal_variance_component`: mean within-location interannual variance;
- `space_time_variance_ratio`: spatial / temporal component.

Predictions:

- larger temporal component -> C;
- larger spatial component and spatial:temporal ratio -> S.

This variance decomposition is the primary conceptual bridge between the biogeographic C/S response and the space-time framing; hypervolume geometry is complementary.

## Statistical design

The conditional organization model is kept aligned with the rebuilt primary analysis:

- unresolved is never recoded as biological absence;
- documentation propensity is estimated from FCP-source count and outcome-independent literature attention;
- documented C-only, S-only and C+S species are compared with stabilized-IPW L2 multinomial models;
- each environmental metric is modeled with geographic extent as a covariate;
- uncertainty is evaluated with 499 family bootstraps;
- pure C vs S family-clustered binomial models are retained as a sensitivity and Holm-adjusted within the seven dynamic metrics.

The seed is `20260827`.

## Prospective decision rule

Seven dynamic metrics are fixed before observing the TerraClimate result. Directional agreement with the prediction plus a family-bootstrap interval that excludes 1 is treated as a candidate signal, not automatic causal proof. A signal that disappears under basic geographic or climate-position sensitivities will be labelled confounded/fragile rather than elevated.

If all seven dynamic metrics remain unsupported, the conclusion will be that broad climatic space-time geometry, at the present sampling scale, does not explain C/S organization. The next mechanistic layer would then move to pollinator turnover, breeding system, dispersal/gene flow and interaction-community structure rather than searching additional climatic summaries.
