# FCP edaphic + geographic-gradient decomposition protocol — 2026-10-07

## Motivation

The previous FCP IBD/IBE analysis used great-circle geographic distance and then a macroclimate PC distance. Two biologically important structures remain unresolved:

1. **edaphic differentiation** — soil pH, fertility, texture, bulk density, coarse fragments and water-holding properties can alter pigment physiology, plant stress and local fitness;
2. **directional geographic/topographic gradients** — latitude and elevation can structure flower colour beyond isotropic geographic distance.

This protocol therefore decomposes flower-colour turnover into five simultaneous components: isotropic distance, latitude, elevation, climate and soil.

## Frozen input layers

### Climate

WorldClim 2.1:
- BIO1–BIO19;
- mean monthly solar radiation from 12 SRAD layers.

### Elevation

WorldClim 2.1 10-arcminute elevation.

Elevation is **not** included in the environment PC. It is retained as its own topographic gradient.

### Soil

SoilGrids aggregated 5-km global mean products. The primary soil representation is a 0–30 cm thickness-weighted mean from the standard depth intervals:

- 0–5 cm, weight 5/30;
- 5–15 cm, weight 10/30;
- 15–30 cm, weight 15/30.

Soil features are:

1. pH in H2O (phh2o);
2. soil organic carbon (soc);
3. total nitrogen (nitrogen);
4. cation-exchange capacity (cec);
5. bulk density (bdod);
6. coarse fragments (cfvo);
7. clay content (clay);
8. sand content (sand);
9. available-water proxy = wv0033 - wv1500.

Silt is omitted because sand + silt + clay form a closed texture composition and silt is therefore redundant once sand and clay are included.

Native SoilGrids numeric scaling is retained because every feature is standardized on discovery before PCA; constant published scale factors therefore cannot affect standardized distances.

Exact source-file SHA256 values are frozen by the outcome-blind fcp-edaphic-elevation-checksum-freeze-20261007 workflow before the flower-colour analysis is run.

## Outcome-independent dimensionality reduction

Three environmental representations are fitted on **discovery observations only** without using flower-colour outcomes:

1. climate PCA;
2. soil PCA;
3. combined abiotic PCA = climate variables + soil variables.

For each representation:
- discovery mean and population SD are frozen;
- PCA is fit by deterministic SVD;
- the primary dimension is the smallest K explaining >=95% cumulative variance;
- the same centering, scaling and loadings are transported unchanged to validation and third cohort.

Fixed sensitivities:
- >=80% variance;
- >=90% variance;
- full standardized feature-space Euclidean distance.

## Pairwise predictors within species

For every unordered photo pair:

- Y = continuous nine-colour Jensen–Shannon flower-colour dissimilarity;
- G = great-circle geographic distance;
- LAT = absolute latitude difference;
- ELEV = absolute elevation difference;
- CLIM = Euclidean distance in the frozen discovery climate-PC space;
- SOIL = Euclidean distance in the frozen discovery soil-PC space;
- ABIOTIC = Euclidean distance in the frozen discovery combined climate+soil PC space.

All distances are unsigned because the response is colour dissimilarity rather than a signed colour axis.

## Main simultaneous decomposition

The primary four quantities use rank residualization with all listed covariates included symmetrically.

### Combined abiotic environment beyond geography/topography

rho_ABIOTIC = partial Spearman(Y, ABIOTIC | G, LAT, ELEV)

### Isotropic geographic distance beyond environment and directional gradients

rho_GEO = partial Spearman(Y, G | ABIOTIC, LAT, ELEV)

### Latitudinal gradient beyond distance, elevation and abiotic environment

rho_LAT = partial Spearman(Y, LAT | G, ELEV, ABIOTIC)

### Elevational gradient beyond distance, latitude and abiotic environment

rho_ELEV = partial Spearman(Y, ELEV | G, LAT, ABIOTIC)

These quantities directly test whether a directional north-south or elevational gradient adds information beyond simple pairwise distance.

## Climate-versus-soil decomposition

To ask whether soil adds information beyond climate, and vice versa:

rho_SOIL_unique = partial Spearman(Y, SOIL | G, LAT, ELEV, CLIM)

rho_CLIM_unique = partial Spearman(Y, CLIM | G, LAT, ELEV, SOIL)

This is the key edaphic test.

## Matched null

Within every species, complete continuous nine-colour vectors are permuted among fixed observation vertices 199 times.

The null preserves:
- exact coordinates;
- latitude and elevation;
- climate;
- soil;
- all covariance among geographic, topographic and environmental predictors;
- sample size;
- exact continuous colour-vector multiset.

Every partial-rank statistic is recalculated in each null world.

Species contribute equally to cohort means.

A component is supported only when:
- observed equal-species mean > 0;
- matched-null upper-tail p < 0.05.

## Replication

Discovery and validation must independently support a component before it is called replicated.

The third cohort is an additional species-disjoint transport and cannot rescue discovery-validation failure.

## Comparisons of interest

The analysis reports, without outcome-dependent switching:

- ABIOTIC versus GEO;
- LAT versus GEO;
- ELEV versus GEO;
- unique SOIL versus unique CLIM;
- combined ABIOTIC versus the earlier macroclimate-only PC result;
- macroclimate-PC versus earlier BIO5-only result.

## Interpretation ladder

1. **SOIL unique unsupported**: adding edaphic heterogeneity does not improve the environmental-sorting account beyond climate/geography.
2. **SOIL unique replicated, ABIOTIC still < GEO**: soil contributes nonredundant flower-colour turnover, but spatial/history structure remains stronger.
3. **ABIOTIC replicated and >= GEO**: measured climate+soil differences explain at least as much unique flower-colour turnover as isotropic geographic distance.
4. **LAT or ELEV replicated**: flower-colour turnover contains a directional geographic/topographic gradient not reducible to pairwise distance and measured abiotic environment.
5. **LAT/ELEV disappear after ABIOTIC adjustment**: their apparent effects are largely mediated by measured climate/soil gradients.

## Hard nonclaims

Even a strong climate+soil result does not establish:
- local adaptation;
- fitness advantage of local colour morphs;
- genetic differentiation rather than plasticity;
- causal soil or climate variables;
- pollinator-mediated selection;
- genetic isolation by environment.

The analysis measures **phenotypic environmental sorting beyond spatial structure**, not adaptation itself.
