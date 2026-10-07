# FCP full geographic–climate–soil gradient partition protocol — 2026-10-07

## Status

Post hoc spatial-process decomposition defined after the BIO5-only and macroclimate-only IBD/IBE-like analyses.

The purpose is not to rescue a preferred adaptive interpretation. It asks whether flower-colour turnover retains independent associations with **geographic separation, latitudinal position, elevation, macroclimate and soil** when these spatially correlated gradients are considered together.

No outcome changes frozen H1/H2/D-spatial decisions.

## Data sources

### Flower colour

Use the exact frozen discovery, validation and third-cohort measured tables used by the preceding spatial analyses.

Response:
- continuous nine-colour flower Jensen–Shannon dissimilarity.

### Macroclimate

Use the same fixed WorldClim 2.1 macroclimate block as the preceding macroclimate-PC analysis:
- BIO1–BIO19;
- mean monthly solar radiation.

Discovery-only standardisation and PCA are transported unchanged to validation and third.

### Soil

Use SoilGrids 2.0 **5-km aggregated mean layers** for a reproducible global macroecological soil representation.

SoilGrids source properties:
- pH in water (`phh2o`);
- soil organic carbon (`soc`);
- total nitrogen (`nitrogen`);
- cation exchange capacity (`cec`);
- bulk density (`bdod`);
- coarse fragments (`cfvo`);
- clay;
- sand;
- volumetric water content at 33 kPa (`wv0033`);
- volumetric water content at 1500 kPa (`wv1500`).

Silt is omitted because clay + sand + silt are compositionally closed; retaining all three would duplicate one degree of freedom.

For every property use the official SoilGrids mean maps for:
- 0–5 cm;
- 5–15 cm;
- 15–30 cm.

The 0–30-cm topsoil value is the thickness-weighted mean:

`(5*x_0_5 + 10*x_5_15 + 15*x_15_30) / 30`.

Mapped integer values are converted to the conventional units documented by SoilGrids before transformation. Positively skewed `soc`, `nitrogen`, `cec`, and `cfvo` are transformed with `log1p`; the remaining soil variables are left on their conventional scales.

Soil standardisation and PCA are fit on discovery observations only and transported unchanged to validation and third.

### Latitude and elevation

- latitude is the photo latitude already frozen in each observation table;
- elevation is WorldClim 2.1 10-minute elevation, the same global surface family underlying the WorldClim climate data.

Pairwise geographic/topographic gradients are:
- great-circle distance;
- absolute latitude difference;
- absolute elevation difference.

Longitude is not entered as a separate predictor because east–west separation is already contained in great-circle distance; latitude is retained explicitly to test a directional broad-scale biogeographic gradient requested a priori.

## PCA rules

### Climate PC

Identical to the frozen macroclimate-PC protocol.

### Soil PC

Fit deterministic SVD to discovery-only standardised soil variables.

Primary soil dimensionality:
- smallest K with >=95% cumulative soil variance.

Sensitivities:
- >=80%;
- >=90%;
- full standardised soil space.

### Combined abiotic PC

Construct a block-balanced climate+soil matrix from discovery standardised variables:

- each climate variable is divided by sqrt(number of climate variables);
- each soil variable is divided by sqrt(number of soil variables).

This gives the climate and soil blocks equal total variance before fitting the combined PCA, preventing the larger climate block from receiving greater weight solely because it contains more variables.

Fit the combined PCA on discovery only and transport its means/transform unchanged.

Primary combined dimensionality:
- smallest K with >=95% cumulative block-balanced variance.

## Model A — unique component decomposition

Within every eligible species and unordered photo pair calculate:

- Y = flower-colour JSD;
- G = great-circle distance;
- LAT = absolute latitude difference;
- ELEV = absolute elevation difference;
- CLIM = Euclidean distance in frozen climate PC95 space;
- SOIL = Euclidean distance in frozen soil PC95 space.

For each focal predictor X, calculate a multiple partial-rank correlation:

`rho_X = cor(residual(rank(Y) | ranks(other 4 predictors)), residual(rank(X) | ranks(other 4 predictors)))`.

Thus each component is the association remaining after all other measured gradients are removed.

The five focal components are:
1. geographic distance;
2. latitude;
3. elevation;
4. macroclimate;
5. soil.

## Model B — combined abiotic environment versus geography/topography

Replace CLIM and SOIL by one block-balanced combined ABIOTIC-PC distance and calculate unique partial-rank effects for:

1. geographic distance;
2. latitude;
3. elevation;
4. combined climate+soil abiotic distance.

This is the primary test of whether measured abiotic environment collectively explains flower-colour turnover beyond broad geographic/topographic structure.

## Matched null

Within each species, permute complete continuous flower-colour vectors among fixed observation vertices 199 times.

The null preserves:
- all coordinates;
- latitude and elevation;
- climate and soil values;
- covariance among all geographic and environmental gradients;
- exact continuous colour-vector multiset;
- sample size and species-wide colour distribution.

Each partial-rank component is recalculated in every null world.

## Cohort inference and multiplicity

Species contribute equally.

For Model A, within each cohort:
- report the equal-species mean partial rho for all five components;
- calculate matched-null upper-tail p-values;
- apply Holm correction across the five component tests.

A component is cohort-supported only when:
- mean rho > 0;
- Holm-adjusted p <= 0.05.

A component is replicated only when it is supported in both discovery and validation. Third cohort is a species-disjoint transport and cannot rescue discovery–validation failure.

For Model B, apply Holm correction across its four component tests under the same rule.

## Comparisons of interest

Predefined comparisons:

1. Does soil retain a positive unique component after distance, latitude, elevation and climate are controlled?
2. Does macroclimate retain a positive unique component after soil and geographic/topographic gradients are controlled?
3. Does elevation retain a component independent of macroclimate and soil?
4. Does latitude retain a component independent of great-circle distance and measured environment?
5. Does combined climate+soil abiotic distance retain a replicated unique component?
6. Is the combined abiotic component larger than the unique great-circle-distance component?

The last comparison uses the species-wise difference in partial rho and a deterministic 9,999-replicate sign-flip test plus the matched permutation distribution.

## Interpretation

Possible outcomes:

- **distance only**: spatial/history/dispersal structure dominates under measured environment;
- **soil or climate + distance**: environmental sorting contributes nonredundantly but does not replace IBD/history;
- **elevation after climate+soil**: a topographic gradient remains that may proxy unmeasured radiation, exposure, pollinator or microclimate processes;
- **latitude after distance+environment**: broad biogeographic structure remains beyond measured abiotic conditions;
- **combined abiotic > distance**: strongest observational evidence in this dataset that measured environmental differences organize colour turnover more strongly than pure separation.

## Hard nonclaims

Even a strong combined abiotic result is not:
- genetic isolation by environment;
- a fitness test;
- local adaptation;
- proof of causal soil or climate selection;
- proof of genetic differentiation rather than plasticity.

SoilGrids is a modelled global soil product and is not a local soil measurement at the photographed plant. WorldClim elevation/climate are likewise macroenvironmental surfaces.

Pollinator assemblages, herbivory, direct UV-B exposure, land use, demographic history and neutral population-genetic structure remain unmeasured.
