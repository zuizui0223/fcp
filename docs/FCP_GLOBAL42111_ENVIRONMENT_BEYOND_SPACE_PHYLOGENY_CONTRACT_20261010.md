# FCP 42,111 original flower colours — Does environment predict BEYOND spatial and phylogenetic structure? (2026-10-10)

## Estimand: trait geography, not an ordinary SDM

A normal species-distribution model estimates where species are present or likely to be observed. This study instead asks whether **the original photographed four-state FLOWER-COLOUR outcome is conditionally associated with environmental covariates, over and above source geography, spatial autocorrelation, species identity, and phylogenetic ancestry**. The worldwide 42,111-taxon one-photo breadth, the 85,337 repeated taxon×cell photographs, and the directly matched LCVP reference tree are distinct supports. They must NEVER be pooled as 127,448 independent taxa or treated as genomic colour states.

**Main incremental heldout estimand** within a FIXED evidence-defined cohort:

`Delta_env = L(M_geography + M_spatial + M_phylogeny) - L(M_geography + M_spatial + M_phylogeny + M_env)`

where lower heldout log loss/Brier is better. The models must be trained on precisely the same original photo labels, species, environmental completeness, folds and spatial/phylogenetic support. A positive increment is **predictive environmental information conditional on modelled space and tree**, NOT heritable adaptive selection. An observed positive score is not independent of sample selection, data-driven feature discovery, or residual species/geography confounding.

### Three orthogonal inference scales

1. **Interspecific atlas / 42,111 species**: 18,457 classifiable original photos, one original photo per taxon. Geography (true photo coordinates and elevation) + spatial smooth or distance kernel; taxonomic genus controls and region blocking assess transport. **A genus fixed effect is NOT a true phylogeny**; LCVP direct tips exist only for a selected 342/872 and 649/1,761 local original taxa. Thus **NO full 42,111 species phylogenetically controlled conclusion is identifiable from the currently available source tree**.
2. **True-phylogeny limited source subset**: the original LCVP direct species tips, actual dated branch covariance from the phylogenetic tree, spatial location/distance covariance and site-elevation controls must jointly enter the baseline. Incremental environmental covariate importance is scored after BOTH spatial and dated phylogenetic baselines. Compare heldout original genera and heldout source geographic cells separately and attempt crossed genus+cell holdout **only if supported**, never claim cross-phylogenetic/geographic robustness from same-group training. Report original direct tip-coverage selection biases and the full-cohort HOLD. Patristic distance as ONE covariate of source dyadic mismatch is an exploratory sensitivity, **not equivalent to fitting phylogenetic residual covariance**.
3. **Within-species across original 85,337 species×cell source photos**: species training-only colour baseline/fixed intercept, spatial source-cell holdout, expanded original site environment. An environmental increment here can test cross-region PHOTO contrast beyond species-invariant differences including their shared phylogenetic heritage, but the species-invariant phylogenetic main effect is mathematically **absorbed by the species intercept** and cannot simultaneously be independently estimated. Species represented in only one original classified region cannot identify their own environmental colour slope and must be excluded from the estimand while preserved in the source denominator. Within-species photo difference != genetic polymorphism.

### Source feature inventory under active extraction

- Original WorldClim 2.1 mean climate **BIO1–BIO19**: temperature means, warmest/coldest, diurnal/seasonality and precipitation volume and seasonality.
- Source-photo 12-month long-term **solar incident radiation, wind speed and water-vapour pressure**: their annual mean and monthly coefficient of variation. Radiation is not actual under-canopy light or the photo date.
- Original 0–30cm modelled **SoilGrids** pH, carbon, nitrogen, clay, water-retention proxy; additional CEC, sand, silt, bulk density and coarse fragments; retain providers' raw scaling and source layer SHA256.
- True observed photographed-source **elevation and geodesic location**; no equal-area-cell midpoint soil/weather sampling.

Seven locked ecological blocks: `elevation`, `temperature`, `precipitation`, `solar_radiation`, `wind`, `vapor_pressure`, `soil`. Report the independent conditional heldout improvement when omitting any block and, in an explicitly exploratory separate table, EVERY named individual variable. **Do not announce the best positive variable alone**. Co-linearity means drop-one-block improvements need not sum, and marginal temperature/soil effects can invert when rainfall/radiation are included.

### Nonnegotiable spatial and phylogenetic controls

- A random row split is INVALID. Different original photo sites can share genus, species, neighbouring climate cells, recorder opportunity and the same WorldClim raster pixel. Train-only normalization and spatial/phylogenetic baselines, genus/region-heldout and joint blocking where feasible are essential.
- Geodesic latitude/longitude and flexible **spatial covariance** must be treated separately. Adding latitude/longitude once is insufficient to claim spatial autocorrelation removed.
- Nominal genus-mean controls are not equivalent to true dated species-tree covariance. The currently frozen direct-tip backbone covers **342/872 (39.2%)** and **649/1761 (36.9%)** of preselected local original taxa and is environmentally selected. This fails the established >=50% whole-source tree gate even if exploratory subset model scores improve.
- Validate residual source PHOTO-COLOUR spatial structure after the environmental model (e.g., residual Moran's I by original geographic locality/micro-neighbourhood) and compare a **local exchangeability null** constrained within 50/100km site neighbourhoods or same-region taxa. Broad within-genus×20-degree-cell photo-label shuffles previously produced p=.010/.015 but geographically constrained shuffles gave p=.175–.880, and colour-blind geographically exchangeable source subset scores also failed. These results are NOT cured by adding seven blocks and selecting a promising coefficient.
- Outcome-exposed exploratory multiple-feature scans require explicit multiple testing control or later truly held-out source photophenotypes if claiming individual physical variables. A single nominal positive interval is NOT a validated mechanism; exact source result includes all null/negative cases.

## Decision statuses

`SUPPORTED_CONDITIONAL_PREDICTION` only if the same named block consistently reduces independent geographic, genus/phylogeny and, where estimable, crossed-block out-of-group loss, with practical sample sizes, robust local null calibration and acceptable residual spatial dependence. Otherwise `HOLD_OR_CONTEXT_DEPENDENT`. `NOT_IDENTIFIABLE` when no same-species replicated photo, direct phylogenetic branch covariance, or cross-group support exists. A photo source prediction remains observational even if passing this rubric; adaptive selection requires genetic colours/pigment/fitness evidence.

The scientific claim should be **"Is the visible flower-colour biogeographic signal environmentally predictive after spatial structure and evolutionary relatedness are accounted for?"**, not simply "which variable predicts photographed flower colour?". No current result verifies the affirmative across all 42,111 taxa. Main, old 1,499 photo-intensive H1/H2, and untouched prospective 2,000+730 taxa remain separate.
