# FCP 42,111 photographed flower species: all-environment evidence after geography, space, species and LCVP (2026-10-10)

## Current scientific conclusion — source-verified

**The original global flower-photo atlas contains modest environmental predictive information beyond coarse geography. Some environmental information also improves prediction of an already-photographed species' colour at a different regional source site beyond its species-specific training colour mean. However the seven full thermal, rainfall, sunlight, wind, vapor, elevation and soil blocks do NOT provide a stable positive out-of-group increment over the combined true-geodesic spatial and dated-LCVP phylogenetic covariance model on its directly tipped, selected original-source species subset.**

The results are not mutually inconsistent: *different denominators, outcome grain and covariance controls* are being tested. These are observed photo categories, not inherited pigment genotypes. This is an exploratory scale-dependent question, not a claim of adaptation.

## Source-backed data grains and committed numerical results

| Investigation / dependent variable | Input denominator | Tested complete-case group | Geography and lineage treatment | Main result |
|---|---:|---:|---|---|
| Global BETWEEN original nominal species; original one-photo 4-colour category | 42,111 original taxa / 18,457 classifiable | **14,118 original species** with full seven environmental blocks | source real photo geolocation, latitude/longitude; group heldout by original genus and independent original equal-area photo cell; **not true branch covariance** | ALL environmental blocks vs geography log-loss +0.007128 genus-heldout / +0.005663 cell-heldout; no *single complete predictor block* robustly unique in both holdouts |
| Same original species at different photographed geography cells; 4-colour source category | 85,337 original taxon×cell photos /39,075 classified | **20,546 evaluated source photos / 5,127 taxa** for climate; **11,136 evaluated original photos** with all soil | species intercept learned on OTHER geographic cells of the same species; heldout full original cells | all abiotic blocks beyond nominal species means+geographical terms +0.004615 climate pool; +0.007988 soil complete, conditional Brier |
| Same-species contrast with additional nonlinear spatial basis | same 85,337 source | matching climate/soil study pools | original species intercept + nonlinear photo-site coordinate trends; original geography cell heldout | all environmental blocks still add conditional Brier +0.002950 climate, +0.005662 soil; residual spatial GP/field NOT modeled |
| BETWEEN different original species: actual dated LCVP tree+photo-distance spatial covariance | fixed photographed 872/1761 parent local congeners | 342 / 649 direct dated LCVP tips, 261 / 518 after all soil fields | **actual Brownian shared ancestral branch covariance**, original photo great-circle exponential spatial kernels 50km and 250km, genus and original photo cell heldout; 5-fold | genuine dated phylo branch covariance improves prediction in many source heldouts; adding all environmental predictors does NOT consistently improve further |

**Source files, successfully verified and committed to exploratory PR #144:**
- [Official WorldClim and SoilGrids original-source extraction receipt](../results/fcp_global42111_extended_climate_soil_photo_20261010/result.json)
- [Global 14,118 complete-case species, 36 individual variables and 7 blocks](../results/fcp_global42111_complete_multienvironment_photo_models_20261010/result.json)
- [Original 85,337 taxon-cell same-species fixed intercept contrast](../results/fcp_global42111_same_species_all_abiotic_20261010/result.json)
- [Same-species nonlinear spatial-basis contrast](../results/fcp_global42111_species_nonlinear_spatial_all_abiotic_20261010/result.json)
- [Seven blocks under TRUE spatial+LCVP covariance](../results/fcp_global42111_all_abiotic_beyond_phylo_spatial_20261010/result.json)

All source data retain exact previously photographed taxon, observation and photo IDs. Original source locations recovered for 100,348 of 100,543 distinct photo IDs and the full expanded source has 100,255 original site photos with solar, wind, vapor pressure. Soil extension has 74,957 original source-site photos with new soil layers. All unclassifiable original photographed flower colours remain missing; no original 1,499 high-depth or reserved independently future-selected 2,000+730 taxa are opened.

## Precise seven-block results for direct phylo+spatial, 500km parent source

Brier improvement: baseline out-of-group Brier error MINUS enlarged-model out-of-group error; **positive** helps predict. All seven blocks enter the enlarged model, so their *block-specific uniqueness* measures the fall in predictive quality when just that block is removed. None of the seven showed nominal positive 95% conditional intervals in **all four** genotype-free geographic/taxonomic heldout conditions (geodesic spatial kernels 50/250km, genus/region holdout). Multiple comparisons are uncorrected and the spatial kernel bandwidth/RIDGE are fixed from the exploratory protocol.

| Source group and heldout | All seven extra over geodesic space+actual dated LCVP | Independent precipitation group | Independent solar radiation | Independent soil (when eligible) |
|---|---:|---:|---:|---:|
| 649 directly tipped source species; 250km photo kernel; genus holdout | **−0.00495** | −0.00460 | +0.00040 | not conditioned on soil |
| 649; 250km; photo-cell holdout | +0.00200 | −0.00361 | +0.00154 | not conditioned on soil |
| 649; 50km; genus holdout | −0.00457 | −0.00373 | +0.00063 | not conditioned on soil |
| 649; 50km; photo-cell holdout | −0.00054 | −0.00313 | +0.00184 | not conditioned on soil |
| 518 fully soil-complete directly tipped source species; 250km; genus holdout | **−0.01331** | −0.00169 | +0.00091 | −0.00330 |
| 518; 250km; photo-cell holdout | **−0.00736** | −0.00235 | +0.00105 | −0.00463 |
| 518; 50km; genus holdout | **−0.01616** | −0.00113 | +0.00107 | −0.00441 |
| 518; 50km; photo-cell holdout | **−0.00821** | −0.00266 | +0.00108 | −0.00396 |

Some small positive solar or thermal contrasts appear in specific region-heldout folds, but none provides a robust **across-space-and-genus** full-block claim. Even the eight nested model/scoring conditions are not eight independent replications. No condition is an independent calibration of the original spatially matched flower-photo label-shuffle null that previously did not support a local precipitation inference.

## Inference classification and remaining scientific work

**SUPPORTED (observational):** Source photographed four-colour composition varies geographically and has some conditional environmental predictability, particularly among different sampled taxa when only coarse geography is controlled, and within some species when photographs from other regions form the training colour baseline.

**NOT ROBUSTLY SUPPORTED:** In the selected exact direct-tree-tip sources, the expanded abiotic blocks yield a consistent positive heldout gain **beyond both geodesic spatial covariance and actual dated Brownian phylogenetic covariance**. The full 872/1761-local-taxon phylogenetic study retains its precommitted **50% source direct-tip coverage HOLD** (actual 342/872=39.2%; 649/1761=36.9%), and photo-site climate differs between tipped and untipped taxa.

**NOT YET IDENTIFIED:** Gene-mediated flower pigment plasticity vs adaptation, source true species-modal flower colours, photograph ROI classifier misclassification, local population co-occurrence and pollen/seed fitness; residual fine spatial autocorrelation after within-species nonlinear spatial model, fully resolved intra-genus phylogenetic branch uncertainty, and source-selection adjustment for incomplete dated tips.

**Actionable next biological identification:** make residual spatial/observer/photo classification diagnostic a hard gate on the >20k within-species heldout photographs; a LOW-RANK SPATIAL RANDOM FIELD and independent local source photophenotype replications would give stronger spatial specificity. A polynomial latitude-longitude spatial basis is not a distance-based Gaussian field. A model with fixed species intercepts does not yield a separately identified species-invariant phylogenetic main effect. Do not call a positive Brier difference causal selection.

All results are *historical, outcome-exposed, multiple-comparison exploratory*. In no circumstance does this file alter source FCP 1499 deeply photographed H1/H2 or independently reserved future 2000+730 taxa.
