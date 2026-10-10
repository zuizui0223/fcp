# FCP original 85,337 photo-cells: exact-source geodesic spatial-kernel results and residual Moran audit

Source date 2026-10-10. CI [Actions 38014630364](https://github.com/zuizui0223/fcp/actions/runs/38014630364) succeeded after synthetic original-site and geographic-nearest-neighbour guards; [machine result](result.json) is source-committed.

## Test: can physical environment predict photographed floral colour beyond named species and photo-site spatial covariance?

Original 85,337 source species×equal-area-cell photographs; 39,075 already four-colour classifiable and 46,262 unclassifiable retained outside the prediction estimator. In the **climate-no-soil** common-case population 38,968 classified original photographed cells, 20,546 original heldout photos of 5,127 different named taxa can be predicted from their OTHER region's source photo labels. With ten SoilGrids variables complete, 27,003 source photo-cell cases give 11,136 heldout photo cells from 3,401 taxa. Original public locations and photo IDs were preserved exactly, no new photographs or colour outcomes.

Fixed fivefold **source geographic-cell heldout** comparisons in exactly the SAME original test photo IDs: species means estimated from other TRAIN source cells only; latitude/longitude controls; add original geodesic photo-site exponential spatial distance kernel estimated with 128 original TRAIN-only, outcome-blind Nyström landmark photos; then add all climate/elevation/solar/wind/vapour source physical blocks (and ten modeled soil traits in soil-complete cohort). Two spatial length scales **100km and 500km** specified and all results reported. These are low-rank spatial source covariate features in a species-centred ridge photo-colour predictor, NOT a fully estimated spatial random field.

| Support / photo geodesic-kernel length | Species+latitude/longitude original photo Brier | +actual site spatial feature basis | +all abiotic blocks | Environmental improvement beyond spatial basis | Original cell-block fixed-fit 95% CI |
|---|---:|---:|---:|---:|---|
| Climate only, 20,546 test photos, 100km | 0.402327 | 0.398295 | 0.394917 | **+0.003378** | [+0.002853,+0.003833] |
| Climate only, same 20,546, 500km | 0.402327 | 0.395443 | 0.393582 | **+0.001861** | [+0.001372,+0.002398] |
| Climate+soil, 11,136 test photos, 100km | 0.427283 | 0.421270 | 0.415888 | **+0.005381** | [+0.004594,+0.006199] |
| Climate+soil, same 11,136, 500km | 0.427283 | 0.416461 | 0.413315 | **+0.003146** | [+0.002382,+0.003984] |

All four source outcomes also have positive *species*-cluster conditional fixed prediction 95% intervals. Thus the within-species photographed between-region **environmental predictive increment survives this particular source geodesic low-rank covariance approximation**. This does NOT imply independent environmental adaptation, allele turnover or fully removed spatial structure.

On the climate-only 500km-kernel cohort, individual block contributions after all OTHER environmental groups: precipitation **+0.000641** (geo-cell 95% [+0.000285,+0.000947]); temperature **+0.000377** ([+0.000181,+0.000616]); wind **+0.000163** ([+0.000028,+0.000331]); solar +0.000118 (cell interval includes 0), humidity/vapor +0.000031 (interval includes 0), elevation +0.000090 (cell interval reaches slightly below zero). In the soil-complete 500km cohort, source soil ten-trait block adds +0.000655 ([+0.000300,+0.001044]); precipitation +0.000982 ([+0.000523,+0.001398]). All these are post-outcome conditional feature checks across correlated blocks, not independent validated mechanisms.

## Residual space audit added afterward: **do not claim all spatial confounding was eliminated**

The original full-cohort and fixed fold forecasts were reproduced *unchanged* in success [Actions 38014630364](https://github.com/zuizui0223/fcp/actions/runs/38014630364), then original heldout photo residuals for all four colour classes were compared among the closest 8 other ACTUAL photographed sites within 100km and 500km. The values below are **mean absolute descriptive four-class Moran I**, a directed near-neighbour statistic without null significance calibration. They are not standardised effect sizes to compare with heldout Brier gain:

| Source cohort, source distance-kernel length, residual neighbour cap | Species+basic geography | +source spatial covariance approximation | +all original abiotic |
|---|---:|---:|---:|
| Climate 20,546; model 500km; neighbour radius 100km | 0.00382 | 0.00465 | 0.00478 |
| Climate 20,546; model 500km; neighbour radius 500km | 0.00356 | 0.00415 | 0.00428 |
| Soil-complete 11,136; model 500km; neighbour radius 100km | 0.00849 | 0.01013 | 0.00961 |
| Soil-complete 11,136; model 500km; neighbour radius 500km | 0.00622 | 0.00754 | 0.00719 |

The residual diagnostic is SMALL but **not identically zero** and its mean absolute local colour residual correlation **did not uniformly decrease** after adding approximate spatial covariance or environmental predictors. There is no spatially constrained null/p-value here. Therefore **even this strengthened model must NOT be described as proving all source space/autocorrelation is removed.** Higher-rank random fields, explicit spatial null and verification of environmental novelty independent of photography/taxon sampling are still required for that stronger identification.

## Cross-scale joint decision

- **Within species in the source repeated species×region photo sample:** environmental data can improve original four-colour photo prediction after training-only species intercepts and an approximate real-geodesic spatial kernel across two prespecified spatial scales. This is a positive descriptive *conditional prediction*, not an environment-driven genetic phenotype mechanism.
- **Between different species with REAL dated phylogeny:** [649 (or 518 soil-complete) direct LCVP tree taxa](../fcp_global42111_full_abiotic_spatial_LCVP_20261010/result.json) do not show stable incremental all-environment benefit beyond Brownian shared-ancestry and geodesic spatial covariance. This is a distinct, environmentally biased subset and cannot replace a global 42,111-tip tree.
- **Observed location vs ancestry:** species intercept absorbs species-invariant phylogenetic effects, so within-species estimation does not independently estimate a phylogenetic covariance component. The two studies have different support, estimands and model families; apparent different results do not prove evolutionary scale shifts.
- **Strong causal claims remain HOLD:** photo colours are one original picture per taxon-region, not genetically observed alleles. Existing full 42,111-species phylogenetic coverage, ecological exposure, colour classification error, spatial matched null and fitness data are insufficient to infer local climatic flower-colour adaptation.

No original 1,499 detailed H1/H2 manuscript, main, or reserved independent 2,000+730 source taxa changed.
