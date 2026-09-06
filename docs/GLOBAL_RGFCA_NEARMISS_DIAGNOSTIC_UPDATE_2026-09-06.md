# RGFCA near-miss diagnostic update — 2026-09-06

## Result first

The positive full-domain G1 concentration excess is geographically fragile in an output-domain deletion diagnostic. The full observed/mean-null ratio is approximately 1.219; deleting the 10 evaluable Afrotropic-assigned grid cells reduces it to 0.824261. These cells carry only 0.8311% of the original aggregate opportunity weight. This is evidence about leverage in the existing field, NOT discovery of an African colour-transition zone, not an independent replication, and not a new significance test.

The original G1 remains unsupported (p = 0.070). No G1/G2 decision, species-disjoint decision, colour outcome, seed, null index, threshold, or kernel was changed.

## Exact scope of the new diagnostic

Contract: `docs/supporting/global_rgfca_frozen_field_realm_jackknife_diagnostic_contract_v1.json`.

Runner: `scripts/analysis/run_global_rgfca_frozen_field_realm_jackknife_diagnostic.py`.

Tests: `tests/test_global_rgfca_frozen_field_realm_jackknife_diagnostic.py`.

Results: `data/derived/global_rgfca_frozen_field_realm_jackknife_diagnostic_v1.csv`.

The diagnostic was fixed before calculating the deletion statistics. It reuses the original observed aggregate field, opportunity weights, evaluable mask, and all 999 existing aggregate null fields. It samples the exact frozen RESOLVE 0.01-degree integer ecoregion raster at the centers of the original 36 x 18 equal-area grid and maps row IDs through the original GPKG REALM field. All eight source realms and UNASSIGNED were included without outcome-based selection.

For each deletion, observed and null fields lose the same output cells; the opportunity-weighted mean and variance are recalculated on the retained domain. No new permutation is generated and no deletion p-value is reported. The original full-domain observed statistic and all 999 null statistics were reconstructed within numerical verification tolerance before reporting deletions; the original p = 0.070 was also reproduced.

### All deletion outcomes

Ratios and percentages below are rounded. The baseline ratio in the accompanying CSV is reported to six decimals; the runner regenerates full numerical precision.

| Deleted output-domain group | Removed cells | Removed opportunity weight | Observed / mean-null variance | Direction of excess |
|---|---:|---:|---:|---|
| None: full domain | 0 | 0% | 1.219004 | Positive |
| Afrotropic | 10 | 0.8311% | 0.824261 | Negative |
| Antarctic | 0 | 0% | Not evaluated | No original evaluable cells |
| Australasia | 40 | 5.1912% | 1.214592 | Positive |
| Indomalaya | 9 | 1.4528% | 1.215562 | Positive |
| Nearctic | 65 | 50.0629% | 1.226713 | Positive |
| Neotropic | 25 | 2.9902% | 1.236324 | Positive |
| Oceania | 1 | 0.0123% | 1.219361 | Positive; very small deletion |
| Palearctic | 82 | 16.3228% | 1.227017 | Positive |
| UNASSIGNED | 135 | 23.1369% | 1.300915 | Positive |

An unassigned grid center is not recoded to a terrestrial realm. Unassigned cells may include ocean/coastal kernel support. There were 367 original evaluable cells.

### What this does not complete

This deletes output-domain cells, not source photos or edges. Kernels from a deleted realm can still contribute to retained cells. Coarse grid-center assignment also does not literally excise the ecological area of a realm. Accordingly, this is NOT the fully refitted leave-one-realm-out analysis, does NOT complete the original strong-stability gate, and does NOT complete leave-one-major-family-out analysis. Those remain separate unresolved tasks. The pre-existing running-consensus strong-stability failure also remains unchanged.

Five focused local tests passed: weighted-variance calculation; re-centering after deletion; no false pass for empty regions; rejection of invalid values/weights; and enforcement of all 999 frozen null fields.

## Upstream measurement and eligibility geography

A separate descriptive audit used only photo ID, observation ID, species ID, coordinates and the classifiability flag from the frozen 50,000-row measured table. It did not fit any colour-environment relation.

| Photo-presence footprint on the same 36 x 18 equal-area grid | Measured frame | Final eligible/classifiable frame |
|---|---:|---:|
| Species | 500 | 369 |
| Photos | 50,000 | 21,424 |
| Cells containing at least one retained photo | 280 | 248 |
| Cells containing at least five distinct species with retained photos | 217 | 170 |

Photo-presence cells are NOT the kernel-supported G1 evaluable cells. The loss from 217 to 170 is 47 cells (approximately 21.7%) of this descriptive multi-species photo footprint. This combines photo classifiability and species eligibility; it does not isolate a causal missingness effect and does not imply that the frozen conditional estimand is undefined. It quantifies an important cost hidden by merely calling the upstream design complete.

Measured-table SHA256: `6d04b9798d3f98a8f3174f6fa2869a0bd8d9b0db216dd45b71bda1cbb47d04ea`.

## Interpretation corrections for the simulation line

- The 200 outer maps reuse species and photos and obey a balanced schedule. They are not 200 independent biological replicates.
- The S1 simulation uses continuous latent-score differences followed by ranks, not the complete four-group JSD/classification pipeline. It is a geometry/rank detectability diagnostic, not automatically end-to-end measurement power.
- Its calibration null has no spatial signal. The f_shared = 0, amplitude > 0 scenarios instead contain independent species-specific spatial boundaries. Their exceedance rates are crucial for separating detectable spatial structure from specifically shared architecture.
- The 369- versus 500-species frame comparison changes species eligibility, photo availability and the sampling schedule together. It does not by itself isolate the effect of missing photographs.
- S4 zero-mean heteroskedastic noise can increase variance; it does not alone create a nonzero expected equal-species mean. Any claimed false-positive inflation must be tied to the actual calibration rule.
- An 80%-power effect of 0.015 from the previous interaction supplement is the first successful tested grid value, not a continuously estimated exact detection boundary. Power conditional on an assumed effect is not evidence that the observed effect is the true effect.

## Source artifacts and execution state

G1 source: Actions run `34002956347`, artifact `9980198765`, `global-rgfca-g1-single-runner-primary-v1`, ZIP SHA256 `65553c98b4a3dbe92bb379a779c181d02b8a03bdd094663e8af0b7af28741912ee`.

RESOLVE source: Actions run `34007596955`, artifact `9982073937`, ZIP SHA256 `727e6dc800f2704a6bb7a58d8574d2fe146c8e462333b20c6c435e47c6b137394`.

The separate design-power recovery (`34038925164`) and environmental heterogeneity inference (`34039210812`) were still computing at the checks made for this update. Neither was rerun or reclassified by this diagnostic. Their final biological or simulation results are not reported as completed here.

## Scientific consequence

The current near miss should not be explained solely as a weak signal evenly repeated across the world. The saved field has substantial geographic leverage concentrated in a small contribution to its opportunity weight. Region-specific structure, effective support, and the distinction between shared and idiosyncratic species patterns therefore need to remain explicit. Partial sharing remains a hypothesis to test, not a conclusion already established by this audit.
