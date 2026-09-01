# FCP image-first atlas execution status

## Current decision

The active mainline is a global image-first flower-colour atlas with a species-free public map and photo bar, but species-conditioned measurement, spatial fields and null inference.

Current gate: **the 200-species live metadata audit, live-geometry environmental coverage and both independent ROI v4 gates have passed; the versioned many-to-many dated-source resolver and final dated-source environmental-coverage rerun remain; every atlas candidate pixel remains closed**.

The fail-closed state is unchanged: bulk atlas image opening is stopped until the exact dated-source, environmental-coverage and locked ROI v4 gates all pass.

The terminal scale-out is eight disjoint cohorts of 25 species x 300 observations: 200 species and 60,000 observations. All eight cohorts are required. They are one experiment, not eight opportunities to stop after a favourable result.

## Immutable earlier evidence

- Six-species Chapter 1: 1,200 photographs, Stage A `p = 0.0113`, Stage B `p = 0.0906`. No atlas change retunes it.
- Literature comparison: 34 species from 25 families. It remains method and supplementary evidence, never the atlas admission frame.
- Automated three-species validation: 717 locked photographs, 306 admitted encounters, zero supported species. This is a frozen negative measurement/spatial validation, not proof of spatial randomness.
- Fifty-species sentinel geometry: 20,200 observations frozen before pixels. Its images remain unopened and it is not silently promoted into the 200-species experiment.

## Frozen branch outcomes before atlas colour

### Geographic shared-boundary concentration — `not_evaluable`

The exact-geometry qualification used 100 repetitions per scenario and 999 within-species permutations. The proposed geographic statistic failed the predeclared signal-recovery gate (`effect = 2` power 0.21; heterogeneous-boundary false sharing 0.16). This is a method failure, not evidence that a common geographic boundary is absent. The statistic is prohibited from reading real scale-out colour.

### Environmental concordance — next evaluable branch

The environmental null recovered its intended calibration before colour. Three independent primary families are frozen on the same 100, 250 and 500 km equal-area grids:

1. WorldClim 2.1 macroclimate;
2. ESA WorldCover 2021 land-cover composition;
3. RESOLVE Ecoregions 2017.

Copernicus terrain access did not complete under the frozen source route and terrain is permanently `not_evaluable`; no replacement DEM may be chosen after colour. Final inference uses within-species continuous Lab transition fields, equal-species surfaces, equal-cohort aggregation and one 9,999-randomization joint maximum across all frozen scales, date sensitivities, flower/background roles and environmental overlays. Only the predeclared 100 km, all-dates, equal-eight-cohort flower tests can promote a primary claim.

### *Bombus* biogeography — `not_evaluable`

The frozen colour-blind source gate required a citable GBIF download with its own realm, opportunity, balance and stability checks. The required authenticated download could not be frozen, so the branch remains `not_evaluable`. It cannot be replaced by a convenient post-colour occurrence source.

## ROI evidence

- ROI v3 (SegFormer) failed its independent JRC development gate: 400 development images, 17 admitted, recall 0.1207. The locked 100-image test remains sealed forever.
- ROI v4 is a prospectively frozen YOLO11n flower detector followed by EfficientSAM using one shared runtime for JRC qualification and all 60,000 atlas measurements.
- The final training executable, materialization receipt, terminal epoch weight and training curves were committed before any development prediction.
- The 400-image development gate passed all frozen criteria: 351/400 images admitted (0.8775), detector precision 0.7760, recall 0.8465 and pooled flower-mask pixels inside the reference-box union 0.9064. The committed evidence passed CI before the locked partition was opened.
- The locked 100-image test then passed every frozen criterion without retuning: 85/100 images admitted (0.85), detector precision 0.7304, recall 0.7956, medium-object recall 0.8296, large-object recall 0.5000, median image mask containment 0.9139 and pooled mask containment 0.8597. The row-level result, executed gate code and exact trained weight are frozen under `data/atlas/qualification/roi_v4_locked_test/`.
- ROI v4 now authorizes the estimator for scale-out. It does not independently authorize atlas image acquisition; the dated-source and final environmental-coverage gates remain closed.

## 200-species metadata and dated-source gate

GitHub Actions run `33405153936` completed the full predeclared 500-species candidate-pool audit successfully. It did not stop after reaching the target. Of 500 queried candidates, 358 passed metadata and geometry admission. The frozen terminal selection contains 200 genus-distinct species in eight disjoint 25-species cohorts, with exactly 300 observations per species: 60,000 unique observations and photos. Retained maxima were 10 records per species-observer, 5 per species x 0.25-degree cell and 10 per species x 0.5-degree cell. All selected rows still declare that candidate pixels were unopened.

The same selected geometry passed the live-feasibility environmental opportunity gate at all three scales. Macroclimate and ecoregion coverage was 99.05%, 99.72% and 100% at 100, 250 and 500 km; land-cover coverage was 100%, 100% and 99.12%. This is pre-colour coverage evidence only. The live API result cannot authorize images, and the coverage gate must be repeated against the exact dated-source reconciliation.

The official fixed resolver is the 2026-08-27 iNaturalist Open Data snapshot:

- object size: 35,093,052,336 bytes;
- computed SHA-256: `c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d`;
- moving `latest` is prohibited.

The official snapshot has observations, observers, photos and taxa but no flowering-annotation table. The pre-image contract therefore has two source stages: one complete API audit selects the exact 60,000 rows, then every photo, observation, observer and taxon must reconcile to the dated snapshot.

The v1 resolver stopped `not_evaluable_dated_source_reconciliation` at selected photo asset `950871`, before any association row or image content was inspected, because it assumed that `photo_id` was a unique table key. iNaturalist documents a many-to-many photo-observation relation and the composite association key `photo_uuid + observation_uuid`. The v1 result is retained as a technical source-schema failure. A separately versioned v2 resolver was frozen before inspecting the association rows. It allows multiple observation links for a photo asset but requires exactly one full-metadata match to the already frozen API observation. Zero matches, multiple matches, duplicate composite association keys or conflicting photo-asset fields yield `not_evaluable_dated_source_m2m_reconciliation`; no replacement or resampling is allowed.

## Image-access authorization

The measurement firewall cannot be built unless all of these pass independently:

1. exact v2 many-to-many dated-source reconciliation for 200 species and 60,000 photo assets;
2. opportunity-cell environmental coverage at 100, 250 and 500 km, with macroclimate plus at least one other primary family;
3. ROI v4 locked JRC test.

The acquisition worker additionally checks the exact sealed-key hash and the firewall receipt. Passing metadata, a visually plausible map, or a development model alone cannot open images.

## Remaining route to a submission result

1. Complete the exact v2 many-to-many dated-source reconciliation and repeat pre-colour environmental coverage against that dated source.
2. If both source gates pass, acquire and measure all 60,000 images location-blind; incomplete measurement stops the coordinate join.
3. Run the frozen environmental joint null once, preserving `supported`, `not_supported` and `not_evaluable` outcomes.
4. Generate the species-free map/photo bar, species-conditioned tables, manuscript, reproduction bundle and CI-complete PR.

No colour-dependent retuning, early stopping, silent row repair or new environmental/pollinator overlay is permitted.
