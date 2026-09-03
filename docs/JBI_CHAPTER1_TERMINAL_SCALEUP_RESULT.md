# JBI Chapter 1 — terminal 200-species scale-out result

## Editorial role

This terminal scale-out is **not a second biological test of the Chapter-1 hypothesis**. It is a prospectively frozen extension designed to determine whether the image-first pipeline could yield an evaluable 200-species colour field before any protected coordinate-colour inference was opened.

For the Journal of Biogeography manuscript, the six-species held-out Stage A/B analysis remains the primary biological evidence. The 200-species experiment belongs in the reproducibility / limitations narrative and Supporting Information because its frozen measurement-completeness gate stopped the cascade before species-conditioned spatial organization.

## Frozen design

The terminal experiment contained eight disjoint cohorts of 25 genus-distinct species, with 300 fixed observations per species: **200 species and 60,000 observations/photos**. The 60,000-record denominator, source roles, estimator, thresholds, seeds, cohort requirements, branch order and `not_evaluable` stopping rule were fixed before terminal outcome inspection. Replacement, resampling, favourable partition selection and early stopping were prohibited.

The formal sequence was:

`measurement completeness -> species-conditioned spatial organization -> shared transition -> environmental concordance -> pollinator biogeographic concordance`

A `not_evaluable` gate never advanced to the next stage.

## Exact execution result

Exact GitHub Actions run: `33592829701`, attempt 2, frozen execution head `aea19a4eff9585f501aa6a833ad44c80080eddcf`.

Exact artifact: `jbi-atlas-terminal-measurement-v5`, artifact `9857419176`, SHA-256 `6bfff229a90215d016b2dd7e2bcca9446474f6d216efe33782f9faac11e69b55`.

The location-blind compute and reassembly layer completed exactly as frozen:

- 256 / 256 compute partitions;
- 16 semantic shards;
- 60,000 unique measurement IDs;
- 60,000 terminal measurement records;
- no coordinate opening during measurement;
- no persisted candidate image pixels in the reassembly evidence.

The exact reassembly status was `pass_exact_256_compute_partition_coverage`, and the exact measurement-bundle status was `pass_complete_location_blind_roi_v4_measurement_v5_bundle`.

The subsequent frozen measurement gate returned:

**`not_evaluable_scaleout_measurement_completeness`**.

Only **58 / 200 species** met the frozen measurement-evaluable rule. Cohort counts were:

| Cohort | Evaluable species |
|---|---:|
| C01 | 3 |
| C02 | 4 |
| C03 | 7 |
| C04 | 8 |
| C05 | 8 |
| C06 | 12 |
| C07 | 7 |
| C08 | 9 |

All eight cohorts were therefore `not_evaluable`. The binding gate state was `coordinate_join_permitted = false` and `coordinates_opened = false`.

## Scientific decision

The terminal confirmatory cascade stopped **before species-conditioned spatial organization**. No shared-transition, environmental-concordance or pollinator-concordance result exists for the 200-species experiment.

This outcome must not be translated into biological absence. In particular, it does **not** show that flower colour lacks spatial organization, that shared transitions are absent, or that environmental or pollinator concordance is unsupported. Those hypotheses were not evaluated because the frozen measurement-completeness requirement was not met.

The scale-out instead identifies a measurement-scaling limit of the current frozen ROI-v4 pipeline: increasing taxonomic breadth and image count did not by itself produce enough measurement-evaluable species in every fixed cohort to authorize coordinate-colour inference.

## Manuscript-ready Methods sentence

A prospectively frozen terminal scale-out independently selected eight disjoint 25-species cohorts with 300 observations per species (200 species; 60,000 photographs) and required complete location-blind ROI measurement before any coordinate-colour join. The predeclared cascade stopped as `not_evaluable` whenever the fixed species/cohort measurement-completeness rule failed; no downstream spatial, environmental or pollinator test could then be opened.

## Manuscript-ready Results paragraph

A prospective 200-species scale-out completed all 256 location-blind compute partitions and reassembled exactly 60,000 unique measurement records into 16 semantic shards, but it did not pass the frozen measurement-completeness gate. Only 58 of 200 species were measurement-evaluable, with 3, 4, 7, 8, 8, 12, 7 and 9 evaluable species across the eight fixed 25-species cohorts. The gate therefore returned `not_evaluable_scaleout_measurement_completeness`, coordinate joining remained prohibited and no species-conditioned spatial or downstream concordance analysis was run. This is a measurement-evaluability result rather than evidence for or against the biological spatial hypotheses.

## Manuscript-ready Discussion paragraph

The attempted terminal expansion also exposes an important scaling constraint. Although exact location-blind measurement was completed for all 60,000 frozen records, the current frozen ROI-v4 pipeline yielded enough evaluable measurements for only 58 of 200 species and for none of the eight cohorts under the predeclared completeness requirement. We therefore stopped before opening coordinate-colour inference rather than relaxing the gate after observing measurement performance. The six-species result should consequently be interpreted as evidence that the held-out pipeline can recover spatial organization in selected taxa, not as proof that the same measurement representation scales uniformly across angiosperms. A future global atlas will require a new prospectively qualified measurement system with broader taxonomic transfer, rather than post-result reanalysis of these 60,000 records under substituted thresholds.

## Claim boundary

Allowed:

> The prospectively frozen 200-species extension was not evaluable because the predeclared measurement-completeness gate was not met, despite exact completion of the 60,000-record location-blind measurement bundle.

Not allowed:

- the 200-species experiment found no spatial organization;
- the 200-species experiment refuted a shared boundary;
- environmental concordance was unsupported;
- pollinator concordance was unsupported;
- the 58 measurement-evaluable species may be analysed as a favourable post-result subset;
- thresholds or metrics may be substituted to rescue the experiment.

## Durable evidence

- `docs/supporting/jbi_atlas_terminal_measurement_v5_receipt.json`;
- `docs/JBI_IMAGE_FIRST_ATLAS_STATUS.md`;
- run `33592829701`;
- artifact `9857419176`.
