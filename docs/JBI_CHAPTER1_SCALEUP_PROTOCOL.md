# JBI Chapter 1 — 12-species scale-up precursor (superseded as the terminal design)

## Status after the image-first atlas redesign

This document preserves PR #21's first prospective scale-up design for auditability. Its 12 literature-selected species × 200 photographs are a validated metadata-method precursor, **not** the final FCP cohort or active mainline.

The active prospective protocol is [`JBI_IMAGE_FIRST_ATLAS_PROTOCOL.md`](JBI_IMAGE_FIRST_ATLAS_PROTOCOL.md). Its 50-species, 20,200-observation cohort is an unopened sentinel; the terminal experiment is the separately frozen v2 design with eight disjoint 25-species cohorts x 300 observations (200 species and 60,000 observations). The literature ledger and frozen 12-species outputs below remain unchanged and may be used to test acquisition machinery, but they cannot select or substitute atlas species.

## Current boundary

The completed six-species analysis is closed and remains unchanged:

- 1,200 photographs split into 480 calibration and 720 held-out evaluation records;
- Stage A supports within-species spatial organization of continuous flower colour (`p = 0.0113`);
- Stage B does not confirm one universal global concentration of transition boundaries at the frozen primary support (`p = 0.0906`).

The scale-up is a new prospective phase. It does not retune, replace or reinterpret the completed six-species evaluation.

## Purpose

The six-species study established that the pipeline can detect spatial organization while preserving species distributions in the null model. It did not provide enough evidence to conclude that independently reconstructed transitions repeatedly occupy the same global cells.

The next question is therefore:

> Does the distinction between repeated within-species organization and non-confirmed shared geography persist when the number of independently sampled colour-variable species is increased under a prospectively frozen, colour-blind admission rule?

## Scale-up cohort v1

The first scale-up cohort contains **12 new species**. Together with the completed development cohort, the cumulative resource will contain 18 species, but inferential roles remain separated:

- **new-cohort confirmatory analysis:** the 12 prospectively admitted species only;
- **cumulative analysis:** all 18 species, reported as an updated integrated estimate;
- **completed development result:** the original six-species result remains reported exactly as frozen.

The new cohort target is:

- 200 photographs per species;
- 2,400 new photographs total;
- 80 calibration and 120 held-out evaluation photographs per species;
- 960 new calibration photographs;
- 1,440 new held-out evaluation photographs.

## Prospective species admission

Species are admitted without inspecting candidate photographs, flower-colour measurements, Stage-A effect sizes or Stage-B transition surfaces.

The candidate order is inherited from `data/global_flower_colour_species_ranked.csv`. A manual evidence ledger records whether the ranked literature match supports naturally occurring within-species floral-colour variation that can plausibly be represented by a still photograph.

Allowed ledger decisions are:

- `eligible`;
- `exclude_taxonomic_overlap`;
- `exclude_breeding_only`;
- `exclude_ontogenetic_transition`;
- `exclude_irrelevant_match`;
- `hold_requires_review`.

An `eligible` record must have:

1. source-supported natural intraspecific floral-colour variation;
2. low known risk that apparent colour is primarily a flower-age transition;
3. no overlap with a taxon already represented in the completed six-species cohort;
4. no evidence limited to a cultivar, crop-breeding line, engineered phenotype or mapping population.

The cohort selector traverses eligible records by frozen rank, applying a maximum of two admitted species per family. It stops after 12 admissions. It does not optimize on observed colour geography.

## Photograph feasibility gate

After species admission, iNaturalist acquisition uses one photograph per observation and the existing outcome-blind metadata controls. Each species must independently provide 200 records after applying the frozen observer, spatial-cell, month and coordinate-quality rules.

A species that fails the photograph gate is not replaced ad hoc. Replacement, if required, follows the next eligible ranked ledger record under the same frozen rule and is recorded in a versioned amendment before any colour measurement begins.

## Measurement gate

The 12-species scale-up uses a common still-image measurement contract frozen before evaluation:

1. flower localization/segmentation is established on the 80-image calibration split;
2. continuous colour features are retained rather than forcing universal discrete morph labels;
3. unresolved or not-evaluable records remain explicit;
4. all species-specific centring and scaling parameters are estimated from calibration images only;
5. the 120-image evaluation split remains unopened until the representation and failure rules are frozen.

The new-cohort confirmatory representation must be defined without using any of the 1,440 held-out colour values.

## Spatial inference

### Stage A — new cohort

For each admitted species:

- build colour-blind spherical k-nearest-neighbour graphs from fixed coordinates;
- primary `k = 5`, sensitivities `k = 3` and `k = 8`;
- calculate within-species continuous-colour edge discontinuity;
- permute complete colour vectors strictly within species;
- use 9,999 permutations;
- combine species-specific statistics with equal species weight.

### Stage B — new cohort

Stage B is entered only if the prospectively defined new-cohort Stage-A gate passes. Detectability is computed from observation geometry only. Shared-transition strength retains the opportunity denominator

\[
A(x)=\sum_i D_i(x)
\]

and cells without sufficient opportunity remain not evaluable.

Spatial support selection is repeated using geometry only and is frozen before observed colour discontinuities are scored. The complete within-species permutation pipeline is then rerun for every null replicate.

### Cumulative 18-species analysis

The original six and new 12 may be combined only after the new-cohort result is frozen. This cumulative estimate is not allowed to replace either cohort-specific result.

## Claims and stopping rules

The scale-up can support one of three outcomes:

1. repeated within-species organization and a supported shared transition concentration;
2. repeated within-species organization without supported shared concentration;
3. no replicated evidence for within-species organization in the new cohort.

No environmental, historical or biogeographic layer enters species admission, colour measurement, Stage A or primary Stage B. Geographic-reference overlays remain post-discovery explanations and cannot rescue a failed shared-boundary gate.

## Versioned files

- contract: `docs/supporting/jbi_ch1_scaleup_contract_v1.json`;
- evidence ledger: `docs/supporting/jbi_ch1_scaleup_species_ledger_v1.csv`;
- selected cohort: `docs/supporting/jbi_ch1_scaleup_cohort_v1.csv`;
- selection manifest: `docs/supporting/jbi_ch1_scaleup_cohort_manifest_v1.json`;
- selector: `scripts/data/select_jbi_ch1_scaleup_cohort.py`.
