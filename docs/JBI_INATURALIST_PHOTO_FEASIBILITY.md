# JBI iNaturalist photo feasibility and calibration gate

## Purpose

This document fixes the operational gate between photograph acquisition and confirmatory global flower-colour spatial inference.

The current development sample contains 6 species with 200 photographs each (1,200 photographs total). Acquisition feasibility and overlap controls have passed, but flower-colour classification has not yet begun. The next valid result is therefore not a spatial boundary result; it is blinded measurement feasibility.

## Split

For each of the six development species:

- calibration: 80 photographs;
- held-out evaluation: 120 photographs.

Totals:

- calibration: 480 photographs;
- held-out evaluation: 720 photographs.

The split must be deterministic, frozen, and hash-manifested before any flower-colour rule is tuned.

### Executable split rule

The repository implements the split in `fcp_pipeline/photo_split.py` and `scripts/data/freeze_jbi_ch1_photo_split.py`.

Within each species, assignment is determined only from a stable photograph identifier:

`SHA256("fcp-jbi-ch1-photo-split-v1", species, photo_id)`.

Photographs are sorted by this digest within species. The first 80 enter calibration and the remaining 120 enter held-out evaluation. Row order, coordinates, observer identity, observation date, environmental metadata, and all later flower-colour measurements are ignored by construction.

The source acquisition manifest is required at:

`data/frozen/jbi_ch1_photo_source_manifest.csv`

The frozen outputs are:

- `data/frozen/jbi_ch1_photo_split_v1.csv`;
- `data/frozen/jbi_ch1_photo_split_v1_manifest.json`.

The machine-readable split contract is:

`docs/supporting/jbi_ch1_photo_split_spec_v1.json`.

The CLI auto-detects common species columns (`species`, `taxon_name`, `scientific_name`, `accepted_species`) and photograph-ID columns (`photo_id`, `image_id`, `media_id`). Explicit column names can also be supplied.

Before any measurement outcome exists, the split is materialized with:

```bash
python scripts/data/freeze_jbi_ch1_photo_split.py \
  data/frozen/jbi_ch1_photo_source_manifest.csv \
  data/frozen/jbi_ch1_photo_split_v1.csv \
  data/frozen/jbi_ch1_photo_split_v1_manifest.json
```

The command fails unless there are exactly six species, exactly 200 globally unique photograph IDs per species, and no downstream colour/visibility/segmentation outcome columns in the source manifest.

CI reconstructs the split from the committed source and requires byte-identical split output plus matching stable hashes. A partial freeze in which only one or two of source/split/manifest are committed fails CI.

## Calibration sequence

Each calibration photograph passes through the following ordered decisions.

### 1. Flower visibility

Assign one of:

- `evaluable`;
- `not_evaluable`.

Failure reasons should be coded explicitly (for example occlusion, distance, blur, overexposure, non-target organ, or insufficient flower area) rather than silently discarded.

### 2. Flower-region segmentation

For evaluable photographs, define the flower region used for colour assessment. A segmentation failure remains a measurement failure and must not be forced into a colour state.

### 3. Species-specific colour coding

Colour-state definitions are calibrated within species. No global binary assumption is imposed. Species may require two or more discrete states, or may remain unsuitable for discrete coding.

### 4. Ambiguity

If the photograph cannot be assigned reliably under the frozen species-specific rules, classify it as `unresolved`. `Unresolved` is not merged with any biological colour state.

## Blindness and leakage control

During calibration, the operator or algorithm used to establish visibility, segmentation, and colour rules must not use downstream boundary outcomes or environmental/geographic explanatory layers to choose thresholds.

The 720-image evaluation set must remain unopened for rule tuning. Any rule changed in response to evaluation-set performance starts a new development version and the affected evaluation set is no longer confirmatory.

The split-generation code rejects a source manifest that already contains recognized downstream outcome columns such as flower-colour state, visibility, segmentation status, evaluability, or unresolved status. This does not prove that a file was never viewed, but it prevents a post-outcome manifest from being silently presented as the original pre-measurement split source.

## Freeze artifact

The split is frozen before calibration starts. After the 480-image calibration is complete, the measurement-rule freeze additionally records:

- species list;
- photograph identifiers assigned to calibration/evaluation;
- split-generation rule and salt;
- visibility codebook version;
- segmentation procedure/version;
- species-specific colour-state codebook version;
- unresolved/not-evaluable rules;
- software commit SHA;
- content hashes of the split and codebooks.

The held-out evaluation set is not opened until this second measurement-rule freeze is complete.

## Gate decision

The calibration gate passes only when the six species have a documented measurement rule that is stable enough to apply unchanged to the held-out photographs. A species that cannot support reliable flower visibility, segmentation, or colour-state coding is retained in the audit and marked `not_evaluable` for the relevant downstream analysis rather than being silently removed.

## Current repository boundary

As of 2026-08-28, the split specification, generator, hashes, tests, and CI contract are implemented on the shared-boundary branch. The acquired 1,200-photo source manifest itself is not present in the remote repository, so the 480/720 assignment has not yet been materialized. No flower-colour labels have been generated, and the evaluation set therefore remains unopened with respect to the new colour-measurement protocol.

## Downstream constraint

Passing photograph measurement feasibility does not itself imply non-random colour geography. Only after held-out colour states are generated can the species-conditioned random-labelling analysis be run.
