# RGFCA Step 8E — 42,111-species breadth measurement protocol

Date: 2026-09-11 JST
Status: prospectively fixed before Step-8E breadth image pixels are opened.

## Primary objective

Measure one outcome-blind, frozen discovery anchor photograph for every species in the 42,111-species RGFCA discovery universe for which the exact frozen anchor can still be resolved. Keep all 42,111 species in the denominator, including transport, acquisition, and measurement failures.

This is a **breadth census of observed floral-colour states**, not a census of each species' modal, canonical, or representative flower colour.

## Frozen species and anchor frame

The species universe is exactly the V1+V2 metadata-only discovery union, deduplicated by `inat_taxon_id`:

- species universe: 42,111;
- no flower-colour, climate, pollinator, literature-polymorphism, or image-pixel outcome was used to define it;
- each species has at least one frozen V1/V2 discovery observation/photo anchor;
- the anchor for each species was selected by the deterministic SHA-256 rule already frozen in Step 8C;
- no anchor may be substituted after Step-8E colour outcomes open.

Required lineage:

- `results/rgfca_42111_tiered_measurement_step8c_20260911/species_anchor_allocation.csv.gz`;
- completed Step-8D 200/200 metadata smoke gate;
- completed full Step-8D metadata resolver before pixel authorization.

## Evaluable transport gate

Before any Step-8E image byte is opened:

1. all 42,111 species remain present in the Step-8D full-resolution table;
2. no unresolved anchor is silently replaced;
3. the full resolver used the exact frozen observation/photo IDs;
4. image pixels and flower colour remained unopened during resolution;
5. at least 95% of the 42,111 frozen anchors must resolve to the exact expected photo and a current image URL for the full breadth atlas to be called transport-evaluable.

If the 95% gate fails, Step 8E may retain a bounded measurement-feasibility result but must not call the measured subset a near-complete 42,111-species breadth atlas. The threshold is fixed before the full-resolution outcome is inspected.

## Measurement target and inherited limitation

Use the existing frozen RGFCA ROI-v4 / fixed-palette measurement implementation and exact model artifacts already used in the global discovery/reserve cohorts. Do not retrain, retune, change thresholds, or choose a new segmentation model using Step-8E outcomes.

The model's independent Monarda localization-agreement gate previously did not pass. Therefore Step 8E output is described conservatively as **model-derived colour from detected visible flower regions in the selected observation**, not verified focal-species petal truth. A 42,111-species scale-up does not repair that measurement limitation.

## Location-blind execution firewall

Before pixel access, split metadata into:

- a worker packet containing only blind record IDs and ephemeral filenames;
- a sealed acquisition key containing exact source URLs and photo identities;
- a sealed metadata join key containing species, `inat_taxon_id`, source cell, geography, and prior-open status.

Measurement workers must not receive species names, taxon IDs, coordinates, observers, literature labels, C/S labels, or global-atlas outcomes. Source-URL surfaces are deleted before colour measurement. Pixels are deleted after terminal partition sealing. Coordinate/species–colour joins occur only after all required terminal partition receipts are present.

## One-run / failure policy

- one frozen anchor per species for the breadth layer;
- no favourable rerun to replace failed images;
- no alternate photo for a deleted, inaccessible, blank, or measurement-failed anchor;
- no post-outcome change to the ROI, palette, admission rules, or morph mapping;
- all failures remain in the 42,111-species denominator with explicit status;
- the already opened 1,000-species discovery/reserve cohort is flagged and is never described as a fresh validation cohort.

## Primary Step-8E estimands

### E1. Species-equal observed colour-state composition

Among successfully measured breadth anchors, summarize the fixed 12-anchor palette and the pre-existing four display morph groups with each species contributing at most one anchor. This is an observation-state composition, not species prevalence of intrinsic flower colours.

### E2. Geographic breadth surface

Join each terminal colour result to its frozen discovery equal-area cell only after measurement completion. Report cell-wise species-equal colour composition with the original 162-cell grid and retain cells with no discovered species as structural absence from the discovery frame.

### E3. Breadth saturation

Using the complete measured breadth table as the internal reference, evaluate outcome-blind hash-based nested species counts:

- 1,000;
- 2,000;
- 5,000;
- 10,000;
- 20,000;
- 30,000;
- 42,111.

Also run 200 deterministic species-equal subsampling realizations at each N. Report separately:

- global palette/four-morph cosine similarity;
- occupied-cell retention;
- cell-wise composition similarity;
- uncertainty from unresolved/failed anchors.

Do not collapse these curves into one universal minimum data requirement.

## What Step 8E cannot estimate

One breadth anchor per species cannot establish:

- whether a species is flower-colour polymorphic;
- species-level Simpson D;
- a species' modal or representative colour;
- local coexistence C*;
- spatial segregation S*;
- within-species geographic colour organization;
- causal environmental or pollinator mechanisms.

Those require the separate depth tiers or the already completed high-depth discovery/reserve cohorts.

## Depth continuation after breadth closure

The pre-colour Step-8C capacity allocation remains:

- 12,985 species at target depth 20;
- 5,316 additional species at target depth 10;
- 6,311 additional species at target depth 5;
- 9,332 additional species at target depth 2;
- 8,167 species at breadth-only depth 1 after applying the one-anchor floor.

These depths are availability strata, not outcome-selected groups. Multi-photo analyses must report detection/estimation uncertainty by depth and must not lower the previously frozen C*/S* high-depth requirements merely to increase sample size.

## Reporting boundary

The Step-8E headline, if its transport and measurement completeness gates are met, is limited to the scale and reproducibility of a species-equal global photographic flower-colour breadth atlas. It does not establish an unbiased census of all angiosperms, a universal flower-colour boundary, or focal-species petal-colour truth.
