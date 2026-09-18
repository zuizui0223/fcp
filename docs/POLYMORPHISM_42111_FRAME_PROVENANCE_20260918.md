# 42,111-species opportunity-frame provenance — 2026-09-18

Status: reporting-only provenance note for the active flower-colour polymorphism paper. This document does not recompute the species frame or alter any biological result.

## Purpose

Document how the 42,111-species global opportunity frame used by the current polymorphism paper was constructed upstream in RGFCA, and distinguish that metadata-discovery universe from the later high-depth inferential cohorts.

## Upstream source state

Authoritative upstream branch:

- `analysis/polymorphism-42111-h1-h2-gates-20260912`

Source Git objects:

- `docs/GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md` — blob `205bf8ec61155cc591bb0ea19f303e7deaffcdc3`
- `docs/REPEATED_GLOBAL_FLOWER_COLOUR_ATLAS_METHOD.md` — blob `8173e418f71fb35ef9ff74167b5ce35ebeb2395c`
- `docs/RGFCA_42111_BREADTH_DEPTH_STEP8B_PROTOCOL_20260911.md` — blob `e496ad013161e893a6b8844bc486b936b9786f59`
- `results/rgfca_42111_breadth_depth_step8b_20260911/result.json` — blob `5c8e1e9d983339ce29c7cfc323534c08aff8e59d`
- `docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json` — blob `3dfbe41375e72b8b59e1b9ad80b7e17357a1dccf`

## Species-discovery design

The original RGFCA metadata discovery used an 18 × 9 equal-area global grid (162 cells). The original one-pass baseline produced 8,989 species and was retained as round 0 rather than treated as the complete species frame.

Before new colour pixels were opened, species discovery was expanded with 20 metadata-only rounds over the same 162-cell grid. Each cell-round attempted one 200-record iNaturalist page using the same research-grade, flowering-annotation, georeferenced, positional-accuracy and licence filters as the original atlas. This produced 3,240 fixed fresh request attempts.

A technical audit showed that adjacent random-page requests could be nearly duplicated by upstream page persistence. The valid V1 discoveries were retained, but V1 odd/even overlap was not treated as evidence of independent repeated discovery. A cache-resistant V2 was therefore frozen using stable ID ordering, deterministic cell-specific pages and explicit exclusion of V1 observation IDs.

The completed V2 manifest records:

- 20 V2 rounds;
- 3,240 request attempts;
- 0 request errors;
- 393,992 API-returned records;
- 218,357 records retained after V1-ID exclusion;
- 205,955 unique V2 observations;
- 30,910 V2 species;
- 11,718 species new beyond V1;
- 30,393 species retained from V1;
- deduplicated V1 + V2 union = **42,111 species**.

No candidate image pixels or flower-colour outcomes were used to define this union.

## Frozen 42,111-species frame

The complete metadata-discovered universe is frozen as:

`data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv`

with SHA256:

`5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`.

The Step-8B census retained all 42,111 species regardless of current image capacity. Species with zero currently eligible images remained in the universe as structural sampling limitations rather than being silently dropped.

The same census recorded 85,337 discovery taxon-cell links and 128 occupied equal-area cells out of 162.

## High-depth opportunity layer

A separate metadata-only capacity scan evaluated how many of the 42,111 discovered species could support fixed high-depth photographic sampling after the observer cap.

Counts included:

- >=1 retained photo: 41,322 species;
- >=10: 18,301;
- >=20: 12,985;
- >=40: 8,753;
- >=60: 6,770;
- >=80: 5,558;
- >=100: **4,730 species**.

Thus the later high-depth opportunity ceiling is:

`U100 = 4,730 species`.

This is an observation-capacity subset, not a biological polymorphism subset.

## Relation to the current paper

The 42,111-species frame is used as a broad outcome-blind opportunity universe. It is not the denominator for estimating global polymorphism prevalence.

The original discovery/reserve polymorphism cohorts came from the earlier frozen RGFCA 1,000-species high-depth resource and are not a random sample of the 42,111 species.

For the later prospective H2 expansion:

- U0 = 42,111 species;
- U100 = 4,730 species;
- exclude legacy high-depth discovery/reserve = 1,000 species;
- P100 = 3,730 species;
- exclude all prospectively selected P500 species = 500;
- third-cohort candidate universe = **3,230 species**;
- deterministically select 500 species outcome-blind;
- fresh metadata yielded 499 authorized species with exactly 100 rows each.

## Claim boundary

The 42,111-species frame supports statements about the metadata-discovered iNaturalist opportunity universe and the availability of high-depth sampling.

It does not establish:

- the number of angiosperm species globally;
- global flower-colour polymorphism prevalence;
- colour frequencies for all 42,111 species;
- absence of polymorphism in low-capacity species;
- unbiased geographic sampling of plants;
- independence from iNaturalist observation processes.
