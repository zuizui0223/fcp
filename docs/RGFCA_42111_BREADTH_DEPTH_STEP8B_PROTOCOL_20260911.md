# RGFCA Step 8B — 42,111-species breadth–depth census

Date: 2026-09-11 JST
Status: frozen before any new Step-8B image pixel is opened.

## Purpose

Treat the complete metadata-discovered RGFCA species universe, not a high-photo-count subset, as the primary sampling frame. The universe is exactly the 42,111 unique iNaturalist species in `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv`, produced by the outcome-blind V1+V2 equal-area discovery program.

No species is removed because it has few or zero currently eligible photographs after the capacity filters. Low-capacity species remain in the census as structural sampling limitations.

## Inputs

- `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv` — exactly 42,111 species.
- `data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv` — one capacity row per discovered species after transport-only recovery.
- `data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz`.
- `data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz`.
- `docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json` — must report combined species = 42,111 and no colour/pixels opened.
- `docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json` — must report discovered species scanned = 42,111 and remaining request errors = 0.

## Step 8B-A — availability census

For every one of the 42,111 species retain:

- species name and iNaturalist taxon ID;
- V1/V2 discovery membership;
- `after_observer_cap` from the recovered capacity scan;
- capacity-scan geographic span;
- deterministic availability bins;
- eligibility at depths 1, 2, 5, 10, 20, 30, 40, 50, 60, 80 and 100 raw photos/species.

Report counts at every threshold. A species with `after_observer_cap = 0` remains a member of the 42,111-species universe and is marked `structural_no_current_eligible_photo` rather than silently dropped.

## Step 8B-B — repeated breadth saturation

Use only the frozen V1+V2 discovery taxon–cell links. Do not open images or use flower colour.

For N in 500, 1,000, 2,000, 5,000, 10,000, 20,000, 30,000 and 42,111, generate 200 deterministic equal-species random orderings. For each replicate and N compute:

1. retained occupied equal-area cells / occupied cells in the full 42,111-species frame;
2. Pearson correlation between the 162-cell species-count vector for the prefix and the full-frame vector;
3. retained unique taxon–cell links / full-frame links.

Seeds are `20260911 + replicate_id`. Sampling is without replacement and independent of photograph count, colour, family, climate, geography, significance or prior literature.

This estimates how quickly geographic representation saturates as species breadth increases. It is a metadata sampling result, not a biological flower-colour result.

## Two-layer measurement interpretation

The later image program must preserve the distinction between:

### Breadth layer

Aim to measure at least one image for every discovered species for which a prospectively fixed eligible image can be obtained. Species with no eligible image remain explicit missingness. One image per species can support broad global colour-composition summaries but cannot estimate intraspecific polymorphism.

### Depth layer

Use higher-photo-count species for increasingly demanding species-level estimands:

- >=2/5 photos: descriptive colour heterogeneity only;
- >=10/20 photos: preliminary species-level polymorphism distribution with uncertainty;
- >=40 photos: existing high-depth within-species spatial/C*–S* estimands become eligible under their separate frozen rules.

Do not infer absence of polymorphism from one or a few photos.

## Claim boundaries

- 42,111 is an iNaturalist metadata-discovery frame, not all angiosperms and not 42,111 known polymorphic species.
- Capacity is an observation-process property, not a plant trait.
- Geographic coverage saturation does not establish colour-pattern saturation.
- No literature-positive selection is used.
- No colour, climate, pollinator or ecological outcome may alter the Step 8B species universe or availability thresholds.
