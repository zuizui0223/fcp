# RGFCA independent global-axis replication — metadata recovery acquisition protocol

Date frozen: 2026-09-09 JST

Parent ecological protocol: `docs/RGFCA_INDEPENDENT_GLOBAL_SIGNAL_AXIS_REPLICATION_PROTOCOL_20260909.md`, frozen at commit `6ea2783970d0e0fdea959a4ce13abaafc4e9b5d9`.

This protocol is a **metadata-support recovery lane only**. It does not modify the ecological estimand, the 300-species floor, the 300-km lattice, the 500-km separation floor, the 6-photo / 3-observer support rule, the measurement gate, or any M1-M3 success criterion. It does not authorize image download or ecological outcome opening.

## Why this recovery exists

The first metadata-only prequalification used the already-frozen 1,000-species candidate pool and failed prospectively: after excluding the 500 discovery/reserve species and currently known opened-frame identities, only 42 species met the frozen spatial support rule, versus the required 300. No candidate image pixels, flower colour, M1-M3 scores, climate or island information were opened.

That failure is a support/acquisition failure, not evidence against M1-M3 recurrence. Thresholds are therefore not relaxed. Instead, this lane uses capacity-selected species that were never queried by the first candidate acquisition and applies identity filtering **before** spatial-support selection.

## Frozen source population

The only admissible source population is:

- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`
- frozen SHA-256: `850501d7293a3968cf6566f8e9d7a396afd3927e8166ab2468e58ed22afcabde`
- source denominator: 4,730 species at the already-selected 100-photo capacity target.

Exclude every `inat_taxon_id` present in:

- `data/frozen/global_monte_carlo_candidate_species_audit_v1.csv`
- frozen SHA-256: `c9e9e76883d18bbb06355e1dec46a83aaf8f9a8680fb9b2c99d3cc811d62b152`
- denominator: 1,000 previously queried candidate species.

The expected pre-span recovery source is therefore exactly **3,730 species**. If the set difference is not exactly 3,730, stop before any API request.

### Necessary span prefilter

The parent ecological protocol requires at least one pair of qualifying spatial-group centroids separated by >=500 km. Therefore a capacity-selected species whose already-frozen metadata-only `maximum_span_km < 500` cannot possibly qualify.

Before requests, exclude only those impossible taxa. This is not a relaxation or enrichment by colour outcome; it is a logical necessary-condition filter using metadata frozen before ecological pixels.

Retain all remaining source taxa. Do not rank or select them by flower colour, M1-M3, climate, island status, taxonomic affinity, hotspot membership or prior ecological result.

## Exact query contract

For each retained taxon, reuse the old candidate-acquisition query semantics exactly unless stated below as transport-only execution detail:

- stable order: `order_by=id`, `order=desc`;
- `per_page=200`;
- page-1 count probe;
- maximum API page 50;
- exactly 3 deterministic candidate pages when at least 3 pages exist, otherwise all available pages;
- page ranking seed `20260917` using the already-implemented SHA-256 page rule;
- no early success stop;
- no extra page after underperformance;
- request retries: 0;
- quality grade: research;
- photos required;
- georeferenced observations required;
- taxon rank: species;
- flowering annotation term `12=13`;
- positional accuracy <=5000 m;
- unobscured coordinates;
- allowed photo licences: `cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc`, `cc-by-nc-sa`.

No image URL is dereferenced in this lane. Only API metadata are parsed.

### Transport rate

Use 8 deterministic shards, assigned after sorting by ascending `inat_taxon_id`, then species, with zero-based row index modulo 8. At most 2 shards may execute concurrently, each with a >=1.30-second request interval. This is a transport-only rate setting and does not alter taxon membership, pages or selection. All shards are required before reduction.

If the aggregate request-error fraction exceeds 0.05, the recovery is not evaluable and no biological or metadata-support conclusion is promoted.

## Identity firewall before spatial support

Before a queried metadata row can contribute to spatial support, exclude known overlaps with all identity sources currently recoverable without image opening.

At minimum exclude:

1. all reserve `species` / `inat_taxon_id`, `photo_id`, `observation_id` and exact `observer_id` from Actions artifact `10025692225` (`rgfca-reserve-measured-photos-v1`, digest `sha256:fbe4e8ab93e665115b5560cd05ddd31d02b9c9138b8a30ce460b61b6c6775ef5`);
2. all H9 fresh `photo_id`, `observation_id`, exact `observer_id` from `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`;
3. all `photo_id` and `observation_id` from `data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv`.

The legacy pre-H7 ledger does not retain observer IDs. Therefore this recovery may establish spatial prequalification but **cannot set the final all-opened-frame observer identity gate to true**. Missing historical observer identity must remain an explicit downstream blocker; it must not be treated as absence of overlap.

After identity exclusion, remove duplicate observation and photo IDs before support calculations.

## Frozen spatial support calculation

Use the exact fixed 0-shift discovery lattice:

- CRS: EPSG:6933;
- cell width: 300,000 m;
- assignment: `gx=floor(x/300000)`, `gy=floor(y/300000)`;
- no random lattice shift.

A support group qualifies only with:

- >=6 distinct photographs;
- >=6 distinct observations;
- >=3 distinct observers.

A species spatially qualifies only with:

- >=3 qualifying groups; and
- maximum pairwise distance between qualifying-group centroids >=500,000 m.

The query and these support rules are applied once. No page, cell, observer or species threshold may be changed after the recovery result is known.

## Retained metadata for qualifying species

For every spatially qualifying species, retain at most 5 qualifying groups using the same deterministic dispersion rule already implemented for the first prequalification: farthest pair first, then farthest-point maximin until 5 groups or all groups are retained.

Within each retained group freeze exactly 6 metadata rows by deterministic observer-balanced round robin, with max 2 photos from any observer in that group.

These retained rows are metadata manifests only. Their image pixels remain closed.

## Global reduction and 300-species decision

After all shards are frozen, pool all spatially qualifying species and rank them exactly as in the parent protocol:

1. descending qualifying spatial-group count;
2. descending identity-clean distinct observer count;
3. descending maximum qualifying-group centroid separation;
4. ascending stable taxon ID.

If >=300 species spatially qualify, retain exactly the first 300 and their pre-frozen support groups/photo metadata.

If <300 qualify, record `spatial_recovery_prequalification_pass=false` and stop. Do not query additional pages, lower support floors, add previously queried taxa, or modify the lattice after seeing the result.

Even if >=300 qualify, keep:

- `final_identity_gate_pass=false` until all historical opened/model-validation observer identities are reconstructed and checked;
- `measurement_gate_passed=false` until the independent human flower-reference gate passes;
- `opening_authorized=false`.

## Forbidden operations

This recovery lane must not:

- download candidate images;
- open candidate image pixels;
- calculate flower/background palette values;
- project M1-M3 scores;
- use climate, island/mainland status, hotspot membership or candidate geography beyond the fixed support lattice;
- use genus/family or known colour biology to choose taxa;
- rescue failure by extra pages or altered thresholds.

The only admissible output is a frozen metadata-support result plus an explicitly sealed candidate manifest for any species that pass.