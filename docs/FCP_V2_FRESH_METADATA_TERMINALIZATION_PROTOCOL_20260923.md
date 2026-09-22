# FCP v2 fresh-metadata terminalization protocol

Date frozen: 2026-09-23 JST  
Status: **PRE-PIXEL / OUTCOME-BLIND**

## Purpose

Convert the already frozen 300-species Panel P and 300-species Panel N queues into exact terminal measurement panels using fresh iNaturalist metadata only.

This gate does not open image pixels and does not read flower colour, morph, D, q_white/W, spatial outcomes or H3 predictors.

## Frozen queue inputs

- Panel P: `results/fcp_v2_species_queue_20260923/panel_P_candidate_queue.csv`
- Panel N: `results/fcp_v2_species_queue_20260923/panel_N_candidate_queue.csv`
- each queue contains exactly 300 species;
- queue order is immutable;
- terminal target is 200 species per panel;
- target rows are exactly 100 fresh photo rows per terminal species.

No queue extension is allowed after metadata retrieval.

## Freshness exclusion union

Before row selection, exclude every observation/photo ID used by the earlier FCP/RGFCA measurement programme.

The exclusion sources are restored read-only from commit
`f7582767da237aeb344dfc072452c6daf92e662b`:

1. `data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv`
2. `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`
3. `results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz`
4. `data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz`
5. `data/frozen/polymorphism_h2_third_cohort_authorized_metadata_v1.csv.gz`

The first four formed the previously frozen 178,462-ID exclusion union. The third cohort added 49,900 source-byte-verified rows that were already disjoint from that union. The expected combined exclusion denominator is therefore **228,362 unique observation IDs and 228,362 unique photo IDs**. Any drift stops execution.

## Frozen metadata query

For every queued species, exactly one bounded metadata request is attempted using the inherited prospective query:

- iNaturalist Research Grade;
- species rank;
- flowering annotation term 12, value 13;
- open/unobscured georeferenced observations;
- positional accuracy <=5,000 m;
- allowed licences: CC0, CC BY, CC BY-SA, CC BY-NC, CC BY-NC-SA;
- page 1 only;
- `per_page=200`;
- `order_by=random`;
- observer cap = 2;
- deterministic geographic maximin reduction;
- target = exactly 100 rows/species;
- request interval = 1.05 s;
- timeout = 45 s;
- retries = 0.

Previously used observation/photo IDs are excluded before observer capping and geographic reduction.

## Panel terminalization rule

For each panel independently:

1. preserve frozen `queue_rank`;
2. mark species that produce exactly 100 fresh eligible rows;
3. among those species, take the first 200 by frozen queue rank;
4. freeze those 200 species and their 20,000 rows;
5. do not use the remaining queue species for replacement after pixel opening.

This is a **pre-pixel availability gate**, not biological replacement.

## Gate

Each panel is evaluable only if:

- request-error fraction <= 0.05 across its 300 fixed attempts; and
- at least 200 queued species produce exactly 100 fresh rows.

Both panels must pass.

Terminal PASS requires:

- Panel P = 200 species / 20,000 rows;
- Panel N = 200 species / 20,000 rows;
- combined = 400 species / 40,000 rows;
- zero cross-panel species/taxon/photo/observation overlap;
- zero overlap with the frozen 228,362 prior-ID exclusion union.

If either panel has fewer than 200 full species, v2 is `UNDERIDENTIFIED_PREPIXEL`; the queue is not extended.

## One-draw rule

This metadata draw may run only on GitHub Actions attempt 1. Completed metadata results are terminal for these queues.

Forbidden after completion:

- rerunning requests to improve capacity;
- reordering either queue;
- extending rank >300;
- changing observer cap or maximin rules;
- lowering the exact-100 requirement;
- replacing species after pixels are opened.

## Outcome firewall

This gate must not read or compute:

- image pixels;
- flower/background colour;
- morph;
- palette fractions;
- D;
- q_white projection;
- W;
- spatial colour organization;
- H3 variables.

Coordinates and observer IDs are used only by the frozen metadata row-selection algorithm.

## Outputs

Persist and hash:

- all fresh metadata rows collected for all 600 queued species;
- species audit for all 600 queued species;
- terminal Panel P 200-species manifest;
- terminal Panel N 200-species manifest;
- combined 400-species manifest;
- exact 40,000-row authorized metadata table;
- transport/capacity result JSON.

A PASS still does not itself authorize biological pixel opening. Pixel opening requires a separate authorization after the exact 40,000-row denominator is committed and read back.
