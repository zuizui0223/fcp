# H2 third-cohort fresh-metadata gate protocol

Date frozen: 2026-09-16 JST
Status: PREOUTCOME — metadata only; image pixels and colour outcomes remain unopened

## Purpose

Freeze the exact species and row denominator for the third-cohort prospective H2 transport test without opening flower-colour outcomes. This gate inherits the P500 metadata design and is not itself a biological H2 test.

## Frozen selected cohort

Input manifest:

`results/polymorphism_h2_third_cohort_selection_20260916/selected_species_manifest.tsv`

Expected ordered-manifest SHA256:

`16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59`

Exactly 500 ranked species must be present, with unique `inat_taxon_id`, unique species identity, and ranks 1..500. Any drift stops execution before requests.

## Freshness firewall

Every observation/photo ID used previously must be excluded before selection. The exclusion union is frozen as:

1. `data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv`
2. `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`
3. `results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz`
4. `data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz`

The first three contribute the already-frozen 128,464 unique observation/photo IDs. P500 contributes 49,998 additional fresh IDs. Because P500 was itself verified disjoint from the earlier union, the required total exclusion union is exactly **178,462 observation IDs and 178,462 photo IDs**. Any union drift stops execution before the gate result can pass.

## Frozen metadata query

Inherit the P500 query unchanged:

- one iNaturalist request per selected species;
- `quality_grade=research`;
- species rank only;
- flowering annotation term `12=13`;
- unobscured/open georeferenced records;
- positional accuracy <= 5,000 m;
- photo license in `cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc`, `cc-by-nc-sa`;
- `order_by=random`;
- `per_page=200`, page 1 only;
- observer cap = 2;
- geographic maximin reduction to exactly 100 rows/species where possible;
- request interval 1.05 s;
- request timeout 45 s;
- request retries = 0.

No target relaxation and no species replacement are permitted.

## Gate

Transport failure: request-error fraction > 0.05.

Capacity failure: after successful transport, fewer than 300 of the frozen 500 species retain exactly 100 fresh rows.

Pass: request-error fraction <= 0.05 and at least 300 species retain exactly 100 fresh rows.

The exact authorized biological denominator is all and only species reaching exactly 100 fresh rows in this single bounded draw. The authorized row denominator must therefore equal `100 * n_authorized_species` and every observation/photo ID must be unique and outside the frozen exclusion union.

## One-draw rule

This is one bounded metadata draw. GitHub Actions run attempts other than attempt 1 are forbidden and must fail before any request. A completed adverse transport/capacity result is terminal for this selected cohort; it must not be rerun, repaired by species replacement, or retuned into a pass.

## Outcome firewall

This gate must not read or compute image pixels, flower colour, morph, palette class, D, H2 vectors, W, structured-null p-values, H3 predictors, or any recovered P500 H2 result. Coordinates are used only by the already-frozen geographic maximin row sampler.

## Outputs

Persist and hash:

- all fresh metadata rows from the frozen 500 species;
- species audit for all 500 species;
- exact authorized species denominator;
- exact authorized metadata rows, 100 per authorized species;
- result JSON with transport/capacity decision and outcome-firewall receipt.

A metadata PASS still does **not** authorize biological pixel opening. A separate one-shot biological execution authorization must be frozen after this exact denominator is durably committed.
