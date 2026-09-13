# Prospective H2 P500 candidate-metadata acquisition protocol

Date frozen: 2026-09-13 JST
Branch: `analysis/polymorphism-42111-h1-h2-gates-20260912`

## Purpose

Acquire a fresh, metadata-only candidate draw for the already frozen outcome-blind P500 species panel. This is the final gate before any new image pixels or flower-colour outcomes may be opened for prospective H2 confirmation.

The prospective target is the already frozen achromatic–chromatic axis (`q_white`) and statistic `W`. Neither `W`, D, morph labels, palette values nor image pixels are used in this acquisition.

## Frozen P500 identity

- parent selection protocol: `docs/POLYMORPHISM_H2_PROSPECTIVE_U100_SELECTION_PROTOCOL_20260913.md`
- U0: 42,111 species
- U100: 4,730 species with metadata capacity `after_observer_cap >= 100`
- legacy discovery + reserve species excluded before selection: 1,000
- prospective remainder P100: 3,730 species
- frozen P500: 500 species
- P500 SHA256: `f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4`
- replacement of P500 species after outcome opening: forbidden

The P500 file must be regenerated from the frozen metadata-only selector and must match this SHA byte-for-byte before any request is sent.

## Fresh acquisition query

Use the same metadata eligibility rules that defined the U100 capacity tier:

- one iNaturalist request per P500 species;
- `per_page = 200`;
- `quality_grade = research`;
- `photos = true`;
- `geo = true`;
- `rank = species`;
- flowering annotation term `12`, value `13`;
- positional accuracy <= 5,000 m;
- unobscured/open coordinates only;
- photo licences: `cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc`, `cc-by-nc-sa`;
- `order_by = random`;
- observer cap = 2 photographs per observer;
- exact target = 100 fresh metadata rows per species;
- target relaxation is forbidden.

After observer capping, if more than 100 rows remain, select exactly 100 using the existing deterministic geographic-maximin rule. This selection uses coordinates only as sampling metadata.

## Prior-ID exclusion

Before observer capping, exclude every observation ID or photo ID present in any of the three preflight-frozen sources:

1. `data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv`
2. `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`
3. `results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz`

The metadata-only preflight found a union of 128,464 observation IDs and 128,464 photo IDs. The acquisition script must recompute these sets and refuse to run if either union count differs.

## Transport policy

- request pacing: 1.05 s minimum interval;
- timeout: 45 s;
- request retries: 0;
- each P500 species receives exactly one frozen request attempt in this draw;
- low-capacity species are not re-queried to obtain a favourable random page;
- no replacement species are admitted.

The request-error fraction is evaluated across all 500 attempts. The frozen transport ceiling is 0.05, inherited from the capacity-scan contract. A larger error fraction makes the draw transport-not-evaluable; it is not a biological negative.

## Pre-pixel gate

Let `N100` be the number of frozen P500 species with exactly 100 fresh selected metadata rows after prior-ID exclusion, observer cap and geographic-maximin selection.

Pixels may be authorized in a later, separate step only when both conditions hold:

1. request-error fraction <= 0.05;
2. `N100 >= 300`.

Verdicts:

- `P500_CANDIDATE_METADATA_GATE_PASS` — both conditions hold;
- `P500_CANDIDATE_METADATA_TRANSPORT_NOT_EVALUABLE` — request-error fraction > 0.05;
- `P500_CANDIDATE_METADATA_CAPACITY_NOT_EVALUABLE` — transport is adequate but `N100 < 300`.

No threshold or species replacement may be changed after seeing this result.

## Outcome firewall

During this step:

- image pixels opened = false;
- flower-colour classification opened = false;
- four-state morph opened = false;
- nine-colour palette opened = false;
- D opened = false;
- H2 `W` opened = false;
- H3 predictors are not used.

Coordinates, observer identity, dates, taxon identity, licences and image URLs are sampling metadata and may be frozen for the later image-acquisition step.

## Outputs

- `data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz` — exact selected metadata rows (up to 100 per P500 species), no pixels;
- `results/polymorphism_h2_p500_candidate_metadata_20260913/species_audit.csv` — reason-coded per-species acquisition audit;
- `results/polymorphism_h2_p500_candidate_metadata_20260913/full100_species.csv` — species that passed the exact 100-row gate;
- `results/polymorphism_h2_p500_candidate_metadata_20260913/result.json` — immutable gate receipt and hashes.

A gate failure is retained as a result; it must not be converted into monomorphism or absence of H2 structure.
