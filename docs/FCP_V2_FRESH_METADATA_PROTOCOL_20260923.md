# FCP v2 fresh-metadata gate — 2026-09-23

Status: **FROZEN BEFORE FRESH-METADATA EXECUTION; PIXELS CLOSED**

## Purpose

Convert the already frozen 300-species Panel P and 300-species Panel N queues into terminal 200-species panels using metadata availability only.

This gate does not test flower colour, D, q_white, W, spatial organization, H3, or any other biological outcome.

## Frozen queues

Queue source:

- `results/fcp_v2_species_queue_20260923/panel_P_candidate_queue.csv`
- `results/fcp_v2_species_queue_20260923/panel_N_candidate_queue.csv`

Expected canonical queue hashes:

- Panel P: `6a3c7171988f0400a053ff181cdbfa9a04a3caecc9f6deaf88ee21a9404c3349`
- Panel N: `f2ee787bda61776ae52983c6b9d03c0232294974870cca7b33868f9bf050bff1`

Each queue contains 300 species. Queue order is immutable.

## Terminalization rule

For each panel separately:

1. make exactly one metadata query for every frozen queue species;
2. exclude all previously used observation/photo IDs before observer capping;
3. apply the same metadata eligibility, observer cap and geographic maximin selection used in the third-cohort fresh-metadata gate;
4. retain at most 100 rows per species;
5. in frozen queue order, take the first 200 species that provide exactly 100 fresh eligible rows;
6. if fewer than 200 species succeed within ranks 1–300, the panel is underidentified and pixel opening remains forbidden.

No queue extension is permitted after seeing fresh-metadata availability.

## Fixed metadata query

Per species:

- iNaturalist API v1 observations;
- Research Grade;
- rank = species;
- photos = true;
- geo = true;
- flowering annotation term 12 / value 13;
- positional accuracy <= 5 km;
- obscuration = none;
- allowed photo licences = CC0, CC BY, CC BY-SA, CC BY-NC, CC BY-NC-SA;
- order_by = random;
- per_page = 200;
- exactly one request per species;
- request retries = 0;
- observer cap = 2;
- final target = 100 fresh rows/species.

Request-error ceiling is 5% across the 600 frozen queue species. Exceeding the ceiling makes the metadata transport gate non-evaluable.

## Frozen prior-ID exclusion sources

All five sources are read only for `observation_id` and `photo_id`.

Source commit:

`4233e893ae6ac323365fb07310fd94222591414a`

Files / Git blob IDs:

1. `data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv`
   - blob `90d2449e948d4af0b2d437d57a30b43a22c1e2b8`
2. `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`
   - blob `bdf2628fb3f0af24b3ad64bc6031ae924f92fb73`
3. `results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz`
   - blob `73c94f98282d895e1ac7dbed2f0d1adaf9aad3c9`
4. `data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz`
   - blob `a50e4e3e2b4e9839fe46380c4e16e75f77fabbca`
5. `data/frozen/polymorphism_h2_third_cohort_authorized_metadata_v1.csv.gz`
   - blob `381ad4df31f6325a96987e31a985c025fce54358`

No colour/palette/result columns from these sources are read.

## Outcome firewall

Before this gate completes, the workflow must not open:

- image pixels;
- morph;
- palette fractions;
- D;
- q_white projection;
- W;
- spatial outcomes;
- H3 predictors.

A metadata PASS only authorizes a **separate later technical-pixel stage**. It does not itself authorize biological colour inference.

## One-shot rule

This metadata draw may execute once.

After execution:

- no API rerun to improve a panel;
- no queue reordering;
- no species substitution beyond the frozen queue;
- no relaxing the 100-row requirement;
- no changing observer cap;
- no changing licence/accuracy/flowering filters.

Transport or capacity failure is a design result, not a biological negative.
