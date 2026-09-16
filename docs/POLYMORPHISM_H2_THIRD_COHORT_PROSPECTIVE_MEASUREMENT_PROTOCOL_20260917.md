# H2 third-cohort prospective pixel measurement and white-axis protocol

Date frozen: 2026-09-17 JST
Status: PREOUTCOME — biological image pixels, colour outcomes, H2 vectors, W, and structured-null outcomes remain unopened

## Purpose

Run one genuinely new species-disjoint and row-disjoint prospective transport test of the already fixed achromatic–chromatic (`white` versus `non-white`) geometry. This protocol does not recover or replay P500 H2 outcomes.

The target was fixed before this cohort was selected:

- palette order: `[white, yellow, orange, red, pink, magenta, purple, blue, bronze]`;
- `q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`;
- per-species contribution: `(u_i^T q_white)^2`;
- cohort statistic: `W = mean_i (u_i^T q_white)^2`;
- construction control: frozen coarse-state-preserving structured null.

No axis refit, rotation, threshold change, species replacement, row replacement, palette change, morph redefinition, or null-model change is permitted after pixel opening.

## Frozen lineage

Third-cohort selection was performed outcome-blind after excluding all 1,000 legacy high-depth species and all 500 P500-selected species.

- selected species before fresh-metadata gate: 500;
- canonical selected CSV SHA256: `4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a`;
- ordered selected manifest SHA256: `16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59`.

Fresh-metadata gate:

- gate execution HEAD: `f8ee9f65b88e5608a07cd9008f34e48a09abff88`;
- immutable metadata freeze commit: `ca69930986e39e3ea9b2d2f97ab247045ba9db0b`;
- workflow run: `35107646642`;
- artifact id: `10451429156`;
- artifact digest: `sha256:0b44c1e3b64c9930a9ccd12a0f140e9cd076d59080b695de435501a1648b6de9`;
- candidate metadata rows across all 500 selected species: 49,999;
- request errors: 0/500;
- exactly 499 species reached 100 fresh rows;
- exactly one selected species reached 99 rows and is excluded before pixels with no replacement and no biological interpretation;
- authorized species: 499;
- authorized rows: 49,900;
- authorized metadata: `data/frozen/polymorphism_h2_third_cohort_authorized_metadata_v1.csv.gz`;
- authorized metadata SHA256: `13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290`;
- authorized species: `results/polymorphism_h2_third_cohort_candidate_metadata_20260916/authorized_species.csv`;
- authorized species SHA256: `36a866b040b835e20539d318b0533bcb905cbe760d23df29e7da6587a054e593`;
- prior observation/photo exclusion union before the draw: 178,462 / 178,462 IDs.

The biological workflow must refuse to open pixels unless all hashes and denominators above match exactly and the parent metadata receipt remains `THIRD_COHORT_METADATA_GATE_PASS` with every outcome-firewall flag closed.

## Pixel-opening population

Only the 499 frozen authorized species may be measured. For each species, open exactly the 100 already frozen metadata rows. There is no re-query, second page, replacement image, replacement species, or target relaxation.

Frozen denominator before pixels:

- species: 499;
- candidate image rows: 49,900.

Download, ROI, segmentation, palette, or other measurement failures become terminal non-classified states under the unchanged measurement machine; they do not trigger replacement.

## Exact inherited measurement machine

Materialize the same location-blind machine from source commit:

`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`

using `scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh`.

This fixes ROI-v4 flower detection, EfficientSAM segmentation, the nine-colour palette, the four coarse biological morphs (`white`, `yellow_orange`, `red_pink`, `blue_purple`), all ROI/flip/palette gates, and all terminal state definitions. `mixed_uncertain` remains non-biological. No image pixels or flower masks may persist after each partition is sealed.

EfficientSAM revision and weight hashes remain the same as the qualified P500 measurement design:

- revision: `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`;
- encoder SHA256: `84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951`;
- decoder SHA256: `a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11`.

## Location-blind firewall

Generate deterministic `measurement_id` values from `photo_id` with salt:

`FCP_H2_THIRD_COHORT_PROSPECTIVE_20260917_V1`

Worker-visible fields are only `measurement_id`, `image_filename`, and `photo_license`. Source URL remains in the acquisition key; species, taxon id, coordinates, observer, observation id, and prospective rank remain sealed until all terminal results are complete.

Partition geometry is unchanged:

- 2 blind batches;
- 32 semantic shards per batch;
- 4 compute partitions per semantic shard;
- 256 terminal partitions total;
- maximum 8 parallel measurement partitions;
- no early stopping.

## Terminal measurement and support gate

Every one of 49,900 frozen rows must receive exactly one unchanged terminal state. Only `classified_four_state_morph` is classifiable.

After all 256 terminal receipts are complete, open the metadata-colour join and compute species support.

- measurement-evaluable species: `n_classifiable >= 40` among its fixed 100 rows;
- require at least 250 measurement-evaluable species.

Fewer than 250 produces `H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT`; it is not evidence against polymorphism or the white/non-white axis.

## H2 admission

Among species with `n_classifiable >= 40`:

1. Primary coarse-state second-frequency threshold = 0.10.
2. Strict sensitivity threshold = 0.20.
3. Normalize the nine palette fractions on admitted classifiable rows.
4. Apply the unchanged deterministic unlabeled two-means in Hellinger coordinates.
5. Require continuous minor-cluster fraction >= the same threshold.
6. Set `Delta_i` to the two cluster-mean composition difference; sign is irrelevant.
7. Require nonzero norm and set `u_i = Delta_i / ||Delta_i||`.

Primary H2 requires at least 20 vectors. Fewer than 20 is `H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT`.

## Fixed statistic and structured null

For the primary admitted set:

`W = mean_i (u_i^T q_white)^2`.

No axis is fitted to third-cohort outcomes.

Use exactly 999 coarse-state-preserving structured-null worlds, with the same algorithm already frozen for P500. Preserve species × coarse-morph row counts, permute normalized nine-colour rows within each frozen coarse morph across the already selected H2 species, refit the unchanged label-free Hellinger two-means, reconstruct unit displacement axes, and recalculate W.

For direct procedural comparability, retain the already frozen seeds rather than selecting new seeds after cohort construction:

- primary 0.10 seed: `20260915`;
- strict 0.20 seed: `20261015`.

Upper-tail Monte Carlo p-value:

`p = (1 + # null W >= observed W) / 1000`.

## Decision rule

Primary 0.10 test only determines the prospective verdict:

- `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED` if support passes, primary vector support >=20, and `p < 0.05`;
- `H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED` if evaluable and `p >= 0.05`;
- `H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT` if fewer than 250 species meet measurement support;
- `H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT` if measurement support passes but fewer than 20 primary vectors exist.

The strict 0.20 result is sensitivity only and cannot rescue or overturn the primary verdict.

## Durable terminal requirement

Before biological authorization, the selected-cohort synthetic qualification has already demonstrated production H2 result writing, disk read-back, stage validation, ZIP packaging, and artifact upload without biological inputs (run `35065814180`).

The biological run is not complete when W/p is calculated in memory. It is complete only when:

1. the final result is serialized;
2. re-read from disk;
3. `SUPPORT_GATE_COMPLETE -> H2_COMPLETE` is validated;
4. the result package is uploaded successfully;
5. immutable result files are committed.

A failure after pixel opening is terminal for this cohort and must not be converted into an untouched prospective confirmation by rerun or H2-only recovery. A pre-pixel technical failure may be repaired only when the execution record proves that no candidate image acquisition or pixel opening began.

## One-shot authorization and nonclaims

A separate authorization marker must be committed after this protocol and all production code/tests are frozen. That marker authorizes one biological execution only.

This test does not estimate global polymorphism prevalence and does not establish pigment chemistry, transition direction, pollinator selection, climate adaptation, or any causal mechanism. A positive result supports only recurrent achromatic–chromatic geometry under the frozen design. A negative evaluable result is evidence against prospective transport of this fixed geometry under this design, not evidence that flower-colour polymorphism lacks structure in every system.
