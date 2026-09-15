# Prospective H2 P500 pixel measurement and white-axis confirmation protocol

Date frozen: 2026-09-15 JST
Branch: `analysis/polymorphism-42111-h1-h2-gates-20260912`

## Status and purpose

This protocol is frozen **before any P500 image pixel or flower-colour outcome is opened**. It converts the already outcome-blind P500 species selection and the already completed fresh metadata gate into one bounded prospective test of the previously fixed achromatic–chromatic (`white` versus `non-white`) axis.

The biological target is not re-discovered here. It was fixed on 2026-09-12 after the legacy high-depth H2 audit:

- palette order: `[white, yellow, orange, red, pink, magenta, purple, blue, bronze]`;
- `q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`;
- per-species statistic contribution: `(u_i^T q_white)^2`, where `u_i` is the unit label-free two-mode displacement vector;
- cohort statistic: `W = mean_i (u_i^T q_white)^2`;
- construction control: the already frozen coarse-state-preserving structured null.

No axis loading, threshold, species replacement, palette, morph definition, or null model may be changed after P500 pixels are opened.

## Frozen lineage

Parent outcome-blind species selection:

- U0: 42,111 species;
- U100: 4,730 species with `after_observer_cap >= 100`;
- legacy high-depth discovery + reserve species excluded before prospective selection: 1,000;
- prospective remainder: 3,730 species;
- fixed P500: 500 species;
- P500 SHA256: `f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4`;
- selection protocol: `docs/POLYMORPHISM_H2_PROSPECTIVE_U100_SELECTION_PROTOCOL_20260913.md`.

Parent fresh metadata gate:

- metadata file: `data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz`;
- metadata SHA256: `a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595`;
- exact-100 species file: `results/polymorphism_h2_p500_candidate_metadata_20260913/full100_species.csv`;
- exact-100 species file SHA256: `568a64e692dc21f4a5f76c5eadb615a63b07d71f01bc9cd91e8496ccd0fa3f9e`;
- 499 species reached the frozen exact target of 100 fresh metadata rows;
- those 499 species contribute exactly 49,900 candidate rows;
- one P500 species contributed 98 rows and is excluded from pixel opening before outcomes, with no replacement and no biological interpretation;
- request errors: 0/500;
- metadata verdict: `P500_CANDIDATE_METADATA_GATE_PASS`.

The measurement step must refuse to run unless all hashes, counts, and the parent gate above match exactly.

## Pixel-opening population

Only the 499 species listed in the frozen `full100_species.csv` may have pixels opened. For each, open exactly the 100 already selected metadata rows. There is no re-query, second random page, row substitution, species substitution, or target relaxation.

The denominator is therefore fixed before pixels at:

- species: 499;
- candidate image rows: 49,900.

A download failure, ROI failure, segmentation failure, palette ambiguity, or other measurement failure does not cause a replacement image. It becomes the corresponding terminal non-classified state under the unchanged measurement machine.

## Exact inherited measurement machine

Use the same location-blind machine that generated the completed legacy global high-depth cohort. The measurement infrastructure is materialized byte-for-byte from source commit:

`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`

using `scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh`.

This fixes:

- ROI-v4 flower detection;
- EfficientSAM segmentation;
- the fixed nine-colour biological palette;
- the four coarse biological morphs `white`, `yellow_orange`, `red_pink`, `blue_purple`;
- `mixed_uncertain` as a non-biological terminal state, never a fifth morph;
- all ROI / flip / palette gates and terminal status definitions;
- image-acquisition and measurement code paths;
- no persistent image pixels or flower masks.

EfficientSAM weights must use the same frozen revision and SHA256 identities as the completed global workflow. The detector weight must pass the hash checks in the materialization script.

## Location-blind firewall

Before pixels open, generate deterministic blinded `measurement_id` values from `photo_id` using the frozen salt `FCP_H2_P500_PROSPECTIVE_20260915_V1`.

The worker-visible surface contains only:

- `measurement_id`;
- `image_filename`;
- `photo_license`.

The source URL is kept only in a sealed acquisition key. Species identity, taxon id, coordinates, observer identity, observation id, and sampling metadata are withheld from colour measurement.

Use the same validated partition geometry as the completed global measurement:

- 2 blind batches;
- 32 semantic shards per batch;
- 4 compute partitions per semantic shard;
- 256 terminal partitions total;
- maximum 8 parallel measurement partitions;
- no early stopping.

The metadata-colour join may open only after all 256 terminal receipts exist and all 49,900 frozen measurement IDs have exactly one terminal result.

## Terminal-result rule

Every frozen candidate row must reach exactly one of the unchanged terminal states:

- `classified_four_state_morph`;
- `not_evaluable_insufficient_flower_pixels`;
- `not_evaluable_no_biological_palette_mass`;
- `not_evaluable_ambiguous_palette_composition`;
- `not_evaluable_roi_or_flip_gate`;
- `image_acquisition_failed`.

Only `classified_four_state_morph` with one of the four biological morphs is `global_classifiable = true`. Every other terminal state is retained as `mixed_uncertain` / non-classifiable rather than deleted or replaced.

The joined table must contain the same nine fraction columns required by the frozen H2 implementation:

- `flower_fraction_white`;
- `flower_fraction_yellow`;
- `flower_fraction_orange`;
- `flower_fraction_red`;
- `flower_fraction_pink`;
- `flower_fraction_magenta`;
- `flower_fraction_purple`;
- `flower_fraction_blue`;
- `flower_fraction_bronze`.

## Postmeasurement support gate

Before calculating any prospective H2 statistic, apply the measurement-support gate inherited from the completed global high-depth measurement:

- a species is measurement-evaluable when `n_classifiable >= 40` among its fixed 100 rows;
- require at least 250 measurement-evaluable species among the 499 measured species.

If fewer than 250 species meet `n_classifiable >= 40`, the prospective H2 test is `NOT_EVALUABLE_MEASUREMENT_SUPPORT`. This is a measurement-support failure, not evidence against flower-colour polymorphism or against the white/non-white axis. No threshold relaxation or additional image acquisition follows.

## Frozen H1/H2 admission within the prospective cohort

After the postmeasurement support gate passes, construct H2 species independently within this new cohort using the already frozen rules.

For each species with `n_classifiable >= 40`:

1. H1 coarse-state admission uses only the four frozen biological morph labels.
2. Primary threshold: second-most-frequent coarse state fraction >= 0.10.
3. Strict sensitivity threshold: second-most-frequent coarse state fraction >= 0.20.
4. Within admitted rows, normalize the nine biological palette fractions and apply the unchanged deterministic unlabeled two-means in Hellinger coordinates.
5. Require the continuous minor cluster fraction >= the same threshold (0.10 primary; 0.20 strict).
6. Define `Delta_i` as the difference between the two continuous cluster mean normalized palette vectors; its sign is irrelevant.
7. Require nonzero `||Delta_i||` and set `u_i = Delta_i / ||Delta_i||`.

The primary prospective H2 test is evaluable if at least 20 species yield primary H2 vectors, matching the minimum already encoded in the frozen label-free H2 validation implementation. Fewer than 20 is `NOT_EVALUABLE_H2_VECTOR_SUPPORT`, not a biological negative. The strict sensitivity is reported when at least 20 strict vectors exist and otherwise marked not evaluable without changing the primary decision.

## Prospective target statistic and structured null

For the primary admitted set compute:

`W = mean_i (u_i^T q_white)^2`.

Do not fit or rotate an axis from the P500 data.

Use exactly 999 structured-null worlds. In each world:

1. condition on the observed prospective H2-selected species set;
2. within the prospective cohort and each frozen four-state coarse morph, permute normalized nine-colour palette rows across selected species;
3. preserve each species x coarse-morph row count;
4. refit the unchanged label-free Hellinger two-means within each species;
5. reconstruct unit displacement axes;
6. recalculate `W` against the fixed `q_white`.

Use the upper-tail Monte Carlo p-value:

`p = (1 + # null W >= observed W) / 1000`.

Random seed: `20260915` for the primary 0.10 test and `20261015` for the strict 0.20 sensitivity.

## Decision rule

The untouched prospective decision is based only on the primary 0.10 threshold:

- `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED` if the measurement-support gate passes, primary H2 vector support is >=20, and structured-null `p < 0.05`;
- `H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED` if the test is evaluable and primary structured-null `p >= 0.05`;
- `H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT` if fewer than 250 species have `n_classifiable >= 40`;
- `H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT` if measurement support passes but fewer than 20 primary H2 vectors are available.

The strict 0.20 result is a frozen sensitivity analysis. It may strengthen the interpretation if it is evaluable and passes, but it cannot rescue a failed primary test and cannot overturn a primary confirmation.

## Outcome firewall and hard nonclaims

Before the 49,900 candidate pixels open:

- P500 identities are fixed;
- measurement rows are fixed;
- `q_white` is fixed;
- `W` is fixed;
- H1/H2 thresholds are fixed;
- the structured null and seeds are fixed;
- measurement and vector-support gates are fixed.

This prospective test does not estimate global polymorphism prevalence. It does not establish pigment chemistry, evolutionary direction (white gain/loss), pollinator selection, climate adaptation, or any other causal mechanism. A positive result supports only a recurrent achromatic–chromatic geometry of within-species flower-colour polymorphism under this frozen measurement/admission design.
