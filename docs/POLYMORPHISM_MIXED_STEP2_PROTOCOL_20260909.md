# Polymorphism mixed-uncertain Step 2 — frozen protocol

Date: 2026-09-09 JST

## Question

When the frozen four-state flower-colour classifier returns `mixed_uncertain`, are those rows mainly technical measurement failures, or do the chromatically evaluable ambiguous rows occupy reproducible intermediate positions between the four admitted colour morphs?

This test is downstream of the Step-1 rejection of direction-specific M1–M2 variance. It is aimed only at Claim 1 (species-level polymorphism measurement), not at rescuing Claim 2.

## Important non-circularity rule

The classifier itself calls a chromatically evaluable photo `mixed_uncertain` when the dominant coarse colour group is <0.50 or the top-two coarse groups differ by <0.10. Therefore simply showing balanced top-two coarse groups would be tautological and is not an admissible result.

The primary geometry test instead uses the full nine biological palette coordinates (`white`, `yellow`, `orange`, `red`, `pink`, `magenta`, `purple`, `blue`, `bronze`) and an independent reserve cohort to define the four admitted-morph reference centroids.

## Frozen cohorts

- Discovery outcome cohort: `data/derived/global_monte_carlo_measured_photos_v1.csv`.
- Independent reference cohort: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.
- The 369-species discovery fingerprint is reconstructed exactly as in Step 1: at least 40 `global_classifiable` rows with morph in `white`, `yellow_orange`, `red_pink`, `blue_purple`.
- No new image acquisition or reclassification is allowed.

## Structural decomposition of `mixed_uncertain`

All discovery rows with `morph == mixed_uncertain` are partitioned by the already frozen `measurement_status`.

1. `not_evaluable_ambiguous_palette_composition`: positive biological palette mass but no four-state assignment. These are the only rows eligible for the intermediate-colour geometry test.
2. all other `mixed_uncertain` statuses: technical/structural measurement missingness (for example insufficient flower pixels or no biological palette mass). These are counted but are not interpreted as colour intermediates.

The result must report the exact fraction of all mixed rows belonging to each status. It is prohibited to describe all `mixed_uncertain` rows as biological intermediates.

## Independent morph geometry

From reserve rows with `global_classifiable == true` and a valid admitted morph, normalize the nine biological flower palette fractions to sum to one and compute one 9D centroid for each of the four morphs.

For each discovery `not_evaluable_ambiguous_palette_composition` row with positive finite nine-palette mass:

- let `d_vertex` be Euclidean distance to the nearest reserve morph centroid;
- for each of the six centroid pairs, project the row onto the closed line segment joining the pair;
- retain the pair with minimum residual distance `d_segment` and its interpolation coordinate `t` in [0,1].

An `interior_bridge` row is frozen as:

- best-segment `t` in [0.20, 0.80]; and
- `d_segment` <= the 95th percentile of reserve classifiable rows' distance to their own morph centroid.

Thus an ambiguous row must be both interior to a pairwise morph bridge and no farther from that bridge than the empirically observed within-morph radius in the independent reserve.

## Primary geometry null

To test whether bridge occupancy depends on the named colour coordinates rather than generic compositional spread, generate 2,000 random permutations of the nine palette-coordinate identities (seed 20260909). For each permutation, apply the same permutation to every discovery ambiguous row while keeping the independently learned reserve centroids fixed, then recompute the `interior_bridge` fraction.

One-sided p-value:

`p_bridge = (1 + number(null_fraction >= observed_fraction)) / 2001`.

The geometry component is supported only if `p_bridge < 0.05`.

## Species-level link to discrete D

For each of the frozen 369 discovery species, using all 100 measured rows as the denominator where available, calculate:

- `mixed_rate`: all `morph == mixed_uncertain` rows / raw rows;
- `ambiguous_palette_rate`: `not_evaluable_ambiguous_palette_composition` rows / raw rows;
- `bridge_rate`: rows satisfying `interior_bridge` / raw rows;
- `technical_mixed_rate`: other mixed rows / raw rows.

Reproduce the previously noted `rho(D, mixed_rate)` descriptively, without treating its exact earlier value as a gate.

Primary species-level statistic: Spearman `rho(D, bridge_rate)` across the 369 species, with a one-sided 20,000-label permutation p-value (seed 20260910). Secondary diagnostics are `rho(D, mixed_rate)`, `rho(D, ambiguous_palette_rate)`, and `rho(D, technical_mixed_rate)`.

## Frozen interpretation rule

The intermediate-colour interpretation relevant to Claim 1 is supported only if BOTH hold:

1. the observed discovery ambiguous rows are enriched on independent reserve-defined pairwise morph bridges (`p_bridge < 0.05`); and
2. species with larger D have larger bridge-row rates (`rho(D, bridge_rate) > 0` and one-sided permutation `p < 0.05`).

If only (1) holds, ambiguous photos can be described as geometrically intermediate overall, but not as the explanation for the positive D–mixed association.

If (1) fails, `mixed_uncertain` must remain structural measurement missingness and cannot be used to strengthen Claim 1.

Even if both pass, the result does **not** license the statement that discrete Simpson D is numerically underestimated. Intermediate continuous photos are not a fifth discrete morph, so any bias claim requires a separately frozen measurement model. The admissible statement is narrower: exclusion of chromatic bridge photos makes D an intentionally conservative summary of the classifiable four-state component.
