# H1 observer-disjoint D reliability — frozen protocol

Date frozen: 2026-09-13 JST

## Question

Can the continuous four-state flower-colour polymorphism score `D = 1 - sum_k p_k^2` be recovered reproducibly from observation sets contributed by completely disjoint observers?

This is a measurement/sampling-stability test. It is not a test of global polymorphism prevalence, biological mechanism, or correctness of the image-level colour classifier.

## Frozen inputs

Two species-disjoint high-depth cohorts are used without changing their existing image measurements:

- discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`;
- reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.

The workflow must verify the already frozen SHA256 fingerprints before analysis:

- discovery: `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`;
- reserve: `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`.

Only the frozen biological states `white`, `yellow_orange`, `red_pink`, and `blue_purple` are used. `mixed_uncertain` is never promoted to a biological morph.

## Species eligibility

For the observer-disjoint analysis, rows without a usable `observer_id` are excluded rather than imputed.

A species enters the split procedure if, among observer-known rows:

1. at least 40 rows are `global_classifiable` and belong to one of the four frozen biological states; and
2. at least two distinct observers are present.

This observer-known eligibility is separate from the existing whole-cohort 369/363 D fingerprints; attrition is reported, not repaired.

## Outcome-blind observer partition

For each cohort, 200 deterministic partitions are generated with seeds `20260913 + 1 ... 20260913 + 200`.

Within each species and partition:

1. observers are kept intact; no observer can occur in both halves;
2. observers are ordered by descending number of all measured rows contributed to that species;
3. ties are broken by SHA256 of `seed | species | observer_id`;
4. observers are greedily assigned to the currently smaller half by all-row count; ties between halves are resolved by the same frozen hash order.

The balancing step therefore uses observer identity and row counts only, not morph labels, D, classifiability, colour composition, geography, or any H2/H3 result.

## Split-half D

For each half independently, retain only `global_classifiable` rows in the four frozen states and calculate

`D_half = 1 - sum_k p_k^2`.

### Primary precision gate

A species contributes to a partition-level primary reliability estimate only if each half contains at least 20 classifiable rows.

### Sensitivity precision gate

Repeat the same summaries requiring at least 15 classifiable rows per half.

Thresholds are frozen before opening observer-disjoint D agreement.

## Primary reliability statistic

For each partition, calculate Spearman correlation across paired species between `D_A` and `D_B`.

The primary cohort summary is the median split-half Spearman rho over the 200 partitions, with the 5th and 95th percentiles reported.

Because each half uses about half of the information in the full estimate, also report the Spearman-Brown projection

`R_full = 2*rho_split / (1 + rho_split)`.

This projection is descriptive; the decision rule is stated equivalently on split-half rho below.

## Frozen decision rule

The reserve cohort is the primary decision cohort. H1 reliability is supported only if, under the primary >=20-per-half gate:

1. the median number of paired species across partitions is at least 100;
2. median split-half Spearman rho is at least `2/3`, corresponding to Spearman-Brown full-estimate reliability >=0.80; and
3. the 5th percentile of split-half rho is at least 0.50, corresponding to full-estimate reliability >=2/3.

Verdicts:

- all three pass: `H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`;
- paired-species median <100: `H1_OBSERVER_DISJOINT_D_RELIABILITY_UNDERIDENTIFIED`;
- otherwise: `H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`.

Discovery is reported as a transport/calibration cohort and cannot rescue reserve failure. If discovery also passes the same numerical rule, append the diagnostic label `DISCOVERY_CONSISTENT`.

## Agreement diagnostics

For every partition also report, without changing the primary decision:

- Lin-style concordance correlation coefficient (CCC);
- mean absolute difference `mean(|D_A-D_B|)`;
- signed mean difference `mean(D_A-D_B)`;
- Spearman-Brown projected reliability.

The >=15-per-half sensitivity cannot rescue failure of the primary >=20 gate.

## Hard non-claims

A positive result would support reproducible ranking/stability of the four-state species-level D score in the existing high-depth design. It would not establish:

- unbiased global prevalence of flower-colour polymorphism;
- that 369/363 species are globally representative;
- perfect image-level state classification;
- absence of geographic morph-frequency structure;
- causal ecological or phylogenetic mechanisms;
- H2 achromatic-chromatic geometry, which has its own frozen structured-null evidence chain.

No threshold, partition rule, state definition, cohort role, or decision criterion may be changed after the observer-disjoint agreement results are opened.