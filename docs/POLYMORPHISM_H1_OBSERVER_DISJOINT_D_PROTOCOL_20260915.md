# H1 observer-disjoint four-state D reproducibility — frozen protocol

Date frozen: 2026-09-15 JST

## Purpose

Test whether the continuous four-state flower-colour polymorphism score

`D = 1 - sum_k p_k^2`

behaves as a reproducible species-level measurement when photographs from the same species are separated by observer identity.

This is a measurement-validation analysis. It does **not** estimate global polymorphism prevalence and it does not test an ecological mechanism.

A historical project summary reported observer-disjoint reproducibility near 0.97, but an exact current receipt tying that number to the present four-state D definition could not be identified in the repository audit. That historical value is therefore not used as H1 evidence here. The present analysis directly reconstructs the current four-state D.

## Data

Use the frozen high-depth measured-photo tables only:

- discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`
- reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

The two cohorts are species-disjoint. Discovery is calibration/history; reserve is the primary replication cohort.

No new image acquisition is allowed.

## Admitted observations

Retain rows only when:

1. `global_classifiable` is true;
2. `morph` is one of the four frozen biological states: `white`, `yellow_orange`, `red_pink`, `blue_purple`;
3. `species` is non-empty;
4. `observer_id` is non-missing.

`mixed_uncertain` and all unresolved/non-classifiable rows remain excluded.

## Outcome-blind observer split

Observer assignment must be constructed without using `morph`, colour proportions, D, palette values, geography, date, or any downstream outcome.

For each species:

1. count admitted rows per `observer_id`;
2. order observer groups by decreasing admitted-row count, with deterministic SHA256 tie-breaking using the fixed salt `fcp-h1-observer-disjoint-d-v1`, cohort, species, and observer ID;
3. greedily assign each whole observer group to the currently smaller split; on an exact size tie, use the fixed SHA256 tie-break;
4. never divide observations from one observer across the two splits.

The use of admitted-row counts is permitted because it contains measurement-depth information only; colour-state identities are not consulted during assignment.

## Species eligibility for reproducibility

A species enters the observer-disjoint reproducibility analysis only if:

- full admitted `n_classifiable >= 40`;
- at least two unique observers contribute admitted rows;
- both observer-disjoint halves contain at least 20 admitted rows after deterministic assignment.

The 20-row half threshold is fixed as one half of the existing 40-row H1 depth gate. Species failing this split-depth requirement are `not_evaluable_for_observer_split`, not biologically monomorphic.

Before D is computed, report the number of eligible species and require at least 100 reserve species. If reserve eligibility is below 100, verdict is `H1_OBSERVER_DISJOINT_NOT_EVALUABLE` and no biological/reliability negative is inferred.

## Four-state D

Within each eligible species and split separately, calculate the four state frequencies `p_k` and

`D = 1 - sum_k p_k^2`.

Also calculate the finite-sample sensitivity

`D_unbiased = D * n / (n - 1)`.

No binary D threshold is introduced.

## Primary reserve reproducibility

Primary statistics in the species-disjoint reserve are:

1. Spearman correlation between `D_A` and `D_B`;
2. Lin concordance correlation coefficient (CCC) between `D_A` and `D_B`.

Use 10,000 species bootstrap replicates for 95% intervals. A two-sided 20,000-label permutation p-value for Spearman is reported as a calibration of nonzero association, not as the reliability magnitude criterion.

The frozen support rule is:

- reserve eligible N >= 100;
- reserve Spearman rho >= 0.80; and
- reserve CCC >= 0.80.

If all three hold, verdict is `H1_OBSERVER_DISJOINT_REPRODUCIBILITY_SUPPORTED`.
If reserve N >= 100 but either magnitude criterion fails, verdict is `H1_OBSERVER_DISJOINT_REPRODUCIBILITY_NOT_SUPPORTED`.

Discovery is reported as calibration and cannot rescue reserve failure.

## Sensitivity and diagnostics

Report without changing the primary decision:

- the same statistics for `D_unbiased`;
- Pearson correlation;
- mean and median absolute `|D_A-D_B|`;
- split sizes and observer counts;
- fraction of the original n>=40 species that are observer-split evaluable;
- discovery-cohort statistics under the identical algorithm.

## Interpretation ceiling

A positive result supports only this statement:

> Under the frozen high-depth photo design, the current four-state D is reproducible across disjoint sets of observers at the species level.

It does not show that the four coarse states are genetically discrete morphs, that the observed D equals true population morph diversity, that 369/363 species are globally representative, or that measurement missingness is ignorable.

## Hard stops

1. Do not choose observers or split assignments using morph identities or D.
2. Do not relax the 20-per-half or 100-reserve-species gates after outcomes are opened.
3. Do not promote `mixed_uncertain` to a fifth state.
4. Do not use discovery to rescue a reserve failure.
5. Do not reuse the historical ~0.97 summary as a substitute for this direct current-D test.
