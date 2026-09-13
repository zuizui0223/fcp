# H1 observer-disjoint D reliability — frozen protocol

Date frozen: 2026-09-13 JST

## Purpose

Test whether the current continuous four-state flower-colour polymorphism score

`D = 1 - sum_k p_k^2`

is reproducible when the observations used to estimate it come from disjoint sets of observers.

This is a measurement-reliability validation of the already frozen H1 quantity. It is not a new estimate of global polymorphism prevalence and it does not alter the H1 admission thresholds.

The previously circulated approximate observer-reproducibility value (~0.971) is **not accepted as evidence for this endpoint unless its provenance can be shown to be the current four-state D recomputed from observer-disjoint samples**. Repository audit did not establish that provenance before this protocol was frozen. Therefore the present analysis is treated as the canonical direct test.

## Frozen inputs

Use the existing measured-photo cohorts only:

- discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`
- species-disjoint reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

The exact source-file SHA256 fingerprints must be checked by workflow before analysis:

- discovery: `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`
- reserve: `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`

## Frozen biological state definition

Reuse the current H1 definition without modification.

Admitted rows must:

1. have `global_classifiable == true`; and
2. have `morph` in exactly:
   - `white`
   - `yellow_orange`
   - `red_pink`
   - `blue_purple`.

`mixed_uncertain` and every unresolved / technical state remain excluded.

For any sample of admitted rows, with four-state frequencies `p_k`, calculate

`D = 1 - sum_k p_k^2`.

No binary threshold is applied to D.

## Observer-disjoint split

The split is constructed independently within each species and uses **no colour-state labels**.

1. Retain admitted rows with a non-missing, non-empty `observer_id`.
2. Count admitted rows per observer. These counts may be used only to balance sampling depth; morph identities are not inspected during assignment.
3. Sort observers by descending admitted-row count. Resolve ties deterministically by `SHA256("h1-observer-disjoint-v1|<species>|<observer_id>")`.
4. Traverse that order and assign the entire observer to the half currently containing fewer rows. If row counts are tied, assign to A before B.
5. An observer can occur in only one half for that species.

This deterministic greedy split is fixed before the split-specific D values are computed.

## Reliability-eligible species

A species enters the reliability analysis only if:

- its full admitted sample has `n_classifiable >= 40` (the frozen H1 depth gate);
- at least two distinct valid observers contribute admitted rows;
- after observer-disjoint assignment, half A has at least 20 admitted rows and half B has at least 20 admitted rows.

Species failing these conditions are `reliability_not_identifiable_at_frozen_depth`; they are not called biologically monomorphic or unreliable.

Report the complete attrition counts by reason for discovery and reserve.

## Primary estimand

For every reliability-eligible species, independently compute `D_A` and `D_B` from the two observer-disjoint halves.

The primary statistic is across-species Spearman correlation

`rho_obs = cor_rank(D_A, D_B)`.

Use a deterministic species bootstrap (`seed = 20260913`, 10,000 resamples) for a percentile 95% confidence interval.

### Frozen measurement-reliability gate

The canonical decision is made on the species-disjoint reserve cohort only. `H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED` requires all of:

1. reserve `rho_obs >= 0.80`;
2. reserve bootstrap 95% lower bound `>= 0.70`;
3. reserve Lin concordance correlation coefficient `CCC >= 0.75`.

These are deliberately demanding engineering-style validation floors for treating D as a stable species-level measurement; they are not claimed as universal biological cutoffs. Discovery is reported as calibration/diagnostic and cannot rescue reserve failure.

If one or more conditions fail, the verdict is `H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`.

## Absolute-agreement diagnostics

Report for each cohort:

- Lin's concordance correlation coefficient;
- Pearson correlation;
- median absolute `|D_A - D_B|`;
- 90th percentile absolute difference;
- mean signed difference `D_A - D_B`;
- least-squares slope/intercept of `D_B ~ D_A`.

These diagnostics do not override the frozen reserve decision rule.

## H1 gate-stability diagnostics

Within reliability-eligible species, independently calculate in each half:

- `second_fraction`, the second-largest four-state frequency;
- primary half-status: `second_fraction >= 0.10`;
- strict half-status: `second_fraction >= 0.20`.

Report percent agreement and Cohen's kappa between A and B for the 0.10 and 0.20 rules. These are diagnostics only because the original H1 gate is defined on the full admitted sample, not on 20-row halves.

## Bias and opportunity diagnostics

Report Spearman association of absolute split discrepancy `|D_A-D_B|` with:

- full `n_classifiable`;
- number of valid observers;
- split sample-size imbalance `|n_A-n_B|/(n_A+n_B)`.

These are descriptive diagnostics and cannot be used to tune eligibility thresholds after outcomes are opened.

## Cohort interpretation

- Discovery and reserve are species-disjoint.
- The observer-disjoint split is fresh for this endpoint, but the aggregate D values in both cohorts have been used in earlier analyses. Therefore this is a direct measurement-reliability validation, not an untouched prospective biological replication.
- The reserve decision is privileged only to prevent using the discovery cohort to choose a favorable reliability interpretation.

## Hard stops

1. Do not use morph labels to assign observers to halves.
2. Do not split photographs from the same observer across halves for a species.
3. Do not lower the 20-row-per-half gate after seeing reliability outcomes.
4. Do not replace Spearman/CCC decision criteria with a more favorable metric after opening results.
5. Do not treat excluded species as evidence of low D or monomorphism.
6. Do not use this reliability result as an estimator of polymorphism prevalence in the 42,111-species sampling frame.
7. Do not resurrect the unverified historical ~0.971 value as canonical evidence after this direct test is opened.
