# H1 observer-disjoint D reproducibility — frozen protocol

Date frozen: 2026-09-13 JST

## Question

Is the continuous four-state flower-colour polymorphism score `D = 1 - sum_k p_k^2` reproducible as a species-level measurement when the contributing observers are completely disjoint?

This is a measurement-validity test. It is not a prevalence estimate and does not test ecological or evolutionary causes of polymorphism.

## Outcome firewall and ordering

1. Run `run_polymorphism_h1_observer_split_preflight_20260913.py` first.
2. The preflight may read only `species`, `observer_id`, and `global_classifiable`.
3. The preflight must not read `morph`, palette columns, D, or any H1 outcome.
4. Observer assignment and the per-half sample-size gate are fixed from this outcome-blind preflight before morph labels are opened.
5. The H1 analysis must validate the preflight firewall receipt before computing D.

Historical reliability values whose exact provenance cannot be tied to the current four-state D definition are not used as H1 evidence.

## Cohorts

- Discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`
- Fresh reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

The reserve cohort is the primary decision cohort. Discovery is calibration/replication context and cannot rescue a reserve failure.

Only species with at least 40 `global_classifiable` rows in the full cohort enter the observer-split opportunity frame. Expected fingerprints are 369 discovery species and 363 reserve species.

## Observer-disjoint split

Within each species, use classifiable rows only.

- Treat each observer as an indivisible block.
- Sort observer blocks by decreasing classifiable-row count, then by a frozen SHA256 tie-break.
- Greedily assign each whole observer block to the half with the lower accumulated classifiable count, with deterministic tie-breaking.
- An observer must never occur in both halves for the same species.
- Morph labels, palette values, coordinates, dates, and all colour outcomes are forbidden inputs to splitting.

The exact assignment produced by the preflight is an immutable analysis input.

## Frozen per-half gate rule

Candidate minimum classifiable counts per half are `{20, 15, 10}`.

Select the **largest** candidate threshold for which both discovery and reserve retain at least 100 species. This rule is evaluated from the outcome-blind opportunity receipt only.

If no candidate satisfies the rule, verdict is `H1_OBSERVER_SPLIT_FEASIBILITY_FAIL` and D must remain unopened for this test.

All candidate thresholds with at least 50 reserve species may be reported as prespecified sensitivity diagnostics, but only the automatically selected threshold determines H1.

## D construction

After the gate is fixed, open only rows assigned to the eligible species/halves and use the frozen four biological states:

- `white`
- `yellow_orange`
- `red_pink`
- `blue_purple`

Rows must be `global_classifiable == true` and `morph` must be one of these four states. `mixed`/ambiguous rows are not biological states and must not be promoted.

For each species and each observer-disjoint half independently:

`D = 1 - sum_k p_k^2`.

Finite-sample sensitivity:

`D_unbiased = D * n / (n - 1)`.

## Primary statistic

Fresh reserve Spearman correlation between independently estimated `D_A` and `D_B` across species passing the frozen per-half gate.

Also compute, without changing the decision:

- species bootstrap 95% percentile CI for Spearman rho (5,000 resamples),
- two-sided species-label permutation p-value (20,000 permutations of half-B values),
- Lin concordance correlation coefficient (CCC),
- median and 90th percentile absolute `|D_A - D_B|`,
- the same Spearman analysis for `D_unbiased`,
- discovery-cohort counterparts as calibration context.

Fixed seeds derive from `20260913`.

## Decision rule

H1 is supported only if the fresh reserve primary threshold satisfies all of:

1. `rho(D_A, D_B) >= 0.80`;
2. bootstrap 95% lower bound for rho `> 0.70`;
3. two-sided permutation `p < 0.001`.

Verdicts:

- all three pass: `H1_OBSERVER_DISJOINT_D_REPRODUCIBLE`
- otherwise: `H1_OBSERVER_DISJOINT_D_NOT_SUPPORTED`
- feasibility rule fails before D opens: `H1_OBSERVER_SPLIT_FEASIBILITY_FAIL`

CCC and absolute-error diagnostics describe calibration but cannot rescue or overturn the frozen decision.

## Hard nonclaims

Even a positive H1 does not establish:

- unbiased global prevalence of flower-colour polymorphism;
- biological independence of photographs;
- absence of image-classification error;
- evolutionary or ecological causes of D;
- direction of colour evolution;
- that the 369/363 validation cohorts represent all 42,111 species.

A positive H1 supports only that the frozen four-state D ranking is reproducible across observer-disjoint photo subsets under the high-depth validation design.
