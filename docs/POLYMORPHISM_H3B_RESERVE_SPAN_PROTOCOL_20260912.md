# H3b reserve sampled-span replication protocol

Date frozen: 2026-09-12 JST

## Status at freeze

The earlier discovery Step 4 species-attribute analysis already reported a positive association between `D` and sampled geographic span (`rho = 0.1810309282`, two-sided permutation `p = 0.00079996`, Holm `p = 0.00479976`). Therefore sampled span is not a fresh discovery hypothesis.

The species-disjoint reserve cohort has not been tested for a `D`–sampled-span association as of this freeze. Reserve is the confirmatory cohort.

`D` itself has already been computed for H3a, but the reserve `D` values have not been joined to or tested against the sampled-span predictor. The predictor panel was frozen before H3a opened `D`.

## Predictor

Primary predictor: `log1p(maximum_span_km_all_measured)` from the frozen pre-outcome H3a sampling-opportunity panel.

The span is the exact maximum pairwise haversine great-circle distance among all 100 measured rows per species with finite latitude/longitude, Earth radius 6371.0088 km. It uses no colour outcome, morph state, or classifiability filter.

This is an observational sampled-span variable, not a true range-size estimate.

## Outcome

Use the same four-state Simpson polymorphism score as H3a:

`D = 1 - sum_k p_k^2`

for `white`, `yellow_orange`, `red_pink`, `blue_purple` among frozen `global_classifiable` rows, with `n_classifiable >= 40`.

Reserve eligible-species fingerprint: 363 species before tree matching. No phylogenetic tree matching is required for the primary species-level permutation test.

## Primary test

For discovery and reserve separately:

1. compute Spearman `rho(D, log1p(span))`;
2. retain the observed predictor vector fixed;
3. permute `D` labels among species 20,000 times;
4. use the two-sided Monte Carlo p-value `(1 + #(|rho_null| >= |rho_observed|))/(20001)`.

Seed: `20260912` with deterministic cohort-specific offsets.

The reserve result alone determines replication support. Discovery is calibration only.

## Prespecified sensitivities

### Finite-sample correction

Repeat the same test with

`D_unbiased = D * n_classifiable / (n_classifiable - 1)`.

### Sampling-opportunity partial rank

Total measured image count is exactly 100 for every species in both cohorts and is therefore not estimable as a control.

Residualize centered ranks of both `D` and `log1p(span)` against an intercept plus centered ranks of:

- `log1p(n_classifiable)`;
- `log1p(n_observers_all_measured)`.

The partial-rank effect is the Pearson correlation between these residuals. For a sensitivity p-value, permute the residualized outcome labels 20,000 times against the fixed residualized predictor and use the same two-sided Monte Carlo formula.

This sensitivity does not convert sampled span into a causal range-size effect.

### Phylogenetic sensitivity

H3a found no broad continuous phylogenetic signal in the fresh reserve cohort under S1-S3. Nevertheless, report a secondary rank-PGLS sensitivity using the exact frozen S1-S3 trees:

`rank(D) ~ rank(log1p(span)) + rank(log1p(n_classifiable)) + rank(log1p(n_observers))`

with `phylolm::phylolm(model = "lambda")`.

This PGLS is diagnostic and cannot overturn the primary permutation decision.

## Decision rule

Reserve supports replication only if:

- observed Spearman rho is positive; and
- the prespecified two-sided raw-D permutation p-value is `< 0.05`.

Labels:

- pass: `H3B_SAMPLED_SPAN_REPLICATION_SUPPORTED`;
- fail: `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

Add `OPPORTUNITY_ROBUST` only if the residual partial-rank effect remains positive with two-sided permutation `p < 0.05`.

The `D_unbiased` and rank-PGLS analyses are sensitivities and do not rescue a failed primary reserve test.

## Mechanistic-trait gate

The existing source-backed trait inventory remains insufficient for a strong mechanistic H3b test without changing gates after seeing coverage:

- discovery pollination mapping previously failed because only `mixed` reached the prespecified level-size gate;
- reserve source-backed final pollination is heavily dominated by `mixed` and has only 10 mapped bee cases;
- strict source-backed self-incompatibility has 49 reserve species (SC 34, SI 15), one species below the pre-existing total-coverage gate of 50;
- life form has no adequate direct source coverage in the current FCP workflow.

These traits are therefore not promoted to the H3b primary family in this analysis. No threshold is relaxed to make them testable.

## Hard non-claims

A supported sampled-span replication would not establish:

- true species range size;
- environmental heterogeneity or divergent selection;
- pollinator-driven polymorphism;
- causal effects of geographic extent;
- a life-form or mating-system mechanism.
