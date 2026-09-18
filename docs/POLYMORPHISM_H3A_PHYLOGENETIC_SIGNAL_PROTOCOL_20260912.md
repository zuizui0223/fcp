# H3a phylogenetic-signal protocol

Date frozen: 2026-09-12 JST

## Question

Does the continuous species-level flower-colour polymorphism score `D` show phylogenetic structure beyond random assignment of the same observed `D` values to species?

This is an evolutionary-structure question, not a mechanism test.

## Independence status

The 369-species discovery cohort is **not an untouched confirmatory cohort for H3a** because the earlier Step 4 association family has already shown significant genus-level taxonomic clustering of `D` (`Holm p = 0.0115`), while family-level clustering was unsupported.

The species-disjoint reserve cohort has not been used in the Step 4 species-attribute association family located in the repository as of this freeze. It is therefore the primary replication cohort for H3a.

Accordingly:

- discovery = calibration / retrospective concordance;
- reserve = primary fresh replication;
- a discovery-only signal cannot support H3a.

## Outcome

Use the existing four-state Simpson-type polymorphism score

`D = 1 - sum_k p_k^2`

among rows passing the frozen `global_classifiable` rule for the four biological coarse states.

Species enter H3a when `n_classifiable >= 40`. `mixed_uncertain` and other unresolved rows are not biological states.

Primary outcome: raw `D`.

Finite-sample sensitivity:

`D_unbiased = D * n_classifiable / (n_classifiable - 1)`.

No binary polymorphic/monomorphic threshold is used for H3a.

## Phylogeny

Use the same dated-megaphylogeny system already used in the repository's 34-species paper:

- `V.PhyloMaker2`;
- time-scaled `GBOTB.extended.LCVP` backbone;
- `nodes.info.1.LCVP` node table;
- placement scenarios S1, S2 and S3;
- package source pinned to remote SHA `7af3fb5152f691af2e4ec9d5e2e467d1b50505e9`.

The tree is constructed from species names before any H3a association statistic is calculated.

### Tree preflight firewall

The preflight may inspect only:

- species identity;
- `n_classifiable` for eligibility;
- genus/family information needed to place taxa;
- V.PhyloMaker2 placement status;
- tree tip coverage.

It must not calculate `D`, `D_unbiased`, Pagel's lambda, Blomberg's K, taxon-wise mean `D`, or any H3a p-value.

### Coverage gate

Before H3a is opened, each of discovery and reserve must satisfy on every S1-S3 tree:

1. at least 90% of the frozen `n_classifiable >= 40` cohort retained as tips; and
2. at least 250 species retained.

Failure is source/placement insufficiency, not evidence of absent phylogenetic signal.

The gate passed before outcome opening: discovery retained 368/369 (99.73%) and reserve retained 341/363 (93.94%) on each S1-S3 tree. The exact trees are frozen by artifact ID/digest and per-tree SHA256 in `results/polymorphism_h3a_phylogeny_preflight_20260912/frozen_tree_manifest.json`.

## Primary signal statistic

Blomberg's `K` is the primary H3a statistic because its significance can be calibrated by tip-label randomization while preserving the bounded and zero-heavy marginal distribution of `D`.

For each S1-S3 tree separately:

1. calculate observed Blomberg `K` for raw `D`;
2. permute the observed `D` values among the fixed tree tips 9,999 times;
3. recompute `K` for each permutation;
4. calculate the upper-tail Monte Carlo p-value `(1 + #null >= observed)/(10000)`.

The permutation changes only the mapping between species and trait value. The observed `D` marginal distribution, tree, branch lengths, species count and missingness pattern remain fixed.

## Secondary signal statistic

Estimate Pagel's `lambda` on each S1-S3 tree for raw `D` and report the package likelihood-ratio test against `lambda = 0` as a secondary diagnostic.

`lambda` is not allowed to overturn the permutation-calibrated K decision by itself.

## Sampling-opportunity sensitivity

A dedicated pre-outcome audit read only `species`, `observer_id`, `latitude`, and `longitude`; it did not read `morph`, `fine_state`, `global_classifiable`, or `D`. The exact panel is frozen by artifact ID/digest and SHA256 in `results/polymorphism_h3a_covariate_preflight_20260912/frozen_covariate_manifest.json`.

The pre-outcome audit found exactly 100 measured images per species in both 500-species source cohorts. Therefore total image count has zero variance by design and is not an estimable control.

The opportunity-adjusted sensitivity is fixed as follows, with no outcome-driven covariate selection:

1. after `D` is computed, retain `n_classifiable` as the estimator-precision/usable-image control;
2. join the frozen pre-outcome `n_observers_all_measured` and `maximum_span_km_all_measured` values;
3. within each cohort, form centered ranks of `D`, `log1p(n_classifiable)`, `log1p(n_observers_all_measured)`, and `log1p(maximum_span_km_all_measured)`;
4. regress centered `rank(D)` on an intercept plus the three centered control ranks;
5. use the residual as the opportunity-adjusted sensitivity trait;
6. repeat the same 9,999-permutation Blomberg-K test on that residual for S1-S3.

Sampled geographic span is explicitly a sensitivity control because it is partly biological as well as observational.

`D_unbiased` is reported as an additional finite-sample effect-size sensitivity, but it does not replace raw `D` as the primary confirmatory outcome.

## Decision rule

### Primary H3a support

Reserve supports H3a only if raw-D Blomberg K has permutation `p < 0.05` on **all three** placement scenarios S1, S2 and S3.

- all 3 pass: `H3A_RESERVE_PHYLOGENETIC_SIGNAL_SUPPORTED`;
- 1-2 pass: `H3A_PLACEMENT_SENSITIVE_UNRESOLVED`;
- 0 pass: `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

Discovery is reported as calibration/concordance but cannot rescue a reserve failure.

### Strong robustness label

Add `OPPORTUNITY_ROBUST` only if the opportunity-adjusted K test also has `p < 0.05` on all S1-S3 reserve trees.

Failure of that sensitivity does not prove observation bias; it means the phylogenetic pattern is not separable from the frozen opportunity controls under this design.

## Interpretation

If supported, H3a permits:

> Species-level flower-colour polymorphism is phylogenetically structured in the independent reserve cohort under all three dated-tree placement scenarios.

If opportunity-robust, add that the pattern persists after adjustment for the prespecified sampling-opportunity panel.

## Hard non-claims

H3a does not establish:

- a causal genetic mechanism;
- adaptive conservatism;
- whether polymorphism was gained or lost;
- a particular pigment pathway;
- a pollination or life-form mechanism;
- global prevalence among the 42,111-species frame.

Pagel's lambda or Blomberg's K are phylogenetic-signal statistics, not mechanisms.

## H3b firewall

No life-form, pollination, climate or other ecological predictor may be chosen, recoded, filtered or discarded based on H3a results. H3b covariates require their own pre-outcome source/coverage freeze.
