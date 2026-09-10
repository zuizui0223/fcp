# Polymorphism species attributes — Step 4 association protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-species-attributes-step4b`

## Status before opening D

The source/coverage preflight is complete and was run behind an outcome firewall.

- FCP discovery denominator: 369 species.
- family: 345/369 resolved prospectively from frozen iNaturalist taxon IDs.
- sampled geographic span: 369/369.
- absolute sampled latitude centroid: 369/369.
- source-backed pollination: 144/369, but the prospective category gate failed (`mixed=137`, `bee=3`, `wind=3`, `other_animal=1`). Pollination is therefore **closed** and receives `p=1` in the six-family multiplicity correction.
- direct GIFT v3.2 `Life_form_1` evidence: 2/369 (`therophyte=1`, `cryptophyte=1`). The prospective life-form gate failed. Life form is therefore **closed** and receives `p=1`.

No alternative pollination mapping, life-form trait, web search, or LLM backfill is permitted after this freeze.

## Frozen outcome

Primary response: the already-frozen four-morph Gini-Simpson flower-colour diversity

`D = 1 - sum_k p_k^2`

from:

`results/polymorphism_directionality_step1_20260909/species_metrics.csv`

Predeclared finite-sample sensitivity:

`D_unbiased = n/(n-1) * D`.

The association runner must first load and validate the pre-outcome covariate panels, then open the D columns. No covariate construction is allowed after D is loaded.

## Primary family 1 — sampled geographic span

Predictor: `log1p(maximum_span_km_after_observer_cap)` from the frozen pre-outcome panel.

Interpretation: sampled geographic span/opportunity, **not true biological range size**.

Statistic: two-sided Spearman rho between D and the predictor.

Null: 20,000 deterministic label permutations of D across species. Two-sided p is

`(1 + #{|rho_null| >= |rho_obs|}) / 20001`.

Sensitivity A: replace D by `D_unbiased`.

Sensitivity B: partial Spearman correlation between ranked D and ranked span after residualizing both on ranked `n_classifiable` with intercept.

Sensitivity C: use the predeclared alternate `log1p(maximum_span_km)` predictor.

Only the primary raw p enters Holm.

## Primary family 2 — absolute latitude centroid

Predictor: `abs(mean_latitude)` calculated prospectively from all colour-blind candidate-photo coordinates.

Statistic and p-value: the same two-sided 20,000-permutation Spearman test.

Sensitivity A: `D_unbiased`.

Sensitivity B: partial Spearman controlling ranked `n_classifiable`.

Sensitivity C: descriptive signed-centroid Spearman with `mean_latitude`; it is not a separate multiplicity slot.

Only the primary raw p enters Holm.

## Primary families 3–4 — taxonomic clustering

These are called **taxonomic clustering**, never phylogenetic signal.

Run the following independently for family and genus.

### Eligibility gate

A taxonomic level is testable only if, among species with nonmissing group labels:

1. at least 10 groups contain at least two FCP species, and
2. those repeated groups jointly contain at least 30 FCP species.

Singleton groups are excluded from the clustering statistic because they contain no within-group pair.

A failed gate receives `p=1` in the six-family Holm correction.

### Statistic

For every repeated taxonomic group `g`, calculate the mean absolute pairwise difference in D:

`W_g = mean_{i<j in g} |D_i - D_j|`.

The omnibus statistic is the **equal-group mean**

`W = mean_g W_g`.

Equal group weighting prevents large families/genera from dominating solely by species count.

Lower W means stronger clustering.

### Null and p-value

Preserve the exact group memberships and group sizes. Shuffle D labels among all species that belong to repeated groups 20,000 times.

Lower-tail p:

`(1 + #{W_null <= W_obs}) / 20001`.

Effect size:

`clustering_gain = 1 - W_obs / mean(W_null)`.

Positive values mean members of the same taxonomic group are more similar in D than expected under exchangeability.

Sensitivity: repeat the observed statistic and permutation p using `D_unbiased`. Only the primary D p enters Holm.

## Primary family 5 — pollination

Closed prospectively because the source-backed category support gate failed. Set raw `p=1`; no D test is run.

## Primary family 6 — life form

Closed prospectively because direct GIFT `Life_form_1` support was 2/369 and the prospective gate failed. Set raw `p=1`; no D test is run.

## Multiplicity

The confirmatory Step 4 family has exactly six slots in this fixed order:

1. sampled geographic span
2. absolute latitude centroid
3. family taxonomic clustering
4. genus taxonomic clustering
5. pollination
6. life form

Apply Holm's step-down correction across all six raw p-values, including `p=1` for closed/failed-gate slots.

A trait family is considered supported in Step 4 only if Holm-adjusted p < 0.05. Raw p-values, effect sizes, and all predeclared sensitivities are still reported.

## Reproducibility

- permutation seed: `20260910`
- permutations per active family: 20,000
- no new images
- no new colour classification
- no post-outcome trait acquisition, remapping, category collapsing, or variable selection
- all output tables retain species identifiers and the exact pre-outcome covariates used

## Claim boundary

Even a supported association is observational. It may support a statement that polymorphism is **associated with** sampled geographic opportunity, latitude, or taxonomic membership. It does not establish adaptation, causal environmental effects, genetic mechanism, or a globally shared flower-colour boundary.
