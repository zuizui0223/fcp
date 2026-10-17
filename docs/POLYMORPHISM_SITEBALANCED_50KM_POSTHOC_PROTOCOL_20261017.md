# Photo-density-balanced FCP geographic 50-km sensitivity — 2026-10-17

## Research question

**Post hoc falsification, not prospective confirmation.** The frozen exact-morph-composition-preserving 50-km null gives positive local flower-colour diversity depletion in three species-disjoint cohorts. Can this result survive after reducing the disproportionate weight given to abundant photographed pairs from highly sampled geographic localities?

## Immutable source and estimand

Read the exact SHA256-verified discovery, validation and third-cohort measured CSVs using the existing distributed-polymorphism source loader. Eligibility is unchanged: one of the four globally classifiable coarse photo-colour labels and >=40 classifiable photographs per named species. Preserve the original geographic 50-km radius, species-wide unordered photo-pair Gini–Simpson diversity, and 199 full-vertex label permutations per species.

Assign original photo coordinates to 0.25-degree geographic cells. Require at least three classifiable photos per contributing cell. Construct an unordered geographic cell-pair block for each combination of cell IDs among conspecific pairs <=50 km apart; include within-cell blocks. Retain blocks with >=4 photo pairs; require >=30 photo pairs within supported cells, >=2 distinct retained blocks spanning >=2 cells. Do not admit a singleton distant photo as independent geographic evidence.

Define the balanced local diversity as the mean of block-level heteromorphic photo-pair fractions, assigning exactly 1/K weight to each qualifying cell-pair block regardless of its photo count. Depletion is the species-wide photo-pair diversity minus balanced local diversity.

The null permutes complete four-state photo labels among all fixed conspecific vertices, preserving the exact photo-colour multiset and positions. For each of 199 matched null worlds, apply the original block weights and average the species depletions equally within a cohort. Compute plus-one upper-tail p, observed mean and median, 1999-iteration species-bootstrap mean 95% interval, positive-species fraction, site/block support, and the original unweighted depletion **recalculated on the same eligible species**.

Cross-cohort nominal support requires positive balanced mean, p<.05, and >=30 eligible taxa in each of the three independent species cohorts. The geographic cell width, minimum support and other parameters must not be tuned after reading the result.

## Interpretation limits

Positive evidence would demonstrate that this photographed 50-km signal cannot be fully explained by very dense photo-pair weighting inside a few geographic blocks. If it disappears or is underpowered, preserve the original numerical results as frozen, and explicitly limit paper claims about geographic colour partitioning.

These 0.25-degree geographic cells are not field-verified populations or individuals; the analysis does not control all observation, camera, measurement, coordinate error, flowering-stage, seasonal or source-provenance biases. It cannot infer selection, genotype, local adaptation, pigment cost or fitness. The separate among-species seasonal climate and white-flower analysis belongs to the mois repository and is not duplicated.

## Executable records

- scripts/analysis/run_polymorphism_sitebalanced_50km_posthoc_20261017.py
- tests/test_polymorphism_sitebalanced_50km_posthoc_20261017.py
- .github/workflows/polymorphism-sitebalanced-50km-posthoc-20261017.yml

Status on creation: synthetic regression tests passed locally; real original-photo outputs are not yet established. GitHub Actions result JSON and full per-species CSV will be the empirical authority. No pre-existing frozen result is edited.
