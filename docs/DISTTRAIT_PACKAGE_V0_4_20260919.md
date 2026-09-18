# disttrait v0.4 package status — 2026-09-19

Status: validated in-repository methods package with alternative-method benchmark coverage.

Package:

`packages/disttrait/`

Version:

`0.4.0`

## 1. New in v0.4

v0.4 adds a comparator benchmark on the same 24 synthetic conditions used by the v0.3 performance surface.

Compared methods:

1. naive pooled pairwise Spearman analysis;
2. equal-species matched-null omnibus;
3. pair-count-weighted matched-null omnibus;
4. one-sided one-sample t-test of species-specific rho values.

Canonical script:

`packages/disttrait/benchmarks/comparator_surface.py`

Frozen results:

- `results/disttrait_comparator_surface_v0_4_20260919/cells.csv`
- `results/disttrait_comparator_surface_v0_4_20260919/result.json`

Regression test:

`packages/disttrait/tests/test_comparator_surface.py`

## 2. Main result

Across the six null cells with strong between-species geographic turnover but no within-species spatial trait effect:

- naive pooled false-positive fraction = **1.00 in every cell**;
- equal-species matched-null maximum false-positive fraction = **0.025**;
- pair-count-weighted matched-null maximum false-positive fraction = **0.025**;
- species-rho one-sample t-test maximum false-positive fraction = **0.05**.

Thus the dominant protection in this benchmark is **species-conditioning**, not one uniquely privileged omnibus statistic.

Under weak signal (effect size 0.8) with 4× species-level observation imbalance:

- equal-species matched-null detection = **0.35–0.675**;
- pair-weighted matched-null detection = **0.325–0.575**;
- species-rho t-test detection = **0.35–0.70**.

Equal species weighting is somewhat more powerful than pair weighting in several imbalanced weak-effect cells, but it is not universally most powerful. The simple species-rho t-test is also competitive.

## 3. Methods interpretation

The benchmark supports a narrower and stronger methodological claim than “our omnibus is best”:

> Separating within-species spatial association from between-species geographic turnover is the central inferential step. Equal-species matched-null aggregation provides a robust default when species sampling depth differs, but alternative species-conditioned summaries can perform similarly when their assumptions are satisfied.

This interpretation is preferable to claiming unique optimality.

## 4. What is still missing for a standalone methods paper

The next comparison layer should include at least one model-based alternative:

- species fixed-effect or hierarchical logistic model;
- mixed-effects spatial model;
- GAM/GP-type within-species spatial formulation where feasible.

The package also still lacks a non-flower **empirical** transport dataset.

## 5. Current validation stack

v0.4 includes:

1. unit tests;
2. synthetic integration tests;
3. exact compact FCP/RGFCA algorithm-equivalence fixtures;
4. optional full-artifact replay hooks;
5. non-flower synthetic demonstration;
6. targeted confounding benchmark;
7. 24-cell FPR/power surface;
8. alternative-method comparator surface;
9. dedicated package CI.

The package remains a reusable inference implementation and methods-paper foundation, not a claim of universal statistical superiority.
