# disttrait v0.3 package status — 2026-09-19

Status: validated in-repository methods package; candidate foundation for a separate methods paper, not yet an external release.

Package:

`packages/disttrait/`

Version:

`0.3.0`

## 1. What v0.3 adds

v0.3 retains the v0.2 layers:

- species-level categorical diversity;
- observer-disjoint reliability;
- deterministic Hellinger two-mode geometry;
- fixed-contrast alignment and construction-preserving structured nulls;
- species-specific spatial organization;
- matched focal-minus-background joint-permutation nulls;
- equal-species spatial omnibus;
- distribution–spatial and partial-rank association;
- frozen FCP algorithm-equivalence fixtures;
- optional raw-artifact replay hooks;
- non-flower synthetic demonstration;
- repeated species-conditioning benchmark.

It adds a **24-cell performance surface** spanning:

- within-species effect size: 0.0, 0.8, 1.6, 2.4;
- species-level observation imbalance ratio: 1×, 4×;
- MCAR observation loss: 0%, 25%, 50%;
- 40 replicate worlds per cell;
- 20 species per world;
- 39 matched spatial-null permutations per species.

Canonical script:

`packages/disttrait/benchmarks/performance_surface.py`

Frozen results:

- `results/disttrait_performance_surface_v0_3_20260919/cells.csv`
- `results/disttrait_performance_surface_v0_3_20260919/result.json`

Regression test:

`packages/disttrait/tests/test_performance_surface.py`

## 2. Performance-surface result

All worlds contain deliberately strong between-species geographic turnover in baseline trait frequency. Under effect size 0, there is no within-species spatial trait dependence.

Across the six null cells:

- naive pooled detection fraction = **1.00 in every cell**;
- maximum species-conditioned matched-null false-positive fraction = **0.025**;
- species-conditioned median rho remains near zero.

Power increases with imposed within-species signal:

- effect 0.8: conditioned detection = **0.225–0.675**;
- effect 1.6: **0.75–1.00**;
- effect 2.4: **1.00** in every imbalance/missingness cell.

Random observation loss mainly reduces power at intermediate effects. Observation imbalance does not inflate the null false-positive fraction in this design because the omnibus gives species equal inferential weight.

## 3. What this benchmark supports

The benchmark supports the following bounded methods statement:

> When strong between-species geographic turnover is present, a species-conditioned matched-null estimand can avoid the false positive produced by naive pooled distance–trait analysis, while recovering increasing power as an imposed within-species spatial effect strengthens across moderate observation imbalance and MCAR observation loss.

## 4. What it does not support

It does not establish:

- universal type-I error calibration;
- robustness to missing-not-at-random observation;
- superiority over every spatial mixed model, GAM, GP or hierarchical alternative;
- optimality of Spearman/Jensen–Shannon geometry;
- empirical performance for arbitrary continuous or multivariate traits;
- independence from species-selection or measurement bias.

## 5. Current validation stack

v0.3 now has five layers:

1. mathematical unit tests;
2. synthetic integration tests;
3. exact compact fixtures against frozen FCP/RGFCA algorithms;
4. optional full-artifact replay hooks;
5. simulation benchmarks, including a 24-cell FPR/power surface.

The dedicated `disttrait package` GitHub Actions workflow installs the package, compiles the source and reruns the full suite.

## 6. Methods-paper readiness

A credible independent methods paper is now plausible, but two major gaps remain.

### Gap A — alternative-method comparison

The current benchmark compares species-conditioned inference mainly against naive pooling. A methods paper should add at least:

- per-species effect meta-analysis without matched nulls;
- pooled analysis with species fixed effects;
- a hierarchical/mixed-effects spatial formulation where computationally feasible;
- sensitivity to alternative within-species dissimilarity measures.

### Gap B — empirical transport

The non-flower demonstration is synthetic. A separate methods paper should include at least one external empirical dataset with repeated georeferenced individual traits that was not used to design FCP/disttrait.

Until those are completed, v0.3 should be described as a validated reusable implementation and benchmarked inference architecture, not a universally superior statistical method.
