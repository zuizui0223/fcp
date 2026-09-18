# disttrait v0.5 package status — 2026-09-19

Status: validation candidate with model-based comparator coverage.

Package:

`packages/disttrait/`

Candidate version:

`0.5.0`

## 1. New in v0.5

v0.5 extends the v0.4 comparator benchmark with a model-based species-conditioned alternative:

- species fixed-intercept logistic model;
- one common nonnegative within-species local-position slope;
- species intercepts profiled out separately for each candidate slope;
- one-sided profile-likelihood ratio test using the 0.5 × chi-square(1) boundary tail.

The benchmark uses the same 24 conditions as v0.3/v0.4:

- effect size: 0.0, 0.8, 1.6, 2.4;
- observation imbalance: 1×, 4×;
- MCAR missingness: 0%, 25%, 50%;
- 40 replicate worlds per cell;
- 20 species per world.

Canonical script:

`packages/disttrait/benchmarks/model_comparator_surface.py`

Frozen outputs:

- `results/disttrait_model_comparator_surface_v0_5_20260919/cells.csv`
- `results/disttrait_model_comparator_surface_v0_5_20260919/result.json`

Regression test:

`packages/disttrait/tests/test_model_comparator_surface.py`

## 2. Main result

Across the six null cells:

- equal-species matched-null maximum FPR = **0.025**;
- pair-count-weighted matched-null maximum FPR = **0.025**;
- species-rho t-test maximum FPR = **0.05**;
- fixed-effect logistic maximum FPR = **0.10**;
- fixed-effect logistic mean FPR across null cells = **0.05**.

For weak signal (effect size 0.8):

- fixed-effect logistic detection = **0.975–1.00**;
- equal-species matched-null detection = **0.225–0.675**.

For effect size 1.6:

- fixed-effect logistic detection = **1.00** in all cells.

## 3. Interpretation

The fixed-effect logistic comparator is deliberately well matched to the binary-logistic data-generating process. Its high power therefore demonstrates what can be gained when a correct trait model is available.

The matched-null approach is less powerful under this correctly specified binary model but has two advantages in the same benchmark:

1. lower worst-cell null rejection;
2. fewer assumptions about the trait response distribution and common parametric effect form.

The preferred methodological claim is therefore not that one method dominates.

> Species-conditioning is the essential protection against between-species geographic confounding. A correctly specified hierarchical/fixed-effect response model can gain substantial power, whereas matched permutation inference provides a more assumption-light calibration route when the trait model is uncertain.

## 4. Boundary

This benchmark does not show that:

- profile-LRT calibration is guaranteed under sparse or separated data;
- a common logistic slope is appropriate for arbitrary traits;
- the matched-null omnibus and logistic model estimate the same biological parameter;
- model-based methods are generally superior;
- matched permutation methods are generally more robust under MNAR missingness or measurement error.

## 5. Methods-paper state after v0.5

The synthetic comparison stack now includes:

1. naive pooled analysis;
2. equal-species matched-null inference;
3. pair-count-weighted matched-null inference;
4. species-level rho t-test;
5. species fixed-intercept logistic profile-LRT.

This closes the immediate “no model-based comparator” gap identified in v0.4.

The largest remaining gap for a standalone methods paper is now **external empirical transport**: at least one repeated, georeferenced non-flower trait dataset not used to design FCP/disttrait.

A second remaining benchmark gap is model misspecification: the logistic comparator should eventually be tested under continuous, multivariate or nonlinear trait processes where its response model is intentionally wrong.
