# disttrait methods-paper figure and table plan — 2026-09-19

Status: source-locked display plan for the v0.12 methods manuscript.

Manuscript:
`docs/DISTTRAIT_METHODS_MANUSCRIPT_V0_1_20260919.md`

Package:
`disttrait 0.12.0`

All display values must be read from frozen repository receipts or frozen fixture files. Figure generation must not rerun biological inference or substitute newly optimized thresholds.

## Figure 1 — pooling answers the wrong level of organization

Purpose: make the between-species-confounding problem immediately visible.

Source:
`results/disttrait_species_conditioning_benchmark_v0_2_20260918/result.json`

Display:
- null worlds: naive pooled false-positive fraction versus species-conditioned matched-null false-positive fraction;
- signal worlds: naive pooled detection versus species-conditioned detection.

Frozen values:
- null naive = 1.00;
- null conditioned = 0.00;
- signal naive = 1.00;
- signal conditioned = 1.00.

Interpretation:
species-conditioning, rather than one particular aggregation statistic, is the first essential inferential step.

## Figure 2 — calibration, power and observation-process boundary

Purpose: distinguish finite-sample power/calibration from observation-process validity.

Sources:
- `results/disttrait_performance_surface_v0_3_20260919/cells.csv`
- `results/disttrait_mnar_observation_v0_10_20260919/cells.csv`

Display A:
a single heatmap of conditioned detection fraction across effect size, imbalance ratio and MCAR fraction.

Display B:
observed matched-null rejection fraction under MCAR, trait-only, position-only and joint trait-by-position selection, at both selection strengths.

Required annotations:
- null performance-surface matched-null maximum = 0.025;
- joint MNAR rejection = 0.80–1.00;
- non-joint MNAR rejection <= 0.05.

Interpretation:
matched randomization can be calibrated conditional on the observed rows but cannot recover a latent biological null after joint trait-by-location selection changes the observed sample.

## Figure 3 — method performance follows the estimand

Purpose: show that model-based efficiency and direction-invariant organization answer different questions.

Sources:
- `results/disttrait_model_comparator_surface_v0_5_20260919/cells.csv`
- `results/disttrait_direction_heterogeneity_v0_7_20260919/cells.csv`
- `results/disttrait_nonlinear_curvature_v0_11_20260919/cells.csv`
- `results/disttrait_multivariate_orientation_v0_12_20260919/cells.csv`

Single-axis scenario plot at the frozen weak/moderate effect levels:
1. binary logistic, shared response;
2. continuous linear, shared direction;
3. continuous linear, 50% reversal;
4. quadratic, shared curvature;
5. quadratic, 50% curvature reversal;
6. multivariate, shared orientation;
7. multivariate, full-circle orientation spread.

For each scenario, display detection fraction for:
- the response-model comparator appropriate to that scenario;
- the matched-null distance/dissimilarity estimand.

Use the no-missingness cell for the main display and report the missingness range in the caption.

Interpretation:
correctly specified response models are highly efficient when species share the modeled response, whereas distance/dissimilarity estimands retain organization when species-specific directions or shapes cancel in a shared signed parameter.

## Figure 4 — flexible species slopes require calibration

Purpose: turn the failed random-effects preflight into a positive methodological result.

Source:
`results/disttrait_random_slope_meta_v0_8_20260919/cells.csv`

Display:
- maximum null rejection for analytic meta mean, analytic heterogeneity, analytic omnibus;
- maximum null rejection for permutation-calibrated mean, heterogeneity and combined omnibus;
- effect 0.4 / 50% reversal detection for common slope, matched-null organization, calibrated slope heterogeneity and calibrated slope-meta omnibus.

Key values:
- analytic heterogeneity / omnibus worst null rejection = 0.30;
- calibrated combined omnibus worst null rejection = 0.05;
- weak-effect 50%-reversal calibrated heterogeneity = 0.75–1.00;
- weak-effect 50%-reversal calibrated omnibus = 0.525–1.00.

Interpretation:
model flexibility and test calibration are separate problems.

## Figure 5 — external empirical transport can be positive or non-supporting

Purpose: demonstrate executable transport rather than a universal ecological result.

Sources:
- `results/disttrait_sf_street_tree_empirical_v0_6_20260919/result.json`
- `results/disttrait_gammarus_empirical_v0_9_20260919/result.json`

Display:
- San Francisco street-tree equal-taxon mean spatial rho = 0.09219, matched-null p = 0.01;
- ShareTrait Gammarus log-metabolic-rate rho = -0.00683, p = 0.68;
- raw-rate Gammarus sensitivity = 0.01093, p = 0.198.

Do not visually imply that the two applications estimate an identical population-level parameter:
- tree value is an equal-taxon aggregate across 20 taxon labels;
- Gammarus value is one-species / three-population organization.

Interpretation:
the same inference layer transports across unrelated systems and can return either support or non-support.

## Table 1 — method / estimand map

Required columns:
- method;
- species-conditioned?;
- trait representation;
- direction-sensitive?;
- model/distribution assumption;
- calibration;
- primary estimand;
- benchmark role;
- main limitation.

Rows:
1. naive pooled pairwise analysis;
2. equal-species matched-null distance/dissimilarity;
3. pair-count-weighted matched-null;
4. species-rho one-sample test;
5. fixed-intercept logistic common slope;
6. fixed-intercept continuous common slope;
7. common quadratic curvature;
8. common multivariate response vector;
9. species-specific signed-slope random-effects summary;
10. permutation-calibrated slope/meta omnibus.

## Generation contract

Canonical generator:
`scripts/analysis/make_disttrait_methods_figures.py`

Canonical output directory:
`docs/figures/disttrait_methods_20260919/`

Expected outputs:
- `disttrait_figure1_pooling.png` and PDF;
- `disttrait_figure2_calibration_boundary.png` and PDF;
- `disttrait_figure3_estimand_alignment.png` and PDF;
- `disttrait_figure4_meta_calibration.png` and PDF;
- `disttrait_figure5_empirical_transport.png` and PDF;
- `disttrait_table1_method_estimand_map.csv`;
- `manifest.json`.

The generator is reporting-only. It must not alter frozen receipts, rerun simulation worlds or change any result contract.
