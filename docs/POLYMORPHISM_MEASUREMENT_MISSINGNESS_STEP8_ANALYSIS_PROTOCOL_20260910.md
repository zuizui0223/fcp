# Polymorphism measurement missingness — Step 8 analysis protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-measurement-missingness-step8`
Preflight parent result: `results/polymorphism_measurement_missingness_step8_preflight_20260910/result.json`

## Preflight result used to freeze the decomposition

The reserve ledger contains 50,000 rows = 500 species x exactly 100 fixed photos/species. Its measurement-status vocabulary is exactly:

- `classified_four_state_morph`;
- `not_evaluable_roi_or_flip_gate`;
- `not_evaluable_ambiguous_palette_composition`;
- `not_evaluable_no_biological_palette_mass`.

The preflight was outcome blind (`D_association_computed = false`).

## Frozen measurement-process classes

For every species in discovery and reserve, using all 100 fixed acquired photos:

1. **technical_failure_rate** = fraction with `measurement_status == not_evaluable_roi_or_flip_gate`;
2. **ambiguous_palette_rate** = fraction with `measurement_status == not_evaluable_ambiguous_palette_composition`;
3. **no_biological_palette_rate** = fraction with `measurement_status == not_evaluable_no_biological_palette_mass`;
4. **classifiable_rate** = fraction with `measurement_status == classified_four_state_morph`.

These four classes must sum to 1 for each species.

`technical_failure_rate` is the only measurement-process control admitted as a clearly operational failure. `ambiguous_palette_rate` and `no_biological_palette_rate` are retained as diagnostics and are not silently treated as technical failures.

## Frozen inputs

Discovery raw measurement ledger:
`data/derived/global_monte_carlo_measured_photos_v1.csv`

Reserve raw measurement ledger:
`data/derived/rgfca_reserve_replication_measured_photos_v1.csv`

Discovery D/species attributes:
`results/polymorphism_species_attributes_step4_association_20260910/species_analysis_table.csv`

Reserve D/span table:
`results/polymorphism_spatial_span_adjusted_step6_20260910/reserve_species_span_adjusted_input.csv`

Discovery observed spatial response:
`results/polymorphism_spatial_subset_step5_20260909/species_membership_and_observed_rho.csv`

Reserve observed spatial responses:
`results/polymorphism_spatial_reserve_step5b_20260909/reserve_species_membership_and_observed_metrics.csv`

Discovery frozen spatial null run: `34088925008`.
Reserve frozen spatial null run: `34178957447`.

## Analysis A — measurement decomposition diagnostics

In discovery and reserve separately report:

- Spearman rho of D with technical_failure_rate;
- rho of D with ambiguous_palette_rate;
- rho of D with no_biological_palette_rate;
- rho of D with classifiable_rate;
- genus clustering of technical_failure_rate;
- genus clustering of ambiguous_palette_rate.

These are measurement-process diagnostics. They are not new biological trait hypotheses.

## Analysis B — genus clustering after pure technical-failure control

Reuse the exact equal-genus weighted-pair statistic from Step 4 / Step 7.

Primary Step-8 genus robustness statistic:

1. rank-transform D, sampled log-span, and technical_failure_rate;
2. regress ranked D on intercept + ranked sampled log-span + ranked technical_failure_rate;
3. apply the unchanged genus clustering statistic to the residual;
4. 20,000 lower-tail permutations among species in repeated genera.

Run in discovery and reserve.

Predeclared sensitivity: replace D with `D_unbiased = n/(n-1) * D` before the same rank residualization.

For comparison only, retain the already-frozen Step-7b span-only and `n_classifiable`-only results; do not retune those analyses.

### Genus decision

If reserve genus clustering remains positive with P<0.05 after `span + technical_failure_rate`, then clearly operational ROI/flip failure is insufficient to explain the replicated genus signal.

If it fails, the genus signal remains measurement-process dependent and should not be promoted as opportunity-robust.

Regardless of outcome, ambiguity-related missingness remains unresolved because ambiguous palette rows are deliberately not declared technical failures.

## Analysis C — D-spatial association after pure technical-failure control

For observed data and every one of the original 999 within-species spatial null realizations, compute:

`partial Spearman(D, species spatial response | sampled log-span, technical_failure_rate)`.

Residualize rank(D) and rank(spatial response) separately on intercept + rank(sampled span) + rank(technical_failure_rate). The observed statistic is their residual Pearson correlation.

Apply to:

1. discovery primary spatial rho;
2. reserve primary spatial rho;
3. reserve matched flower-minus-background differential.

Use the frozen 999 within-species colour randomizations. No spatial permutations are regenerated. Directional P value:

`(1 + number(null partial >= observed partial)) / 1000`.

Sensitivity: repeat with D_unbiased.

### Spatial decision

The pure technical-failure alternative is insufficient to explain the D-spatial relationship only when the adjusted effect remains positive and exceeds the geometry-preserving spatial null.

The reserve matched-background differential is the preferred flower-specific diagnostic because it subtracts matched background spatial structure.

## Claim boundary

Step 8 can separate clearly operational ROI/flip failures from other missingness. It cannot prove that ambiguous palette composition is biological polymorphism, nor can it recover excluded morph states. Successful robustness therefore permits only:

> the reported taxonomic/spatial relationships are not explained by sampled geographic span and the rate of clearly technical ROI/flip failures alone.

It does not establish phylogenetic signal, genetic determination, adaptation, selection, causal range-size effects, population-genetic differentiation, or a shared global colour boundary.
