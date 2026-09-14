# Prospective P500 H2 white measurement-control protocol

Freeze ID: `polymorphism-h2-p500-white-measurement-control-20260913`

Status: **PRE-OPENING FROZEN VALIDITY CONTROL**

Base commit: `c55fb858579074d2e3a4ce1a71b5a2d266e11011`

Parent candidate-metadata protocol: `docs/POLYMORPHISM_H2_P500_CANDIDATE_METADATA_PROTOCOL_20260913.md`

Parent preflight result: `results/polymorphism_h2_p500_candidate_metadata_20260913/result.json`

## 1. Purpose and estimand separation

The prospective P500 H2 analysis estimates the already frozen white-vs-nonwhite alignment. This measurement-control gate asks a different question: whether the prospective white state is detectably coupled to digital highlight/exposure failure within species.

This control **must not** redefine, rotate, optimize, select, or tune the H2 axis. It is a validity gate on interpretation, not an alternative H2 estimand.

The control is frozen before any P500 image bytes, flower-colour state, palette, D, H1, or H2-W outcome is opened.

## 2. Pre-opening receipt

At the base commit, the candidate-metadata gate passed with 499/500 species meeting the >=100-record criterion and zero request-error species. The following P500 fields were still unopened and are required to be `false` before first image opening:

- `image_results_requested`
- `filtered_specimens_opened`
- `image_urls_opened`
- `image_bytes_opened`
- `flower_colour_opened`
- `palette_opened`
- `D_opened`
- `H1_opened`
- `H2_W_opened`

Any mismatch fails closed.

## 3. Allowed technical measurements

Technical diagnostics are computed only from decoded RGB pixels after the frozen flower mask is available. They are deliberately separated from biological colour variables.

For 8-bit decoded RGB flower-mask pixels, let `M = max(R,G,B)` per pixel. Freeze the following diagnostics:

1. `clip_fraction = mean(M == 255)`.
2. `near_clip_fraction = mean(M >= 250)`.
3. `luminance_q99`, the 0.99 quantile of linear-light relative luminance after the same sRGB transfer function used by `fcp_pipeline.photo_first_measurement.srgb_to_lab`; luminance is `0.2126729 R_lin + 0.7151522 G_lin + 0.0721750 B_lin`.

If a source is decoded at a different integer bit depth, values are normalized to [0,1], exact clipping means the channel ceiling, and near clipping means >=250/255 of the ceiling. The decode bit depth must be recorded.

These are **digital highlight diagnostics**. They do not identify physical sensor saturation, exposure settings, or flower reflectance.

### Forbidden artifact predictors

The validity gate must not use HSV/HSL saturation, Lab chroma, palette distance, white fraction, or any derivative of the biological white-vs-nonwhite classification as a technical predictor. Those quantities are circular with the phenotype under test.

Background or mask-contamination diagnostics may be added only if they already exist as objective, response-blind outputs of the frozen segmentation pipeline. They cannot replace the primary highlight diagnostic.

## 4. Response-blind high-clipping set

The primary technical predictor is continuous `near_clip_fraction`.

Before opening the white/nonwhite response, define a response-blind high-clipping exclusion set from technical values only:

`high_clip = near_clip_fraction > max(0.01, q95(near_clip_fraction))`

where `q95` is computed across evaluable P500 images without access to any colour class or H2 result. This set is frozen before the biological response is joined to the technical table.

This rule is a sensitivity filter, not a new primary analysis cohort definition.

## 5. Primary coupling analysis

The primary artifact-coupling analysis is within species.

- Response: prospective image-level white (`1`) versus nonwhite (`0`) state from the frozen colour pipeline.
- Predictor: `near_clip_fraction`, standardized within species among species containing both response states and nonzero predictor variance.
- Nuisance structure: species fixed effects / conditional stratification by species. No between-species association may substitute for this analysis.
- Report: odds ratio (OR) per one within-species SD increase and its two-sided 95% confidence interval.

The operational negligible-coupling interval is frozen as **OR 0.80 to 1.25 per within-species SD**. This is an interpretation tolerance for the validity gate, not a biological equivalence claim.

Strict `clip_fraction` and `luminance_q99` are secondary diagnostics. They cannot turn an indeterminate or flagged primary gate into CLEAR.

## 6. H2 sensitivity analysis

Run the prospective H2 analysis exactly as preregistered on:

1. the full eligible P500 cohort; and
2. the same cohort after removing the response-blind `high_clip` images.

No H2 threshold, axis, null model, species inclusion rule, or statistic may be re-tuned for the sensitivity run.

The sensitivity run is evaluable only if it retains >=90% of the species that were evaluable in the primary prospective H2 analysis. Otherwise the measurement-control decision is INDETERMINATE.

## 7. Frozen decision rule

The only valid gate states are `CLEAR`, `FLAGGED`, and `INDETERMINATE`.

### CLEAR

Return CLEAR only if all are true:

1. the primary coupling model is estimable;
2. its entire 95% CI lies inside OR [0.80, 1.25];
3. the high-clipping sensitivity retains >=90% of primary H2-evaluable species; and
4. the frozen prospective H2 support/non-support decision is unchanged after high-clipping exclusion.

A nonsignificant p-value is never sufficient for CLEAR.

### FLAGGED

Return FLAGGED if either is true:

1. the entire primary 95% CI lies above 1.25 or below 0.80; or
2. a primary prospective H2 SUPPORT decision changes to NON-SUPPORT after the response-blind high-clipping exclusion while the sensitivity retains >=90% of primary H2-evaluable species.

A FLAGGED result means that a positive H2 result must be described as measurement-confounded rather than as prospective replication.

### INDETERMINATE

Return INDETERMINATE for every other case, including unavailable diagnostics, non-estimable within-species coupling, an equivalence CI that overlaps a boundary, insufficient sensitivity retention, missing receipts, or any chronology violation.

This is intentionally fail closed.

## 8. Chronology firewall

The required order is:

1. candidate-metadata/preflight freeze;
2. this measurement-control protocol freeze;
3. implementation/tests and machine-readable freeze receipt;
4. freeze commit SHA / PR provenance;
5. first P500 image/pixel opening and response-blind technical-table construction;
6. freeze the `high_clip` membership before joining biological colour outcomes;
7. prospective white classification and primary coupling analysis;
8. primary and high-clipping-exclusion H2 analyses;
9. measurement-control decision;
10. prospective H2 interpretation.

No step may be backfilled from a later outcome. If chronology cannot be demonstrated, return INDETERMINATE.

## 9. Claim firewall

This freeze produces **no P500 biological result**. It authorizes no statement that H2 replicated, failed, or was free of artifact. At this stage the only permissible claim is that the validity-control rules were frozen while P500 image/colour outcomes remained unopened.
