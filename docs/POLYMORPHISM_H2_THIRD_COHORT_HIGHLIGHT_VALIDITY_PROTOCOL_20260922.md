# Third-cohort digital-highlight validity control — one-shot post-confirmatory execution

Date frozen: 2026-09-22 JST  
Status: **FROZEN BEFORE THIRD-COHORT IMAGE REACQUISITION FOR THIS CONTROL**

## Role

This control is a **post-confirmatory measurement-validity analysis**. It does not replace, rerun, rescue, invalidate by fiat, or redefine the frozen prospective H2 decision `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

The diagnostic procedure is inherited from the pre-existing P500 white measurement-control protocol frozen on 2026-09-13. Its application to the successful third cohort is being frozen here before any third-cohort image is reacquired for this control.

Whatever terminal state is obtained — `CLEAR`, `FLAGGED`, or `INDETERMINATE` — will be reported without changing thresholds, predictors, subsets, the white axis, or the decision rule.

## Fixed population and source identity

The control uses the complete authorized third-cohort denominator, not the 158 H2-positive species subset.

- authorized species: **499**
- authorized rows: **49,900**
- target rows/species: **100**
- authorized metadata SHA256: `13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290`
- source metadata workflow artifact: **10451429156**
- original successful biological workflow artifact: **10496492307**
- immutable biological result commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
- exact inherited measurement source commit: `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`

No species, observation, photo, or failed row may be replaced. No new iNaturalist query is permitted. The exact frozen `photo_id` / `photo_url_large` rows are reacquired once.

## Response-blind technical stage

The technical stage must complete before the frozen biological response table is opened.

For each of the 49,900 authorized rows:

1. reacquire the exact frozen `photo_url_large`;
2. run the exact ROI-v4 / EfficientSAM measurement infrastructure from source commit `9fae6cc...`;
3. use the returned original-resolution flower mask;
4. calculate only the following technical diagnostics from flower-mask RGB pixels:
   - `clip_fraction = mean(max(R,G,B) == 255)`;
   - `near_clip_fraction = mean(max(R,G,B) >= 250)`;
   - `luminance_q99`, the 0.99 quantile of linear-light relative luminance under the sRGB transfer function.
5. seal only technical outputs and blinded row IDs; do not calculate colour palettes, morphs, D, H2, or use species identity in the measurement workers;
6. destroy downloaded image pixels, flower masks and model files after each partition is sealed.

The response-blind high-clipping rule is inherited unchanged:

`high_clip = near_clip_fraction > max(0.01, q95(near_clip_fraction))`

The threshold and the complete high-clip ID set must be serialized and uploaded before the biological outcome artifact is downloaded in the outcome-join job.

Rows lacking an admitted flower mask have unavailable highlight diagnostics and remain outside the high-clip exclusion set. Their counts are reported; they are never replaced.

## Biological join and image-identity audit

Only after the technical table and high-clip membership are sealed may the workflow open the immutable third-cohort biological artifact.

The join uses the frozen photo IDs through a sealed ID map. Reacquired image SHA256 values are compared with the original biological measurement `image_sha256` values when both are available.

Any nonzero count of successfully reacquired rows whose bytes differ from the frozen biological image SHA is reported as source drift and forces the final control state to `INDETERMINATE`. Network/acquisition failures are reported separately and are not replaced.

## Primary coupling analysis

Response: original frozen image-level coarse state, coded white = 1 and the other three biological states = 0. Only originally classifiable rows with an available `near_clip_fraction` enter the model.

Predictor: `near_clip_fraction`, standardized **within species**.

Species eligibility for the coupling model:

- both white and nonwhite original frozen responses are present;
- at least two model rows;
- within-species predictor variance is nonzero.

Primary model: conditional logistic regression stratified by species, with one common coefficient for within-species standardized `near_clip_fraction`.

Report:

- model rows and species;
- coefficient;
- odds ratio per one within-species SD;
- two-sided 95% confidence interval.

Frozen negligible-coupling interval: **OR 0.80 to 1.25**.

`clip_fraction` and `luminance_q99` are descriptive secondary diagnostics only and cannot rescue the primary decision.

## H2 high-clip sensitivity

Starting from the immutable original third-cohort measured table:

1. remove rows whose blinded technical IDs belong to the response-blind high-clip set;
2. do not replace rows;
3. apply the same original classifiability rule, n_classifiable >=40 support rule, 0.10 primary coarse-state threshold, 0.10 continuous minor-cluster threshold, frozen q_white, W statistic, and 999 structured-null algorithm;
4. use the same primary null seed `20260915`;
5. do not refit, rotate, optimize or rename the white axis.

The frozen primary H2 reference is:

- primary vector species: **158**
- W = **0.5172457461053418**
- structured-null p = **0.001**

Sensitivity species retention is:

`N_sensitivity_vectors / 158`

and must be >= **0.90** (at least 143 primary-tier vectors) for the sensitivity to participate in `CLEAR` or the H2-flip branch of `FLAGGED`.

## Frozen terminal decision

### CLEAR

Return `CLEAR` only if all are true:

1. the conditional logistic coupling model is estimable;
2. the entire 95% OR interval lies within **[0.80, 1.25]**;
3. source-drift count is zero;
4. sensitivity primary-vector retention is >=0.90;
5. the frozen primary H2 support decision remains supported after high-clip exclusion.

### FLAGGED

Return `FLAGGED` if either is true:

1. the entire coupling-model 95% OR interval lies above 1.25 or below 0.80, with zero source drift; or
2. the original supported H2 decision becomes non-support after high-clip exclusion while vector retention is >=0.90 and source drift is zero.

A FLAGGED state means the paper must describe the white-axis result as measurement-confounded.

### INDETERMINATE

Return `INDETERMINATE` otherwise, including:

- non-estimable coupling model;
- confidence interval overlapping a frozen equivalence boundary;
- source-byte drift;
- H2 vector retention <0.90;
- incomplete terminal technical census or chronology failure.

A nonsignificant p-value alone is never sufficient for CLEAR.

## One-shot / no-expansion rule

This execution authorizes exactly one reacquisition-and-control run over the frozen 49,900-row denominator.

After the terminal result exists:

- do not change the clipping threshold;
- do not add saturation, Lab chroma, palette distance or white fraction as technical predictors;
- do not switch to the 158 H2 species before defining the high-clip set;
- do not alter the coupling model, OR interval, sensitivity threshold, q_white, H2 null or H2 seed;
- do not add a second exposure control because the first result is inconvenient.

The terminal result will be incorporated into the manuscript regardless of state.
