# FCP v2 biological-pass protocol

Date: 2026-09-23 JST  
Status: **PREP ONLY — MUST NOT EXECUTE BEFORE RESPONSE-BLIND TECHNICAL SEAL**

## 1. Purpose

FCP v2 biological Pass B opens flower-colour outcomes only after the complete response-blind technical table from Pass T has been sealed, hashed and read back.

The question is measurement transport:

> Holding the biological flower fixed, how much do exposure, background context and ROI localization change the measured colour state, species-level D and alignment with the already frozen q_white direction?

This pass adds **no new ecological null model, no new environmental predictor and no new colour axis**.

## 2. Mandatory prerequisites

Pass B is unauthorized unless all are true:

1. the terminal fresh-metadata denominator remains exactly 400 species × 100 rows = 40,000 images;
2. Panel P and Panel N each remain exactly 200 species × 100 rows;
3. the heavy-counterfactual subset remains exactly 20 metadata-hash-selected rows/species = 8,000 rows;
4. Pass T has produced exactly one technical row per frozen measurement ID;
5. the Pass-T technical artifact has been sealed and its SHA256 committed;
6. Pass T records `biological_outcomes_opened = false`, `species_opened = false`, `morph_opened = false` and `q_white_or_W_opened = false`;
7. the frozen technical strata have already been serialized;
8. a separate Pass-B execution authorization pins the exact technical artifact ID/digest.

Failure of any prerequisite stops before Pass-B image reacquisition.

## 3. Source-byte identity

Pass B reacquires the same frozen image URL for each measurement ID.

For each successfully reacquired image:

- calculate Pass-B SHA256;
- compare with the Pass-T `source_image_sha256`;
- any mismatch is `source_byte_drift`;
- source-drifted rows are not replaced and do not enter paired same-image counterfactual inference.

Pass B must report the exact source-drift denominator before species-level results.

## 4. Frozen image measurement

The biological classifier is imported unchanged from source commit:

`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`

Required source modules:

- `fcp_pipeline/flower_roi_v4.py`
- `fcp_pipeline/flower_roi_v4_runtime.py`
- `fcp_pipeline/photo_first_measurement.py`
- `fcp_pipeline/photo_first_atlas_v2.py`
- `fcp_pipeline/photo_first_atlas.py`
- `fcp_pipeline/shared_transition_surface.py`

Frozen colour contract:

- `docs/supporting/random_photo_first_measurement_contract_v1.json`

The biological palette remains:

white, yellow, orange, red, pink, magenta, purple, blue, bronze.

Green, brown and black remain nuisance anchors excluded from the biological denominator.

The four coarse states remain:

- white;
- yellow_orange;
- red_pink;
- blue_purple.

Classification thresholds remain:

- minimum dominant coarse-state fraction = 0.50;
- minimum margin over the second group = 0.10;
- otherwise = `mixed_uncertain`.

No palette anchor, threshold or nuisance rule may be fit to v2 outcomes.

## 5. Biological measurements

### 5.1 Baseline — all 40,000 rows

For every byte-matched reacquired image:

1. run frozen ROI-v4;
2. if ROI/flip gate is admitted, classify the frozen flower-mask pixels with the frozen palette;
3. retain the complete nine-dimensional biological palette fractions;
4. otherwise retain `mixed_uncertain` as structural measurement missingness.

### 5.2 Fixed-mask exposure counterfactual — all 40,000 rows

Use the **baseline flower mask unchanged** and apply linear-light exposure shifts:

`EV = {-1.0, -0.5, 0, +0.5, +1.0}`.

For each EV:

- do not rerun detector or segmenter;
- classify transformed RGB values only inside the unchanged baseline flower mask;
- retain morph, classifiability and nine-dimensional palette fractions.

This is the primary all-row exposure-sensitivity estimand.

### 5.3 Heavy subset full-pipeline exposure — 8,000 rows

For metadata-frozen heavy rows only:

`EV = {-1.0, -0.5, +0.5, +1.0}`.

For each EV:

- transform the complete image in linear light;
- rerun the unchanged detector/segmenter;
- apply the unchanged ROI/flip gate;
- classify the resulting admitted flower mask.

### 5.4 Heavy subset background neutralization — 8,000 rows

Using the baseline response-blind flower mask:

- replace all non-flower pixels by sRGB (128,128,128);
- rerun detector/segmenter;
- apply the same ROI/flip gate;
- classify the resulting admitted flower mask.

### 5.5 Heavy subset deterministic ROI-jitter — 8,000 rows

Use the already frozen seven-prompt set:

- unperturbed detector box;
- x translation ±5%;
- y translation ±5%;
- scale ±10%.

Detector boxes are held fixed; EfficientSAM is rerun under each prompt.

Each jitter variant returns:

- classifiability;
- morph;
- nine-dimensional palette fractions.

No additional jitter magnitude or random jitter is allowed after biological opening.

## 6. Pass-B worker firewall

Partition workers may receive only:

- `measurement_id`;
- blind image filename;
- photo licence;
- heavy-counterfactual flag;
- expected Pass-T source-image SHA256.

Workers may not receive:

- species;
- inat taxon ID;
- photo/observation ID;
- coordinates;
- observer identity;
- prior D;
- prior H2/W;
- spatial outcome;
- H3 variables.

Species/coordinates open only after all 256 biological partitions have one terminal row per frozen measurement ID and source identity has been validated.

## 7. Primary v2 reporting estimands

No omnibus “v2 passed/failed” label is required. Component effects are reported directly.

### MV1 — image-level invariance

For each counterfactual, report both:

1. **complete-state disagreement**, where `mixed_uncertain` is an explicit fifth measurement state;
2. **morph-flip probability among rows classifiable in both baseline and counterfactual**.

Also report:

- baseline-to-counterfactual transition matrix;
- change in classifiability;
- Hellinger distance between nine-dimensional palette fractions;
- equal-image and equal-species summaries.

The primary small-perturbation exposure contrasts are EV -0.5 and +0.5.

EV ±1.0 are stress contrasts and are not substituted if ±0.5 is inconvenient.

### MV2 — species-level D invariance

For each species and each all-row condition:

`D = 1 - sum_k p_k^2`

over the same four biological states, requiring at least 40 classifiable rows under that condition.

Report:

- number of common evaluable species;
- baseline D;
- counterfactual D;
- signed and absolute D difference;
- Spearman rank correlation;
- Lin concordance correlation coefficient;
- calibration intercept and slope.

Panel P additionally compares **fresh baseline v2 D** with the already frozen earlier D for the same species. That comparison is a transport diagnostic, not a redefinition of H1.

### MV3 — q_white geometry invariance

The axis remains fixed:

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`.

The v2 measurement question is the stability of observed W under technical perturbation.

For each all-row condition report:

- coarse-state gate species;
- continuous-vector species;
- observed W;
- paired per-species squared-projection contribution;
- paired ΔW from baseline among common vector species.

No new structured null family is introduced in v2.

The original frozen H2 null remains provenance only; v2 asks sensitivity of the already defined estimand to measurement perturbation.

Heavy-subset full-pipeline/background/ROI W-like summaries are explicitly labelled **20-row/species measurement sensitivities**, not replacements for the all-row H2 estimand.

### MV4 — spatial measurement sensitivity

Only after MV1–MV3 are terminalized:

- calculate the existing species-specific spatial colour-organization statistic on baseline v2 measurements;
- compare it with the same statistic after response-blind technical-stratum exclusions;
- retain flower-minus-background sensitivity where available.

No new climate, pollinator, habitat or range predictor enters v2.

## 8. Technical-stratum analyses

The strata frozen in Pass T are used exactly as serialized:

- high flower near-clipping;
- high background near-clipping;
- low flower/background Lab separation;
- fixed ROI instability;
- heavy prompt instability.

These are not optimized against D/W.

Report baseline biological estimands both:

- on the complete measurement-evaluable set; and
- after the pre-frozen technical-stratum exclusion relevant to that estimand.

No new technical threshold is selected after outcome opening.

## 9. Interpretation

v2 distinguishes three possibilities without forcing one global verdict:

1. **measurement-stable** — biological estimands change little under technical perturbations;
2. **measurement-sensitive but structurally persistent** — image-level states shift, but species-level D/q_white geometry remain similar;
3. **measurement-sensitive and structurally unstable** — technical perturbations materially change species-level D or q_white geometry.

These descriptions must be supported by reported effect sizes; they are not post hoc significance labels.

## 10. Hard nonclaims

v2 must not claim:

- independent-source replication;
- calibrated biological reflectance;
- pigment chemistry;
- direction of evolutionary transitions;
- adaptive causation;
- that a nonsignificant technical association proves absence of bias;
- that technical sensitivity invalidates every biological observation;
- that one successful counterfactual rescues a failed one.

## 11. Relationship to FCP v3

v2 remains same-source measurement validation.

FCP v3 is the independent measurement-system test and remains separately gated. A different iNaturalist subset does not qualify as v3.
