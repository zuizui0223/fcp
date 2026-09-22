# FCP v2 — measurement-validity programme

Date: 2026-09-23 JST  
Status: **DESIGN FROZEN BEFORE FRESH-COHORT PIXEL OPENING**

## 1. Why v2 exists

FCP v1 has already spent substantial inferential effort on observer separation, structured nulls, gate re-application, phylogenetic/sampled-span alternatives and species-specific spatial randomization. The terminal direct-highlight control changes the priority.

The remaining dominant uncertainty is now **measurement validity**, not null-model sophistication.

The current evidence is asymmetric:

- H1 shows that species-level D is reproducible under observer-disjoint resampling, but reproducibility of one image-measurement system is not proof of biological validity.
- H2 shows excess achromatic–chromatic alignment above a coarse-state-preserving null, but the direct highlight control found within-species coupling between near-clipping and white classification (OR 1.444, 95% CI 1.389–1.502).
- H2 remains supported after response-blind high-clip exclusion (W 0.503, p 0.001), so the result is not explained away by the most strongly clipped images.
- D–spatial organization replicates and survives matched flower-minus-background stress tests, but it still inherits the same upstream image/ROI measurement system.
- H3a/H3b test broad explanatory alternatives to D; they cannot establish that D itself is a biologically valid measurement.

Therefore **v2 does not add another ecological null model**. It tests whether the measurement system is invariant to technical perturbations that leave the biological flower unchanged.

## 2. Programme hierarchy

### FCP v1 — inference validity

Already complete for the present paper:

1. observer-disjoint H1 reliability;
2. recurrent H2 geometry and structured null;
3. prospective species-disjoint H2 transport;
4. replicated D–spatial organization;
5. broad phylogeny and sampled-span filters;
6. post-confirmatory direct highlight control.

No additional v1 null family is authorized by this programme.

### FCP v2 — measurement validity

Primary question:

> If the biological flower is held fixed, how much do exposure, background context and ROI localization change the measured image-level colour state, species-level D and recurrent H2 geometry?

The biological object is held fixed using **within-image technical counterfactuals**.

### FCP v3 — independent measurement-system transport

Only after v2:

> Does the same D / H2 / spatial estimand transport to an independently generated image source or calibrated field measurement system without retuning the measurement model?

v3 is not permitted to reuse same-source success as a substitute for an independent source.

## 3. Fresh-cohort structure

v2 has two metadata-selected panels.

### Panel P — paired-species fresh-image transport

Purpose: separate species identity from image-set/observer/technical effects.

- candidate species come from already used high-depth species;
- selection must be deterministic from metadata-only eligibility and a committed hash ranking;
- prior D, W, morph composition, spatial statistic and H3 outcomes must not enter selection;
- all photographs are fresh photo IDs not used in any previous FCP biological analysis;
- target: 200 species × 100 fresh rows/species;
- no species replacement after pixel opening.

This panel supports direct fresh-image transport of D and species-specific technical sensitivity.

### Panel N — novel-species measurement generalization

Purpose: ensure technical conclusions are not specific to previously successful species.

- candidate species come from the 42,111-species opportunity frame after excluding every species previously used in discovery, reserve, P500 and the third cohort;
- selection uses only metadata availability and deterministic hash ranking;
- target: 200 species × 100 rows/species;
- no replacement after pixel opening.

The complete target is therefore **400 species / 40,000 fresh images**, subject to a pre-pixel metadata availability gate. If either panel cannot meet its fixed denominator before pixel opening, v2 is underidentified and stops rather than replacing species after biological measurement.

The target of 100 rows/species preserves direct comparability with the v1 high-depth measurement design. The equal 200/200 panel split is a design choice, not a biological prevalence estimate.

## 4. Measurement-first chronology

The v2 chronology is mandatory:

```
metadata-only selection
        ↓
source-byte acquisition + SHA seal
        ↓
technical diagnostics and counterfactuals
        ↓
technical table sealed
        ↓
predefined technical-quality strata frozen
        ↓
ONLY THEN biological colour states/palettes opened
        ↓
image-level invariance
        ↓
species-level D invariance
        ↓
H2 geometry invariance
        ↓
D–spatial measurement sensitivity
```

No species-level biological outcome may be joined before the technical seal.

## 5. Technical channels measured before biological opening

Every successfully decoded image must persist the following measurement-only information.

### 5.1 Source identity and decode

- blinded measurement ID;
- source image SHA256;
- decoded width and height;
- source file bytes;
- MIME/codec;
- EXIF orientation;
- EXIF exposure fields when present, retained as technical metadata only.

Source URL, photo ID, species and coordinates are held outside technical workers.

### 5.2 Exposure

On the frozen flower mask and matched background mask:

- `clip_fraction`;
- `near_clip_fraction`;
- channel-specific clip fractions;
- linear-light luminance q01/q10/q50/q90/q99;
- dynamic-range span q99-q01;
- flower/background luminance contrast;
- fraction of mask pixels at exact black.

These are saved as continuous variables. No white/nonwhite response is available at this stage.

### 5.3 Background/context

- background effective pixels;
- background L/a/b mean, SD and quantiles;
- background coarse palette fractions;
- flower-minus-background Lab distance;
- flower/background luminance contrast;
- background clipping metrics;
- proportion of flower-mask boundary adjacent to high-luminance background.

### 5.4 ROI stability

Retain:

- detector box count/confidences;
- flower-mask effective pixels;
- image fraction covered by the mask;
- horizontal-flip mask IoU;
- horizontal-flip colour delta-E;
- deterministic prompt-jitter mask IoU distribution;
- deterministic prompt-jitter flower-colour delta-E distribution;
- mask-boundary fraction;
- number of retained flower instances.

## 6. Within-image technical counterfactuals

These are the centre of v2 because the flower is biologically unchanged.

### 6.1 Fixed-mask photometric counterfactual

Using the original admitted flower mask, convert sRGB pixels to linear light and apply exposure multipliers corresponding to:

`EV = {-1.0, -0.5, 0, +0.5, +1.0}`.

After clipping to [0,1], convert back to sRGB and run the frozen palette/classification mapping on the same flower pixels.

This isolates **colour-classifier sensitivity to exposure** from ROI changes.

No transform level may be added or removed after colour outcomes are opened.

### 6.2 Full-pipeline exposure counterfactual

For the same five EV levels, apply the photometric transformation to the full image and rerun the unchanged detector + segmenter + colour measurement.

This estimates combined exposure sensitivity of:

`detector → ROI → palette → morph`.

### 6.3 Background-neutralization counterfactual

Using the original response-blind flower mask, replace non-flower pixels with fixed mid-grey linear-light background and rerun the unchanged full localization/measurement pipeline.

Primary fixed neutral background: sRGB (128,128,128).

This tests whether background context changes the upstream localization/measurement result while the flower pixels are unchanged.

### 6.4 ROI prompt-jitter counterfactual

For every retained detector box, apply the fixed deterministic prompt set:

- x/y translation: ±5% box width/height;
- width/height scale: ±10%;
- unperturbed prompt.

No random jitter is permitted.

The same EfficientSAM weights and mask-selection rule are used for every prompt. The resulting mask IoU and flower-colour delta-E quantify ROI sensitivity.

## 7. Primary v2 estimands

v2 reports effect sizes rather than searching for the most favourable threshold.

### MV1 — image-level classification invariance

For each technical counterfactual:

- coarse-state flip probability from the unperturbed state;
- continuous palette Hellinger distance;
- flower Lab delta-E;
- transition matrix among the four coarse biological states;
- changes in classifiability.

Results are equal-image and equal-species summaries.

### MV2 — species-level D invariance

For each species and counterfactual condition:

- D under the original measurement;
- D under the counterfactual;
- absolute D change;
- rank correlation and concordance across species.

Panel P additionally compares v2 fresh-image D with the earlier frozen D for the same metadata-selected species.

### MV3 — H2 geometry invariance

The already frozen q_white direction is never refit.

Report W under:

1. original v2 measurement;
2. fixed-mask exposure counterfactuals;
3. full-pipeline exposure counterfactuals;
4. background-neutralized counterfactual;
5. ROI-stable-only subset defined before biological opening.

The purpose is measurement sensitivity, not another confirmatory H2 discovery claim.

### MV4 — spatial measurement sensitivity

For species with sufficient georeferenced support, compare the original flower-colour spatial statistic with:

- background-neutralized measurement;
- ROI-stable-only measurement;
- matched flower-minus-background response.

No new environmental or pollinator predictor is authorized in v2.

## 8. Technical seal and outcome firewall

Before colour outcome opening, the technical artifact must contain:

- one row per frozen image;
- all technical metrics;
- counterfactual availability/status;
- source image SHA;
- no species, taxon, coordinates, observer, D, morph, q_white projection, W, spatial outcome or H3 variable.

The seal must publish:

- exact row denominator;
- source-byte success/failure counts;
- technical missingness by channel;
- fixed technical-stratum definitions;
- SHA256 of the complete technical table.

Only after read-back validation of this seal may the biological join occur.

## 9. No-rescue rules

After biological outcome opening, v2 forbids:

- new exposure transforms;
- new clipping thresholds;
- alternative neutral-background values;
- alternative ROI-jitter magnitudes;
- species replacement;
- changing q_white;
- changing H2 thresholds/nulls;
- dropping a technical channel because its result is inconvenient;
- adding an ecological predictor to explain a measurement-validity failure.

A failed or unstable measurement channel is itself a result.

## 10. How H1–H3 are used

v2 does **not** rerun H1–H3 mechanically as another battery of significance tests.

- H1 motivates D transport/invariance across fresh images and technical counterfactuals.
- H2 motivates geometry invariance under exposure/background/ROI perturbation.
- D–spatial motivates matched background and ROI-stability sensitivity.
- H3a/H3b are downstream only: they may be revisited only if D measurement validity is acceptable. They are not allowed to rescue v2.

Thus the programme changes the ordering from:

`measurement → inference → more inference`

to:

`measurement → measurement validation → inference transport`.

## 11. v3 admission gate

v3 independent-source transport is opened only after v2 produces a terminal measurement-validity report.

The independent source must differ in image-generation process, not merely provide another iNaturalist species subset. Preferred designs include:

1. paired field photographs with locked exposure and a colour reference plus ordinary auto-exposure JPEGs of the same flowers;
2. a curated external field-image collection with documented acquisition and colour standards;
3. another community-science source only if its acquisition/processing pipeline is demonstrably independent.

The v3 target is **same estimands, different measurement system**. Species overlap is useful rather than disqualifying because it separates measurement-system transport from biological composition.

## 12. Current-paper firewall

The present New Phytologist manuscript is not reopened to add v2 biological results.

The current paper ends with the terminal direct-highlight result:

- exposure coupling detected;
- high-clip-exclusion H2 support retained;
- executable measurement-validity gate INDETERMINATE.

v2 is a new prospective measurement-validity programme and must be reported separately unless a future editorial request explicitly requires otherwise.
