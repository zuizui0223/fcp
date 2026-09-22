# Third-cohort H2 post-confirmatory validity audit — 2026-09-22

Status: **post-confirmatory diagnostic only; does not replace or reopen the frozen prospective H2 decision.**

Canonical machine-readable audit:

- `results/polymorphism_h2_posthoc_validity_20260922/result.json`

Reproduction script:

- `scripts/analysis/audit_polymorphism_h2_posthoc_validity_20260922.py`

Third-cohort source artifact:

- workflow artifact ID: `10496492307`
- immutable biological result commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`

The sealed third-cohort artifact contains the 49,900-row measured table and H2 result products. Image pixels and flower masks were intentionally not persisted, and no background palette was stored for this cohort. Direct exposure/highlight diagnostics therefore cannot be reconstructed from the sealed result alone.

## 1. What the frozen H2 test actually confirms

The primary third-cohort statistic was:

- species: **158**
- observed `W = 0.5172457461`
- isotropic 8-dimensional expectation: **0.125**
- coarse-state-preserving structured-null median: **0.4571428150**
- observed minus structured-null median: **0.0601029311**
- observed / structured-null median: **1.131475**
- upper-tail p: **0.001**

The structured null preserves each selected species' coarse-state counts while permuting normalized nine-colour rows among species within the same coarse state, then refitting the unchanged label-free two-means construction.

Therefore the prospective result is **not** that the full difference between W = 0.517 and isotropic 0.125 was newly established by the third cohort. Much of the expected achromatic–chromatic alignment is already implied by the measured coarse-state composition and is present in the structured null.

The confirmatory increment is narrower:

> **conditional on the measured coarse-state composition, continuous within-species displacement is more aligned with the frozen white-versus-nonwhite axis than expected under the construction-preserving null.**

This is the claim that the structured-null p-value tests.

## 2. The primary W signal is strongly concentrated in white-containing species

Among the 158 primary H2 species:

- **137 / 158 = 86.7%** have white as one of their two most frequent coarse states;
- mean per-species `(u_i · q_white)^2` when white is in the top two states: **0.5622**;
- mean when white is not in the top two states: **0.2239**.

This decomposition is expected to some degree for a statistic explicitly targeted at the white-versus-nonwhite axis, but it matters for validity interpretation: any technical process that changes whether an image is classified as white can contribute directly to the measured target.

## 3. Post hoc gate-reapplied null

The frozen structured null fixes the 158 species selected after the observed continuous minor-cluster gate. It does not reapply that gate in every null world.

Because this asymmetry could in principle make the null anti-conservative, a post hoc robustness analysis started from the **185 species** that passed the coarse second-state gate and, in every permuted world, refitted two-means and reapplied the continuous minor-cluster threshold.

Using 299 null worlds with the frozen primary seed:

- retained species per world: **166–183**, median **176**;
- null median W: **0.447495**
- null mean W: **0.447583**
- 95% interval: **0.428986–0.469580**
- null values >= observed W: **0 / 299**
- plus-one upper-tail p: **0.00333**

The gate-reapplied null is centered below, not above, the frozen structured null (0.4475 versus 0.4571). Thus this particular gate asymmetry does **not** explain the prospective support; if anything, the frozen fixed-species null is conservative relative to this post hoc alternative.

This diagnostic is not promoted into a new confirmatory test. It is reported only as robustness of the frozen null construction.

## 4. Exposure/background coupling remains unresolved

A direct digital-highlight control had previously been specified for the P500 programme before image opening. Its non-circular technical diagnostics were:

- `clip_fraction`
- `near_clip_fraction`
- `luminance_q99`

with white/nonwhite coupling assessed within species. That protocol was not executed for the successful third cohort.

Because the third-cohort workflow intentionally destroyed image pixels and masks after sealing measurement results, these highlight metrics cannot now be reconstructed without reacquiring the images and rerunning the frozen segmentation path.

A separate post hoc background-bearing sample provides an indirect warning signal. Its summary is recorded in:

- `docs/supporting/POLYMORPHISM_H2_BACKGROUND_WHITE_PROXY_20260922.json`

Across 461 species:

- median background white fraction for white-classified flower images: **0.0692**
- median for nonwhite-classified images: **0.0553**
- **66.4%** of species showed the same direction
- paired Wilcoxon p = **4.6 × 10^-12**

This is consistent with image/background or exposure coupling, but it is not a direct exposure test. It cannot distinguish among:

- overexposure or digital highlight clipping;
- genuine ecological differences in background brightness or habitat openness;
- framing differences;
- flower-mask/background leakage;
- other image-acquisition processes.

The row-level source for this separate proxy summary is not committed in this repository, and it is not the third prospective cohort. It therefore serves as a limitation/diagnostic, not as a correction factor or replacement analysis.

Most importantly, the frozen structured null cannot diagnose a technical process that creates the measured coarse white state itself, because coarse state is explicitly conditioned on by the null.

Accordingly, the paper-level interpretation ceiling is:

> **the third cohort supports excess continuous alignment with the frozen achromatic–chromatic axis conditional on the measured coarse-state composition; exposure-specific validity remains unresolved.**

## 5. disttrait is not a bitwise reproducer of the frozen H2 pipeline

`disttrait` was developed after the biological analysis as a generic implementation of the broader distributional-trait framework. A real-data audit on the same third-cohort rows shows small numerical non-equivalence:

- frozen biological pipeline: **W = 0.5172457461**
- `disttrait.two_mode_axis` specieswise route: **W = 0.5183899565**
- `disttrait.structured_alignment_null` observed route: **W = 0.5157982982**

The differences arise from nearly tied farthest-point initializations. Re-normalization and floating-point reduction order can change which almost-equidistant row is selected as a deterministic seed, which can lead two-means to a different local solution.

Sensitive species in the specieswise route include:

- *Cirsium vulgare*: absolute axis cosine to frozen solution **0.6333**
- *Vicia benghalensis*: **0.9830**

In the structured-alignment observed path:

- *Geranium reuteri*: **0.6867**
- *Vicia benghalensis*: **0.9830**

These differences do not alter the qualitative H2 conclusion, but they mean that `disttrait` must not be described as an exact numerical reproducer of the paper's frozen H2 values.

The manuscript should therefore report W at biologically meaningful precision (three decimals in prose) and identify the study-specific frozen pipeline as the authoritative numerical implementation.

## 6. Claim consequences

The frozen machine-readable verdict remains:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

because this post-confirmatory audit does not modify the prespecified axis, cohort, thresholds, statistic or null.

However, the prose claim is narrowed in two ways:

1. **increment, not total alignment** — the supported effect is the excess above the coarse-state-preserving null, not the whole distance from isotropy;
2. **measurement validity boundary** — direct exposure/highlight coupling was not tested in the successful third cohort, and an indirect background-white proxy motivates explicit caution.

No exposure, pigment, adaptive, evolutionary-direction or independent-source claim is licensed by this audit.
