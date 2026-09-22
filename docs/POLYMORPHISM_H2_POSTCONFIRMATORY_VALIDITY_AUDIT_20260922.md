# Post-confirmatory H2 white-axis validity audit — 2026-09-22

Status: **post-confirmatory validity/robustness audit; frozen prospective H2 decision unchanged.**

Authoritative frozen H2 result:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`
- verdict: `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

This audit asks what that verdict establishes, and what measurement-validity risk remains. It does not alter the frozen axis, cohort, thresholds, statistic, null, or decision rule.

## 1. The confirmed quantity is an increment over the construction baseline

For the primary 0.10 tier:

| Reference | W |
|---|---:|
| Isotropic 8-dimensional reference | 0.125 |
| Coarse-morph-composition-preserving structured null, median | 0.457143 |
| Observed | 0.517246 |

The prospective upper-tail probability is `p = 0.001`, and observed W is 1.1315 times the structured-null median.

The inferential reference is the **structured null**, not isotropy. The structured null preserves each selected species' coarse-morph row counts and the global mapping from coarse states to continuous nine-colour rows. It therefore already contains most of the white-versus-nonwhite alignment induced by coarse morph composition.

The prospective result supports the narrower statement:

> **Continuous within-species colour displacement is more strongly aligned with the frozen achromatic–chromatic axis than expected after preserving coarse morph composition and the corresponding global palette construction.**

It should not be described as though the test discovered the white/nonwhite direction from an isotropic starting point.

## 2. White involvement is concentrated among the 158 primary H2 species

Recalculation from the frozen primary delta-vector table gives:

- H2 species: **158**
- white appears as the primary or secondary coarse morph: **137 / 158 (86.7%)**
- mean species-level `(u_i · q_white)^2` when white is involved: **0.5622**
- mean when neither of the two leading coarse morphs is white: **0.2239**
- overall W: **0.5172457461**

This is descriptive decomposition, not a separate significance test. It shows that the positive H2 statistic is heavily concentrated in species whose coarse morph composition contains white.

## 3. Reapplying the continuous-minor gate in every null world is conservative relative to the frozen null

The frozen prospective null fixes the 158 species selected by the observed continuous-minor gate. Because this selection is data dependent, we ran a post-confirmatory robustness analysis beginning from the **185** species that passed the coarse second-morph >=0.10 gate and reapplying the continuous minor-cluster >=0.10 gate independently in each null world.

Using 299 null worlds and the frozen primary seed:

- median null W: **0.447495**
- mean null W: **0.447583**
- 95% empirical interval: **0.428986–0.469580**
- retained species per null world: **166–183**, median **176**
- observed W: **0.517246**
- upper-tail p: **0.00333**, the minimum attainable with 299 null worlds

The original frozen null median was **0.457143**, which is higher than the re-gated null median. The frozen conditional null is therefore not anti-conservative with respect to this gate asymmetry; if anything, it is slightly more conservative in this diagnostic.

This analysis is post-confirmatory and does not replace the 999-world frozen prospective test.

## 4. Background-white proxy identifies an unresolved measurement-validity risk

The third-cohort measured table does not contain a matched background palette. A separate pre-existing reserve resource does. Using the frozen reserve row table at commit `f14186590c11ac24c95e1985077908b732132e96`, we computed for each classifiable image:

`background_white_fraction = background_palette_count_white / background_effective_pixels`.

For each species containing both white and nonwhite classifiable images, the median background-white fraction was calculated separately for white and nonwhite images. Across **461** paired species:

- median of species-level white-image medians: **0.069210**
- median of species-level nonwhite-image medians: **0.055301**
- species with white > nonwhite background-white fraction: **66.38%**
- median within-species difference: **0.017029**
- paired Wilcoxon two-sided p: **4.55 × 10^-12**

This association is compatible with brightness/exposure coupling, but it is **not a direct exposure test**. At least three alternatives remain:

1. white flowers may occur or be photographed against biologically brighter backgrounds;
2. image framing/background choice may covary with flower colour;
3. ROI or mask leakage may couple flower and background white content.

Accordingly, the proxy establishes a measurement-validity concern that must be disclosed; it does not establish that the H2 result is an exposure artifact.

## 5. The previously frozen digital-highlight protocol was never executed

A P500-specific protocol was frozen on 2026-09-13 with non-circular digital-highlight diagnostics:

- `clip_fraction`
- `near_clip_fraction`
- `luminance_q99`

and with HSV/HSL saturation, Lab chroma, palette distance, and white fraction explicitly forbidden as technical predictors.

That protocol was **not executed**. It also belongs to the earlier P500 chronology and cannot be retrospectively treated as a prospective validity gate for the successful third cohort.

The third-cohort artifact retains photo IDs and the historical authorized metadata retains image URLs, so the exact images could be refetched and the same non-circular highlight diagnostics computed. Any such analysis now would be a **post hoc measurement-validity analysis**, not part of the frozen prospective H2 decision.

## 6. Structured-null scope

The structured null protects against a construction artifact in which coarse morph composition alone generates the apparent continuous-axis concentration.

It does **not** protect against a measurement process that helps create the coarse white classification itself. If brightness/exposure causes both the coarse white state and the continuous palette displacement, preserving coarse morph state in the null also preserves the confounded state.

Therefore the correct interpretation is:

- **supported:** excess continuous white-axis alignment beyond the coarse-composition construction baseline;
- **unresolved:** whether image brightness/exposure, biological background differences, segmentation leakage, or true floral colour biology contributes to the white state itself.

## 7. disttrait numerical boundary

`disttrait 0.12.0` is a later generalization of the analytical architecture, not the engine that generated the manuscript estimates.

An actual-data audit found small differences among the frozen H2 implementation and current disttrait paths at floating-point/two-means tie boundaries. The inferential conclusion was unchanged, but the implementations are not guaranteed to be bitwise or numerically identical for every real-data configuration.

Accordingly:

- manuscript numbers remain those from the study-specific frozen implementation;
- disttrait should not be cited as a numerical replay of the frozen H2 result;
- W is reported to three decimals in prose, with exact values retained in machine-readable results.

## 8. Reproducibility

Machine-readable audit:

- `results/polymorphism_h2_postconfirmatory_validity_audit_20260922/result.json`

Recomputation scripts:

- `scripts/analysis/audit_polymorphism_h2_posthoc_validity_20260922.py`
- `scripts/analysis/audit_polymorphism_background_white_proxy_20260922.py`

The background proxy was independently regenerated in GitHub Actions run `35673190275`, job `106574032963`.

## Claim ceiling

This audit does not overturn the prospective H2 decision. It narrows its interpretation and adds a required measurement-validity limitation.

The strongest defensible wording is:

> **The third cohort prospectively supports excess continuous alignment with the frozen achromatic–chromatic axis beyond a coarse-morph-composition-preserving baseline, while post hoc diagnostics leave brightness/exposure-related white classification as an unresolved measurement-validity alternative.**
