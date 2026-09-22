# H2 post-confirmatory validity diagnostics — 2026-09-22

Status: **POST-CONFIRMATORY VALIDITY / ROBUSTNESS AUDIT; FROZEN H2 VERDICT UNCHANGED**

Reproduction script for the artifact-reconstructable components:\n\n- `scripts/analysis/audit_polymorphism_h2_posthoc_validity_20260922.py`\n\nThis note records diagnostics performed after the third-cohort H2 result was opened. They do **not** replace, rerun, rescue or redefine the prospective H2 decision. The frozen terminal verdict remains `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`. Their purpose is to state more precisely what that verdict establishes and what measurement validity remains unresolved.

## 1. What the confirmatory statistic actually adds

For the primary 0.10 tier:

- isotropic expectation in the eight-dimensional zero-sum subspace: **0.125**;
- coarse-state-preserving structured-null median: **0.457143**;
- observed W: **0.517246**;
- observed / null median: **1.1315**;
- frozen upper-tail p: **0.001**.

The structured null is therefore already strongly aligned with the white-versus-nonwhite contrast. The prospective test confirms the **increment above that construction-controlled baseline**: continuous within-species displacement is more aligned with the fixed white axis than expected after preserving species × coarse-state counts and the coarse-state-specific palette distribution.

It should not be described as if the entire difference between W = 0.125 and W = 0.517 were newly established by the prospective test.

## 2. White-containing coarse modes dominate the admitted set

The immutable primary delta-vector artifact contains 158 species. Recalculation from that artifact gives:

- **137 / 158 (86.7%)** have `white` as either the primary or secondary frozen coarse morph;
- mean per-species `(u_i · q_white)^2` among those 137 species: **0.5622**;
- mean among the 21 species whose two leading coarse morphs do not include white: **0.2239**.

This decomposition is descriptive. It clarifies why the structured null itself has high white-axis alignment and why the inferential target is the residual increment above that null.

## 3. Gate-reapplied structured-null robustness

The frozen structured null conditions on the 158 species that passed the observed continuous minor-cluster gate; it does not reapply that gate in each null world. Because this selection could in principle affect the null distribution, a post-confirmatory diagnostic was run from the 185 species passing the pre-continuous coarse gate.

For each of **299** structured-null worlds, normalized nine-colour rows were permuted within coarse morph as in the frozen null, deterministic two-means was refitted, the continuous minor-cluster fraction >=0.10 gate was reapplied, and W was calculated over the retained null species.

Result:

- observed W: **0.517246**;
- gate-reapplied null median: **0.447495**;
- null mean: **0.447583**;
- 95% interval: **0.428986–0.469580**;
- plus-one upper-tail p: **0.00333**;
- retained vector species per null world: **166–183**, median **176**.

The diagnostic does not weaken H2. Its median is lower than the frozen conditional-null median (0.457143), so the preregistered frozen null is slightly more conservative for this comparison. This result is robustness evidence only; it does not alter the prospective test definition.

## 4. Unresolved white-state exposure / background-context confounding

A digital-highlight validity protocol had previously been frozen for P500 using non-circular image diagnostics (`clip_fraction`, `near_clip_fraction`, and `luminance_q99`). That control was not executed to a biological result for P500 and was not run on the third cohort.

The third-cohort terminal artifact contains flower palette fractions and image hashes but **does not persist image pixels or background palette fractions**. Therefore the prespecified digital-highlight metrics cannot be reconstructed from the confirmatory artifact alone.

A separate post hoc proxy analysis on an auxiliary background-available real-data sample (461 species) found:

- median background white fraction for images classified as white flowers: **0.0692**;
- median for nonwhite-classified images: **0.0553**;
- **66.4%** of species showed the same within-species direction;
- paired Wilcoxon **p = 4.6 × 10^-12**.

This association is consistent with an exposure/scene-brightness or ROI-contamination contribution, but it is **not proof of overexposure**. White-flowered individuals may occur in brighter/open habitats, background composition may covary biologically with morph, and segmentation bleed can also produce flower/background coupling.

Critically, the frozen structured null cannot diagnose an artifact that acts upstream by creating or enriching the coarse `white` state: the null conditions on coarse-state composition. Thus the H2 result supports excess continuous alignment conditional on the measured coarse states, but it does not establish that those coarse white states are artifact-free.

The row-level source used for this auxiliary background proxy is not contained in the immutable third-cohort artifact, so the proxy is reported as a post hoc validity diagnostic rather than promoted into the confirmatory evidence chain.

## 5. disttrait numerical-equivalence boundary

A full-data post hoc audit also tested the later generic `disttrait` implementation against the frozen third-cohort pipeline.

- frozen study-specific W: **0.5172457461**;
- `disttrait.two_mode_axis` over the same frozen 158 species: **0.5183899565**.

The difference arises because `disttrait` defensively renormalizes rows that the frozen FCP loader had already normalized. In near-tied deterministic initializations, floating-point perturbations can change the farthest-point argmax and hence the converged two-means partition.

Two notable cases are:

- *Cirsium vulgare*: |axis cosine| = **0.6333**; frozen cluster sizes 5/45 versus renormalized 11/39;
- *Vicia benghalensis*: |axis cosine| = **0.9830**; frozen 5/44 versus renormalized 6/43.

The scientific decision is unchanged, but `disttrait` is not a bitwise reproducer of the frozen biological H2 pipeline. It should be described as a later general-purpose implementation sharing the estimand and algorithmic structure. The frozen study-specific pipeline and immutable result artifact control the numerical values reported for this paper.

## 6. Reporting consequence

Main-text W values should be rounded to sensible precision (e.g. **0.517**, null median **0.457**) rather than implying ten-decimal measurement precision.

The strongest defensible H2 statement is:

> In a species-disjoint prospective cohort, continuous colour displacement showed **excess alignment with the pre-frozen white-versus-nonwhite axis relative to a coarse-state-preserving structured null**. This increment was robust to reapplying the continuous-cluster gate within null worlds, but digital exposure/background-context confounding of the coarse white state remains unresolved.

