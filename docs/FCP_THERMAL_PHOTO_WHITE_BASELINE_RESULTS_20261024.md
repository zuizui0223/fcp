# FCP: original photo-white coexistence modifies the warmth-associated photographic colour signal

**24 October 2026 — source-authenticated original 447,500-photo observational analysis. Source ZIP, model numeric results and WorldClim BIO all SHA256-checked.**

## The ecological question and what is actually measured

The earlier shared thermal-level model (source PR #269) found a repeated small positive species-equal mean BIO1/BIO5 temperature-axis association with the photographed nine-palette **white fraction**, but only about 54% of source species showed positive individual coefficients. This analysis joins original photo-white *baseline measurements* to those immutable taxon-wise coefficients and asks **which species carry the stronger PHOTO colour-climate relationship**.

Original immutable source 447,500 photographs, 4,475 named taxa, exactly 227,107 classifiable original nine-palette ROI images (220,393 photos remain without classified hue). The baseline predictor and thermal response are computed from the **SAME original source photos**. White refers to the original photographic source *coarse white colour label*, not an independently verified heritable white petal genotype.

Four predeclared source-photo categories:
- no_white_label: no original source photos classed white;
- rare_white_label_1to4: one to four original source white-classed photos;
- both_white_and_nonwhite_ge5: >=5 white-classed photos AND >=5 nonwhite-classed photos;
- nearly_all_white_label: >=5 white-classed photos with fewer than five nonwhite-classed photos.

Photographic colour classification is not random. It can vary with camera, image background, floral age and organ, cultivation or season, and can produce circular descriptive associations when baseline and slope share photos.

## Result A. The stronger thermal-white photo pattern is in species with both source photo colour classes

| Source-photo white visibility | Source identifiable taxa | Species-equal mean partial-slope-like BIO1/BIO5 warmth–photo-white coefficient | 999 whole-genus bootstrap marginal 95% |
|---|---:|---:|---|
| No photograph classed white | 173 | +0.002239 | [−0.001476,+0.006183] |
| Rare white label: 1–4 photos | 670 | +0.001475 | [−0.003382,+0.006242] |
| **Both white and nonwhite classes: >=5 photos each** | **1,359** | **+0.013606** | **[+0.007757,+0.019206]** |
| Mostly/all white-classed | 635 | +0.005391 | [+0.000867,+0.010421] |

Source-colour-blind one-photo-per-observer sensitivity yielded:
- no white: 171 taxa, +0.003059, CI [−0.001040,+0.007025];
- rare white: 662 taxa, +0.001971, CI [−0.002525,+0.006387];
- **both source PHOTO colours**: **1,309 taxa, +0.013243, CI [+0.006967,+0.019262]**;
- mostly/all white: 629 taxa, +0.004786, CI [+0.000093,+0.009851].

These are group means with non-independent source-colour category assignment; intervals are *unadjusted for multiple exploratory subgroup tests*. The 4-class labels are not true genetic morph counts.

## Result B. Does photo sampling and climatic breadth account for the group contrast?

We adjusted the source-taxon conditional white/photo warmth coefficient by four jointly modelled photo colour group indicators (no-white source photo label as reference), original classifiable photo count among each 100 photographs, original climate raster cells, photographed BIO5 climatic q90–q10 range in Celsius, residual unique shared thermal climatic variance and pre-existing five species-disjoint cohort fixed effects.

| Descriptive observational contrast (both-source photo categories versus no source white photos) | Point estimate | 999 genus bootstrap marginal 95% |
|---|---:|---|
| **Original photos** (2,837 taxa) | **+0.008141** | **[+0.000615,+0.015923]** |
| One original photo per observer (2,771 taxa) | +0.008030 | [−0.000203,+0.015631] |

The positive contrast remains **exploratory**: in observer-thinned data the genus interval *slightly overlaps zero*, so it would be incorrect to claim that the photo group effect is securely established under every resampling design or multiplicity choice.

Photo classifiability, original WorldClim raster cell count, photographed BIO5 q90–q10 climate breadth and identifiable residual climate fraction had covariate-adjusted genus intervals spanning zero in both panels. This suggests that **no simple linear explanatory relationship with these source-opportunity proxies** explains the slope on its own, but does **not** demonstrate source selection/camera/colour-error has been removed.

## Result C. The continuous white palette shows a possible non-linear association

The original PHOTO-ROI *mean nine-palette white fraction*, rather than coarse 4-state count, was divided into **exploratory quartiles** among the 2,837 eligible original source species:

| Mean original PHOTO white-fraction quartile | n species | Mean thermal-white photo coefficient | 999 genus marginal 95% |
|---|---:|---:|---|
| Q1, ~0.009–0.119 | 710 | +0.002099 | [−0.001619,+0.005949] |
| Q2, ~0.119–0.396 | 709 | +0.005555 | [−0.001738,+0.012835] |
| **Q3, ~0.396–0.784** | **709** | **+0.019081** | **[+0.011051,+0.027675]** |
| Q4, ~0.784–0.977 | 709 | +0.006110 | [+0.000853,+0.011095] |

A similarly elevated Q3 source-photo group occurred in the observer-one-photo panel (n=692, +0.017703, genus95 [+0.008299,+0.026412]). This is a **candidate non-linear photographic colour opportunity pattern**, not a confirmed frequency-dependent biological selection optimum. The range of a variable bounded between 0 and 1 mechanically constrains its possible thermal responses near the extremes, and conditioning on observed photographs can itself create this shape.

## Result D. Photo classifiability and real climate breadth are not simple monotone correlates

The observed original BIO5 temperature breadth was evaluated from the authentic site-level WorldClim BIO5 normals for classifiable, valid-photo locations, as q90–q10 **in °C**.

The positive thermal/photo-white mean was observed across photographic climate breadth groups, but the four equal-size broad-temperature groups showed no simple monotonic increase or decrease:

- Narrowest q90–q10 width up to 4.93°C: n=710, mean +0.012321.
- Next 4.93–6.61°C: n=709, mean +0.004079.
- Next 6.61–8.28°C: n=709, mean +0.010199.
- Widest over 8.28°C: n=709, mean +0.006231.

Original source photo classifiability rate likewise showed non-monotone thermal slope means: from +0.009740 (lowest classifiability quartile) to +0.014917, +0.005233, +0.001710 (highest photo classifiability quartile). The thermal signal cannot be treated as independent of photo missingness from these numbers. Classification success depends on observed and potentially biological/technical colour properties and photo settings.

## Interpretation appropriate for an ecological manuscript

> The warmth-associated white-palette signal varied among taxa: in this photographic resource it was strongest where both white and nonwhite **photographic** classes were repeatedly observed, and in the intermediate-high rather than extreme baseline white-palette range. The source observation is consistent with stronger *measurable photo colour redistribution opportunity* in species showing multiple photographed white/nonwhite states. It cannot establish genetically maintained flower colour polymorphism, stronger climatic selection on actual flower petals, or an anthocyanin production cost–benefit trade-off. The association could partly reflect bounded hue fractions, source image classification and shared-photo dependence; botanical validation and independent image/background controls are required to support a biological morph claim.

This is a stronger, testable *species-level ecological effect-modifier observation* than simply claiming every species becomes white in warm habitats, but it remains an **observational photographed trait result**.

## Source provenance and scripts

- Original photo source measurement CSV SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6, and original WorldClim BIO archive SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf.
- Frozen source-model seven WorldClim group-axis run 40066224301, original GitHub artifact 12549840336, SHA256 3f290b8028bbfcdb6fa2d1a4ae85f3155b21786d7c8d48175af4805c78c099a0.
- This analysis computes baseline whiteness from **the actual 447,500 original source photo records**, joins frozen 4,475-taxon white/temperature partial slopes by source inat_taxon_id, and preserves explicit unclassified photo data with no hue imputation.
- Script: scripts/analysis/fcp4475_thermal_white_baseline_photo_opportunity_20261024.py.
- Workflow: .github/workflows/fcp-4475-thermal-white-baseline-opportunity-20261024.yml.
- Authenticated source CI run: https://github.com/zuizui0223/fcp/actions/runs/40068289947 — **PASSED**.
- Source results: results/fcp4475_thermal_photo_white_baseline_verified_20261024.json and two original 4,475-taxon derived source-photo opportunity CSVs in the Actions artifact.

**Limits from previous independent validations remain:** the previous flower-minus-background climate test on a different original species-disjoint set did not support a flower-specific BIO5 effect (rho −0.00328, p=0.82). This new source analysis does not correct photograph backgrounds or isolate fresh petal organs, and does not use genetically verified native population morphs or flower fitness.
