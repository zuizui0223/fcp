# FCP: thermal white–chromatic photo-colour direction survives moisture + solar adjustment, not fully separable from other thermal normals

**24 October 2026. Source-authenticated original 447,500-photograph analysis.** This is an exploratory analysis of original PHOTO ROI hue compositions and WorldClim long-term climatic normals, not independently botanically verified fresh flower petals, genetic flower-colour morph frequencies or adaptive fitness.

## Main finding and scientific boundary

The original nine-colour × seven-climate analysis in PR #267 revealed that the cross-species mean PHOTO hue-climate direction aligns strongly with a white-versus-chromatic axis for warmest-month maximum temperature (BIO5), but the annual-precipitation (BIO12) direction was less white-aligned and involved a yellow/bronze photographic hue contrast. Because climate normals covary, PR #268 directly tests whether those two PHOTO colour responses contain distinct climate information.

**Result:** After conditioning BIO5 on annual precipitation (BIO12), mean solar radiation, broad geographic XYZ site coordinates and observation calendar month, the BIO5–white photo fraction slope-like association remains positive across source taxa. The result persists with colour-blind one-photo-per-observer thinning, in the original five disjoint acquisition cohorts and within the original coarse PHOTO labels. On the other hand, residual variation unique to BIO5 after **all six other WorldClim/solar variables** is scarce in many species, and that more strictly identified individual BIO5 effect is unresolved. The positive BIO12 yellow/bronze photographic contrast remains directionally present as a point estimate but its climate×chromatic simultaneous intervals overlap zero after joint adjustment.

The ecological inference is that a recurring warm-climate photographic colour direction is identifiable **apart from the broad precipitation and solar gradients**, whereas uniquely attributing it to BIO5 rather than related thermal indices (BIO1/BIO4) is **not resolved** with these source photo/climate data. This is *not* evidence that high temperatures cause genetically heritable white petals or that rainfall causes yellow/bronze anthocyanin/carotenoid biochemical changes.

## Exact original sources and methods

- Immutable original source: 447,500 measured iNaturalist photographs, 4,475 named plant species, of which 227,107 retain classifiable original flower ROI nine-hue palette measurements. Original source CSV SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6.
- Original WorldClim 2.1 BIO1, BIO4, BIO5, BIO12, BIO14, BIO15 and 12-month mean SRAD, pinned source archive SHA256 BIO 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf, solar c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.
- Per species, center and scale climate by observed original-source within-species SD, remove linear XYZ geographical site coordinates and original photo calendar month. Next isolate the residual of the target climate variable with respect to other climate axes **named in the model**:
  (A) BIO5/BIO12/SRAD conditional on the *other two* (not on all temperature/precipitation normals), versus
  (B) each of seven climate variables conditioned on **all six others**.
- At least 40 classifiable complete-geographical/climate original photos, 20 original photo observers and 5 occupied WorldClim 10-arcminute cells; an environmental axis must retain at least 5% of the original standardized climate variance after all specified linear controls. Explicit not-estimable labels for taxa failing support and environmental identifiability. Every original 4,475 source species remains in the support CSV, with no unclassifiable photo hue imputation.
- For the nine original colour palette fractions, compute the signed PHOTO colour/climate covariance minus the analytically exact calendar-month-conditional permutation expectation, then divide by the **unique residual climate variance**. This is a descriptive conditional partial-slope-like statistic in species-original-climate-SD units, **not a °C slope or causal response**.
- Exactly partition each per-species covariance into original categorical PHOTO-label between and within components; original labels require ≥5 photos for at least two photographic hue categories for a between-label estimate. Source labels and floral organs have not been individually biologically or genetically validated.
- Compute signed white coordinate and a fixed nine-colour q=(8,-1,-1,-1,-1,-1,-1,-1,-1)/sqrt(72) white-versus-other direction, and its eight-dimensional orthogonal chromatic-only residual. Bootstrap complete original taxonomic genera (999 replicates); apply simultaneous 7-climate max-T limits for q-directed photographic colour in a given model/source component. Separately, an independent, artifact-SHA-locked analysis applies max-T to 24 climate×eight-chromatic tests for the core model and 56 for the seven-climate model, per original/photo-component/source panel.
- Five species-disjoint original photographic cohorts are descriptive replication checks. Both original source photos and 1-photo-per-observer thinning share the same iNaturalist source and ROI segmentation assay.

## Result 1. Thermal BIO5–white photo signal persists after BIO12 and solar adjustment

| Source model | Eligible species (BIO5) | Mean species-equal white-palette fraction conditional slope-like association | White-axis squared alignment of mean nine-hue vector | q-axis genus 7-climate simultaneous 95% |
|---|---:|---:|---:|---:|
| BIO5 conditional on BIO12 + SRAD, original photos | **3,210** | **+0.009110** | **94.4%** | **[+0.006466,+0.012860]** |
| Same, one photo per observer | **3,137** | **+0.008881** | **95.1%** | **[+0.006249,+0.012591]** |
| BIO5 conditional on **all other six climatic variables**, original | 622 | +0.004006 | 49.9% | [−0.007708,+0.016206] |
| Same full seven conditional, one photo per observer | 606 | +0.005488 | 58.4% | [−0.006453,+0.018095] |

The q-axis interval (last column) is for the **fixed scaled white–chromatic projection**, not the numerically different unit of the unscaled white palette fraction slope in column 3. It is not a real-temperature, whiteness-percent or genetic colour-morph effect.

The core thermal sign occurs among all FIVE pre-existing species-disjoint acquisition cohort means (discovery, validation, prospective, wave4, and new_U100), although taxon-level positive frequencies are only **55.7%** original and **55.5%** observer-thinned. Thus this is a recurring **cross-species mean photographic direction**, not a positive BIO5 response in every taxon. After fully conditioning on the other six axes, only 622/606 species retain enough unique BIO5 source-climate variance for a coefficient; this comparison cannot establish a uniquely attributable causal BIO5 response.

## Result 2. At the same original species, the thermal association spans both coarse-photo-label distribution and within-label hue

Use only original species with at least two original photographic colour classes, each backed by ≥5 photos, where all three algebraic components are separately estimable for BIO5.

| BIO5 slope-like white fraction, adjusted for BIO12 + SRAD | Original photo species | Observer-thinned |
|---|---:|---:|
| Exact same-taxon denominator | 1,543 | 1,462 |
| Total photo-white conditional slope | +0.013719 | +0.013407 |
| **Between original coarse PHOTO labels** | +0.010560 | +0.010545 |
| **Within original coarse PHOTO labels** | +0.003159 | +0.002862 |
| **Between percentage of the signed mean partial-slope-like association** | **77.0%** | **78.7%** |
| 999 genus bootstrap 95% for that algebraic ratio | [65.2%,88.3%] | [66.1%,90.6%] |

These three components must be compared **on the same original species**. Comparing between-label on 1,543 species with total on 3,210 species would be statistically invalid. The ratio is an algebraic source-PHOTO-class composition measure, **not** the fraction of genotypically discrete colour morphs or genetically inherited pigmentation changes. Within-class photo hue includes lighting, aging and organ/ROI uncertainties, not just biological continuous expressivity.

## Result 3. The proposed precipitation chromatic direction is still exploratory, not a unique confirmed climate response

Under BIO12 conditional on BIO5 + SRAD, the photographed BIO12–white association is not reliably separated from zero (mean source +0.001706 for n=3,180 original taxa, observer-thinned +0.001738 for n=3,113). The pooled **chromatic-only** BIO12 source coefficients preserve the earlier exploratory pattern (yellow positive, bronze negative), but whole-genus simultaneous 24-test confidence intervals overlap zero:

| BIO12 chromatic-only photo hue after BIO5+solar adjustment | Original n=3,180 | 24-test genus simultaneous 95% | Observer thinned n=3,113 |
|---|---:|---|---:|
| Yellow | +0.001132 | [−0.001049,+0.003313] | +0.001064 |
| Bronze | −0.001518 | [−0.003282,+0.000245] | −0.001768 |

Under the more stringent all-seven predictor model, yellow and bronze **still have the same point-estimate signs** but their 56-test familywise intervals also overlap zero; the effect cannot be claimed as uniquely precipitation-specific. This is a meaningful boundary on the original exploratory photo-colour axis comparison, not evidence the original unconditional vector measurements were fabricated or false.

## Result 4. Matched-source climatic directions under two distinct conditioning levels

To avoid different species support driving the white-axis energy comparison, pair *exactly* taxa with identifiable photo-colour response vectors for both BIO5 and BIO12.

| Same-taxon photo direction comparison | n taxa | BIO5 white-axis share | BIO12 white-axis share | Paired difference genus 95% |
|---|---:|---:|---:|---:|
| BIO5/BIO12/SRAD core conditional, original | 3,123 | 94.9% | 20.1% | [+23.3,+94.7] percentage points |
| Core, observer-thinned | 3,053 | 95.5% | 19.7% | [+26.3,+94.9] percentage points |
| Full seven axes mutually conditional, original | 504 | 38.7% | 47.5% | [−68.2,+59.5] percentage points |
| Full seven, observer-thinned | 486 | 24.3% | 39.9% | [−67.0,+63.6] percentage points |

All percentages here describe **squared alignment of a cross-species mean photographic hue covariance vector** with a fixed visual white axis, never explained genetic variance or pigment-cost allocation. In the fully conditional model, residual BIO5 climate information is sparse and the orientation ratios are unstable; do not treat this subset as proof of an opposite climatic evolutionary mechanism.

## Claim ceiling

- **No independent botanical flower-petal review** for the full 447,500 original photographs; 220,393 original photo records could not be classed into the frozen hue groups. Visual screening among original model bimodality candidates already identified bud/fruit/sepal/petal blotch and stage/exposure concerns.
- **No full matched-background control** on this exact 447,500-source climate model. The earlier source-disjoint validation flower-minus-background BIO5 sensitivity failed (rho −0.00328, p=0.82). This remains important and is not overridden by statistical photo covariance.
- The original environmental predictors strongly covary through geography, evolution and selection of source photo localities. A unique thermal driver is **not** established, although the white photo response persists when controlling annual precipitation and solar radiation.
- WorldClim norms are long-term source locality climatic context, not observed temperature at anthesis. Geographic XYZ linear controls do not remove all unobserved spatial/phylogenetic history or sampling/cultivation.
- The two proposed climate response model families, source photo label components and cohort bootstrap comparisons are exploratory and not globally adjusted across all analyses.
- There is no causal local adaptation, within-population genetic morph frequency, pigment maintenance cost, biochemical gene expression or pollinator-mediated fitness evidence in these data.

## Source / reproducibility

- Original source climate conditional computation: scripts/analysis/fcp4475_joint_conditional_climate_hue_20261024.py.
- Original photo and climate archive SHA CI, including exact source number, original photo colour labels and synthetic known confounding tests: https://github.com/zuizui0223/fcp/actions/runs/40061576683 — PASS, original full taxon coefficient/source manifest artifact ID 12547264579, ZIP SHA256 ce145c38144d833be48b1a00a21414f6fb9e6d234e7101c5183e0f5d9e5dede9.
- Independent original source artifact-SHA-locked 24-/56-hue and matched species-photo-label genus bootstrap: scripts/analysis/audit_fcp4475_joint_conditional_chromatic_20261024.py, CI https://github.com/zuizui0223/fcp/actions/runs/40062146684 — PASS, artifact 12548011674.
- Machine-readable complete genuine-source conditional photo-climate results and provenance: results/fcp4475_joint_climate_hue_source_verified_20261024.json.
- Machine-readable independent chromatic conditional/identical-taxon audit and provenance: results/fcp4475_joint_climate_chromatic_audit_verified_20261024.json.
- Underlying four complete all-4,475-taxon support CSVs and four original 4,475×7×3×9 measured-photo numerical arrays are retained as an exact GitHub Actions artifact with no original raw photo IDs, geographical coordinates or observer identities.

## Scientifically bounded manuscript statement

> Across the original within-species photographic nine-hue data, the warmest-month temperature-associated visual shift toward white remained directionally positive after accounting for annual precipitation and average solar radiation, and most of the corresponding signed photo-white covariance was attributable to variation among coarse photographed colour categories. The same data did not isolate a distinct BIO5 effect when all other correlated climatic normals were held constant, and the previously observed precipitation-related yellow/bronze photographic direction was not independently resolved under climate-adjusted multiple-hue uncertainty. These findings identify repeated climate-associated photographic colour geometry, rather than climate selection on genetically inherited floral colour morphs or pigment biosynthesis costs.
