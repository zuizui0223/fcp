# FCP: does better photographed climate opportunity rescue repeatability of species-specific thermal/white slopes?

**24 October 2026. Independent authentic source-photo reliability audit, stacked on PR #271.**

## Why this check

The original 447,500-photo, 4,475-taxon FCP atlas (227,107 classifiable nine-colour photo ROI measurements) showed a repeated positive *across-species mean* photographed white-palette response to shared habitat warmth. However, observer-disjoint 2-fold estimates of individual taxon PHOTO white-climate coefficients were nearly uncorrelated (r=+0.0181 on 1,919 taxa; r=−0.0061 for observer-thinned 1,843 taxa). Those two properties concern different ecological estimands.

The new question is whether this lack of taxon-specific slope agreement is explained by low photographic counts or too little identifiable thermal climate variation, and whether taxa with high opportunity carry reproducible slopes suitable for follow-up biology.

## Source, proxy and tests

Authenticate original source GitHub Actions run 40080001765, artifact ID 12554693263, provider ZIP SHA256 50d5d15a2d145b0db7a176f945435701251768ba614f659d9dac496011b33177. Verify source-derived original JSON SHA256 40b8be90523f876aaf755354fd9527d7429782d8cb42414a69783dedbfa4e4c9 and 17,900-row original photo crossfit CSV SHA256 079241c220a5f5993a17d361a2ea5fc276c7d13e99a168111cc14a3f51ce7253.

The original photo source and original BIO/SRAD climate archives were checked by SHA in the parent PR and are not altered here. All 4,475 taxa remain in both all-observer-photo and one-original-photo-per-observer 4,475-taxon result tables; invalid original photo/climate slopes remain missing, NOT zero.

For named taxa with both observer-disjoint warmth/white PHOTO slopes estimable, define the outcome-blind climatic/photo opportunity proxy I_i as the MINIMUM across opposing source folds of (number of valid heldout source photo records multiplied by the fraction of shared thermal-climate variance remaining after the parent photo model's geography/month/four-other-climate adjustments). This proxy is NOT exact Fisher information or a species-specific slope SE. Neither the original photo colour nor its estimated thermal slope determines its four approximately equally populated opportunity quartiles.

Separately in original source photos and observer-thinned source responses, compute:
- exact same-taxon correlation of the two independently estimated photographic thermal-white slopes;
- fraction of species with concordant slope signs and expected concordance from the marginal positive/negative rates;
- excess sign concordance beyond the empirical chance expectation;
- forward/reverse species-equal mean coefficient and mean difference;
- all these statistics for the entire source matched species set and for four fixed, photo-outcome-blind information quartiles;
- 1,999 whole-genus bootstrap intervals and the difference in slope-correlation and excess sign-agreement between highest versus lowest information quartiles, preserving each genus's taxon membership in a given bootstrap replicate.

Report source five species-disjoint acquisition cohort raw correlations as exploratory checks. The five taxa cohorts used the same photo provider/measurement assay; the two swapped source photo halves are sensitivities, not independent instruments or genetic population replicates.

## Scientific interpretation rules

If higher photographic climate opportunity increases coefficient repeatability, that would identify a photo-based source measurement opportunity condition, not the genetic basis of colour response. If even high-information taxa fail to show reproducibility, then species-level ecological slope rankings in this original photo resource should not be promoted to reliable climate-response phenotype classes.

Do not infer absence of any true biological species-specific thermal response from poor photo coefficient repeatability. Sample sizes per species may be insufficient, residual geographic/camera and measurement biases may remain, and sampling of climate localities is nonrandom. The source photo segmentation has not been fully checked against true petal-specific reflectance; approximately 49.25% original photos were not classifiable. An earlier flower-minus-background BIO5 validation failed (rho −0.00328,p=0.82), and this check does not replace it.

No genuine within-population genetically inherited white flower morph frequencies, actual pigmentation cost, photosynthesis, pollinator fitness or local adaptation were measured. The source study's strongest supported inference remains a small taxon-mean photographic warmth/white composition association, not individually reproducible ecological selection or shared molecular mechanism.

## Reproduction

Source code: scripts/analysis/audit_fcp4475_photo_slope_repeatability_by_climate_opportunity_20261024.py.
CI: .github/workflows/fcp-4475-thermal-fold-reliability-opportunity-20261024.yml.

Outputs: source-authenticated machine-result JSON, 4,475-taxon original and one-photo-per-observer derived CSVs without individual source photos, precise coordinates or observer identifiers. All earlier original PR result files are left unchanged.
