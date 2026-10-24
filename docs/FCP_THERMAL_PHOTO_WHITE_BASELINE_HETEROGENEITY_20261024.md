# FCP: ecological heterogeneity of shared warmth–photographic white colour

**24 October 2026. Source-locked analysis stacked on PR #269. The underlying source-photo measurements and previously published inference are unchanged.**

## Biological motivation

Across the 447,500 genuine measured photos, 4,475 named plant species and 227,107 originally classifiable nine-hue photo vectors, the source verified shared BIO1–BIO5 warmth axis has a positive species-equal mean photographic white-colour slope in 2,837 source taxa. However, only approximately 54% of those individually identifiable species have a positive coefficient. A positive mean can conceal sharply different photographic strategies.

The new question is: **Which named species have stronger photographed white–warmth association, and how does the association differ with existing white photo visibility, photographic colour diversity, classification probability and actual observed climate breadth?**

## Source traits, calculated independently of model slope selection

From ALL 100 original photos for every source taxon, record the original number whose photo colour was classifiable. From the classified originals only, calculate the mean original nine-palette white ROI fraction, original photo coarse-white category counts, source-class photographic colour diversity and 4-class counts each with >=5 photo records. None is assumed to measure an individual wild genetic floral morph.

Four mutually exhaustive original photographic white-visibility strata:
- no_white_label: zero classified white-category source photos;
- rare_white_label_1to4: one to four white-category photos;
- both_white_and_nonwhite_ge5: five or more original white-category photos and five or more original nonwhite-category photos;
- nearly_all_white_label: five or more white-category photos, fewer than five nonwhite-category photos.

Use the same immutable WorldClim 2.1 BIO archive to measure the actual original photo-sampled warmest-month maximum temperature breadth within each species: q90-minus-q10 and standard deviation, expressed in Celsius. Source locations are sampled citizen-science photo localities, not the native range or thermal physiology of individual flowers.

The original thermal model already measured each species' sampled number of distinct 10-arcminute raster climate cells and remaining, uniquely identifiable shared thermal climate variance. These are diagnostic covariates, not biological mechanisms.

## Analysis and inferential checks

Join the old source PHOTO warmth/white partial-slope coefficients on exact inat_taxon_id. Keep a complete 4,475-taxon ledger even where photo/classification/climate support is inadequate; do not impute zero effects.

For both original and colour-blind one-photo-per-observer panels:
- Report each photographic white-category group's exact source species counts, equal-species mean/median warmth-white coefficient, positive-species fraction and 999 whole-genus bootstrap marginal 95% intervals.
- Summarize exploratory equal-size quartiles of continuously measured baseline photo white fraction, photo classifiability among the original 100, observed BIO5 sample breadth in Celsius, raster-cell sampling breadth and unique residual warmth variance.
- Fit a descriptive species-equal regression of photo warmth-white coefficient on the three photographic white-category indicators (reference no_white_label), while adjusting for original photo-classifiability, number of raster cells, climatic breadth, remaining identifiable shared temperature variance and the original five taxon-acquisition cohorts. Bootstrap entire genera for 999 re-fits; do not count photographs as independent replicates.
- Distinguish source photo traits from genetic inheritance, and different model denominators from true zero/absent biological variation.

**Nonindependence caveat:** Photographic white phenotype grouping and the thermal response use the SAME original photographs. Boundaries of the hue fraction, photo selection/observer/camera/lighting, flower age, botanical organ, cultivated taxa and missingness can produce group associations. This analysis tests an observable, potentially informative ecological modifier; it DOES NOT prove genetically discrete morphs, thermal pigment adaptation, local phenotypic selection, anthocyanin investment cost or causal fitness effects.

## Authentic source and execution

Original photo source compressed CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6.
Original WorldClim 2.1 BIO archive SHA256: 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf.
Original warmth source-model run 40066224301, GitHub artifact 12549840336, provider artifact ZIP SHA256: 3f290b8028bbfcdb6fa2d1a4ae85f3155b21786d7c8d48175af4805c78c099a0.

A GitHub Action downloads the immutable photo/BIO release and authenticates source photo/BIO bytes. For the upstream photo warmth model, the workflow verifies the *provider's original ZIP digest* and separately all five immutable extracted source CSV/NPZ/JSON file checksums; no reconstructed ZIP is claimed byte-identical.

Script: scripts/analysis/fcp4475_thermal_white_baseline_photo_opportunity_20261024.py.
CI: .github/workflows/fcp-4475-thermal-white-baseline-opportunity-20261024.yml.

Output must include complete original 4,475-source-taxon species summaries and machine-readable aggregate results, without redistributing original individual photo IDs, observer identities or exact photographic coordinates. Previously unsupported source flower-minus-background temperature test remains an explicit limitation rather than an artificial primary acceptance gate.
