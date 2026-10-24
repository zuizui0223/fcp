# FCP: source-verified nine-hue response direction differs between temperature and precipitation

**24 October 2026 — original-photo empirical results, stacked on PR #266. The original 447,500-source photographs and 4,475 named taxa are unchanged.**

## Result in one paragraph

Nine-hue continuous photographic colour covaried with BIO5 (warmest-month maximum temperature) chiefly in the white-versus-coloured direction: **92.7%** of the squared length of the cross-species **mean nine-palette covariance vector** lay along the pre-fixed white axis (observer-thinned 92.1%). The mean white-source fraction covariance with BIO5 was **+0.003667**, based on **3272** photo/climate-eligible named taxa (observer-thinned **+0.003550**, n=3203). The mean had a positive sign in all five original species-disjoint acquisition cohorts, but only **57.3%** of individual eligible taxa had positive signed white covariances. These estimates are not per-°C slopes, causal changes or inherited genetic flower morph frequencies.

BIO12 (annual precipitation) showed a different exploratory photographic direction: the white-axis projection was not well separated from zero, whereas yellow and bronze exhibited opposing **chromatic-only** associations even after the white coordinate was mathematically removed and after observer thinning. This is a preliminary pattern of source photograph colour composition, **not established petal pigment adaptation to water availability**.

## Authenticated originals and methods

- Original measured source: 447,500 distinct photographs, 4,475 named species, 227,107 originally classified flower-region nine-palette vectors. Exact source photo SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6. WorldClim BIO and 12-month SRAD archives checked against frozen SHA256; original 7 predictor definitions retained.
- An entire source photographic colour vector is compared with a species-standardized historical climate variable after simple linear XYZ geographical and calendar-month adjustment. For signed covariance, the *exact analytic expectation* of permutations restricted within calendar month is subtracted. This estimator is not an individual biological temperature manipulation, common-garden experiment or population-genetic analysis.
- The fixed compositional unit axis q=(8,-1,-1,-1,-1,-1,-1,-1,-1)/sqrt(72) projects the nine measured hues into white-versus-other and an orthogonal chromatic-only eight-hue subspace. These are photographic palette coordinates, not independent biological pigments.
- Original categorical photo labels can be algebraically decomposed into between-label and within-label covariance; neither term is an independently verified genetic-morph frequency or within-genotype plasticity component.
- A species must have >=40 classifiable source photographs with complete original climate normal and valid geography, >=5 occupied WorldClim grid cells, >=20 observers and >=5% original species-level climate variance remaining after broad geography/month projection. The source genus-bootstrap uses 1,999 whole-genus redraws. All 4,475 taxa remain in the output even if photo-climate estimation is impossible.
- Simultaneous 7-variable white-direction genus intervals, 7x9 hue-component intervals, and 7x8 chromatic-only intervals are distinct **exploratory** families; no omnibus adjustment across all component/source panels. The colour direction was examined after prior FCP results, so its new climate-specific patterns are post hoc.

## Seven-climate photographic white-component direction

The last column is the **simultaneous genus-bootstrap 7-climate interval for the white-axis projection** (scaled relative to the white coordinate by the fixed unit vector); do not confuse this with a per-unit temperature slope.

| Climate normal | Source taxa | White-coordinate mean covariance | Squared white-axis share of cross-species mean hue vector | q-projection genus simultaneous 7-climate 95% interval |
|---|---:|---:|---:|---:|
| BIO1 | 3038 | 0.002388 | 89.9% | [0.001500, 0.003566] |
| BIO12 | 3266 | -0.000530 | 15.8% | [-0.001563, 0.000440] |
| BIO14 | 3245 | -0.000268 | 9.9% | [-0.001240, 0.000672] |
| BIO15 | 3184 | -0.000297 | 28.4% | [-0.001265, 0.000634] |
| BIO4 | 2943 | 0.001709 | 91.2% | [0.000887, 0.002740] |
| BIO5 | 3272 | 0.003667 | 92.7% | [0.002873, 0.004906] |
| SRADmean | 2382 | 0.000193 | 5.8% | [-0.000677, 0.001086] |

For BIO5, the individual nine-hue source components with a 7x9-hue max-T interval excluding zero were: **pink -0.000849; purple -0.000758; white 0.003667; yellow -0.001106**. This is the observed multivariate photo-palette response, *not evidence that any plant genus uniformly changes its petals among these four colours*.

## Matched-species comparison of **two different climatic hue directions**

For **3245 identical original eligible taxa**, the squared white-axis share of the *source mean climate–hue vector* was **92.9% for BIO5** versus **18.1% for BIO12**. The paired difference was 74.8%, whole-genus exploratory 95% interval [0.369312, 0.924871] in fractional units (original): expressed in percentage points, **74.8 [36.9, 92.5]**.

After one original photo per observer, the identical-climate eligible set comprised **3178 taxa**: BIO5 white-axis share 92.4%, BIO12 23.0%, paired difference **69.4 percentage points [32.5, 91.3]**. This is a comparison of squared mean vector alignments, **not percentages of genetic variance or environmental causation**. The BIO12 mean hue vector is weak and its stand-alone orientation is uncertain, so the magnitude and interpretation should not be overclaimed.

For BIO12, the chromatic-only source components whose 7-climate × 8-hue simultaneous interval excluded zero were **bronze -0.000676 (simultaneous56 [-0.001299, -0.000053]); yellow 0.000955 (simultaneous56 [0.000185, 0.001724])**. The same signs were present after observer thinning. The yellow/bronze contrasts were also present within the original coarse photographic colour bins, so they are not solely an artefact of changing coarse category frequencies; image colours, flower organ and photo conditions remain unresolved.

## Five species-disjoint source cohorts: descriptive check

Entries are signed source photo-hue/climate covariance means, not slope effects. Yellow and bronze are **chromatic-only within original photo labels**, after removing the algebraic white-axis projection.

| Original sampling cohort | BIO5 taxa | BIO5 photo white | BIO12 within-label yellow | BIO12 within-label bronze |
|---|---:|---:|---:|---:|
| discovery | 368 | 0.003902 | 0.000684 | -0.000411 |
| validation | 357 | 0.004746 | 0.000200 | -0.000499 |
| prospective | 376 | 0.003080 | 0.000593 | -0.000348 |
| wave4 | 340 | 0.001885 | 0.001291 | -0.001211 |
| new_U100 | 1831 | 0.003861 | 0.000640 | -0.000802 |

**All five cohort means** had positive BIO5 white and BIO12 within-photo-label yellow directions, and negative within-photo-label bronze directions. This sign pattern was retained after one-photo-per-observer thinning. But cohort-level genus confidence intervals often included zero, especially for BIO12; the cohort chromatic directions were identified *after* inspecting the full 7×8 family, so these are not untouched prospective tests. Wave4 BIO5 white was positive but its cohort-only 95% genus interval includes zero.

## Claim ceiling and important negative evidence

- The 9-hue palette describes segmented PHOTO ROIs, with no independent botanical adjudication across the full source. Previously retrieved authentic photos from 27 putative two-peak species already demonstrated possible petals vs buds, persistent sepals/fruit, nectar-guide-like basal blotches, lighting and taxonomic concerns. Thus this does not certify genetically maintained colour polymorphisms or anthocyanin states.
- Historical photo-classification failure is ~49.25%. A photo-source hue can covary with climate through observer geography, lighting, camera white balance, cultivation or floral developmental stage as well as real pigmentation. Five taxon cohorts are disjoint in species but share the same source and colour-classification assay.
- The earlier independent validation flower-minus-background BIO5 test was **not supported** (mean rho −0.00328, p=0.82); this follow-up does not introduce an independent matched-background or experimentally calibrated petal-reflectance control.
- The BIO1/BIO4/BIO5 and precipitation variables are correlated. Their contrasting *photo palette* covariance directions do not establish distinct environmental-selection mechanisms. Historical climate normals are not temperatures at anthesis and there are no individual phenotype-linked genes or fitness data.
- The 1°-region observer-thinned null remains in the separate PR #266 as a supplementary spatial-scale sensitivity; it is neither silently deleted nor made the sole pass/fail gate for this broader photo-colour climate question.

## Reproducibility and machine receipts

- Original-source nine-palette direction script: scripts/analysis/fcp4475_nine_hue_climate_direction_20261024.py
- Independent original-CSV and five-cohort audit: scripts/analysis/audit_fcp4475_five_cohort_climate_hue_direction_20261024.py
- SHA-pinned original source run: https://github.com/zuizui0223/fcp/actions/runs/40060153598; artifact 12547087651 (ZIP SHA256 c23f0c7452019940bdbcca88b249be854c3f236e39ce189aa8c77119f021f778). The source has all 4,475 taxon-level nine-hue × 7 climate × between/within coefficients.
- Exact independent species-disjoint cohort re-read: https://github.com/zuizui0223/fcp/actions/runs/40060576651; audit artifact 12547212826 (SHA256 a158ac3d300884513c643f959c8ed4025c0607671be5b8f06c239fcbc3d34e38).
- Original full numeric source JSON: results/fcp4475_climate_nine_hue_direction_verified_20261024.json
- Independent complete cohort and matched-temperature-precipitation comparison JSON: results/fcp4475_climate_hue_five_cohort_verified_20261024.json

## Safe manuscript result text

> In the original within-species continuous nine-colour photographic data, warm-season temperature was associated chiefly with the recurring white–chromatic colour direction: approximately 93% of the squared cross-species mean BIO5–colour covariance vector aligned with a fixed white contrast, although only about 57% of individual taxa showed a positive white association. In contrast, annual precipitation was accompanied by a weaker but directionally distinct chromatic-only source signature, involving yellow and bronze palette components, after removal of the white component. These source-derived relationships are exploratory photographic ecological patterns, not evidence of causal climatic colour selection, stable genetic morph frequencies, shared pigment metabolism or biosynthetic costs.
