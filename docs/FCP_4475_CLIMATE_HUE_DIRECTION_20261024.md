# FCP 4475: shared directions of climate-related continuous photo hue

**24 October 2026 — source-locked exploratory analysis, stacked on PR #266.**

## Question

Does the repeated temperature association of continuous flower-photo colour primarily involve the **white-versus-chromatic axis**, or also a common redistribution among the *eight nonwhite hues*? Which photo-colour directions replicate across the source five species-disjoint cohorts?

Earlier PR #265 estimated positive photographed white-palette/BIO5 covariance after broad geography/month controls on 1,578 source species with two adequately supported original photo labels. It decomposed 75.8% of signed mean covariance between the originally measured FOUR **photographic** classes. The geographically conditional work in PR #266 supports a positive 5° association and keeps the finer 1° sensitivity as supplementary, not as a universal acceptance rule.

## Inputs and identifiability

- Authenticate original **447,500 photo records**, **4,475 source taxa**, **227,107 classifiable nine-palette** flower ROI vectors by the frozen SHA256 and unchanged measurement model.
- The nine photo palette fractions are white, yellow, orange, red, pink, magenta, purple, blue, bronze. They sum to one.
- Seven climatic measurements are immutable WorldClim 2.1 BIO1, BIO4, BIO5, BIO12, BIO14, BIO15, plus the 12-month mean solar radiation, each sampled at the photographed coordinates from the original SHA-verified release.
- Minimum 40 complete-classified/valid-coordinate source photos, 5 occupied WorldClim raster cells, 20 distinct original observers; original photo month known; a variable must retain ≥5% of source climatic variance after species-level geographic/month adjustment. All **4,475** source taxa remain in the output, with nonidentifiable species explicitly recorded, not zero imputed.
- Run once with all available originally classified photographs and once with colour-blind one-photo-per-observer subsampling; five source-disjoint taxon cohorts reported separately but still share the same photo provider and assay.

## Statistics

Each species' climate vector is standardized using observed within-species climate SD; project out full calendar-month dummies and broad global Earth-centred X/Y/Z coordinates. For a climate axis, calculate signed covariance with all original nine continuous photo palette fractions and subtract the **exact analytical expected covariance under all within-calendar-month random climate reassignments**. This exact source permutation expectation is validated against exhaustive synthetic permutations; no p-value is assigned to the exact expectation alone.

The nine-dimensional covariance vector always sums to zero. A pre-fixed unit white–other axis q=(8,−1,−1,−1,−1,−1,−1,−1,−1)/sqrt(72) defines a mathematically exact orthogonal split. Project each climatic colour covariance vector along q (signed white-versus-nonwhite association) and into the orthogonal eight-nonwhite-hue residual (a contrast among nonwhite colours independent of the white coordinate).

Separately decompose photographed covariance according to the pre-existing coarse source photo classes into (a) between original photographic classes and (b) within original photographic classes. The between component requires ≥2 original photo colour classes with ≥5 original observations each. These terms are **not heritable morph frequency versus genetic plasticity**.

Estimate species-equal mean vectors and bootstrap complete genera (1,999 resamples), with a 7-climate max-T simultaneous family for signed white direction, a 7-climate × 8-hue max-T family for chromatic-only residual and a 7-climate × 9-hue simultaneous family for individual source hue coefficients. Source-original/observer-thinned and three photographic label components are distinct exploratory analyses rather than a single jointly corrected omnibus test.

For each climate, report the mean signed photographic white coordinate, within-species median, fraction of species with positive direction, and taxon-cohort direction replication, not just statistical significance. Report the squared white-axis share of the **cross-species mean source climatic-photo covariance vector's energy**; it is neither a proportion of inherited phenotypic variance nor the fraction of plant pigment investment explained by climate.

## Scientific limits

This is **photo-derived coloured ROI composition**, not botanically reviewed fresh open petal colour for all photos. Prior direct-photo screening shows flower bud/fruit/sepals/nectar-guide blotches and camera lighting can contribute. Missing source classifications are common (~49.25%), the same iNaturalist provider is used, and WorldClim historical normals are not manipulated actual developmental temperatures. Prior flower-minus-background BIO5 validation was **not supported** (rho −0.00328, p=0.82); this new study cannot erase that. Geographic and climatic predictors are correlated and controls remain broad, not nonparametric spatial/neutral-ancestry adjustment.

Retain the prior 1° observer-thinned region sensitivity as **supplementary provenance** rather than making an arbitrary region width the paper's main pass/fail rule. This study does not claim genetic morphs, adaptive costs/benefits, pollinator preferences, local fitness, shared anthocyanin machinery or causal climatic differentiation.

## Code and reproducibility

The new analysis is implemented in the Python file scripts/analysis/fcp4475_nine_hue_climate_direction_20261024.py, and source-exact GitHub Actions workflow .github/workflows/fcp-4475-climate-hue-direction-20261024.yml.

The CI must confirm exact original input checksums, exact permutation-expectation mathematical fixture, full 4,475-taxon source-denominator coverage, source-cohort coverage, non-imputed nine-palette mass, genus-bootstrap and covariance compositional zero-sum. It emits one source-evidence JSON, 2 complete taxonomic-support CSVs and 2 species-by-nine-hue-by-seven-climate data tables **without precise coordinates, original photo IDs or observer identities**.

This is an additive PR. The parent source photo data, models, original P1–P4 conclusions and all earlier sensitivity test results are left unchanged.
