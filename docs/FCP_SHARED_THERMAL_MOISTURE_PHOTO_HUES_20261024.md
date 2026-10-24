# FCP: shared thermal and moisture climatic directions of nine-hue photo colour

**24 October 2026. Source-locked additive ecological analysis stacked on PR #268.**

## Current scientific question

The original FCP 447,500-photo, 4,475-species data contain a small shared positive direction between warmer observed habitat climate and white-rich photographic palettes. In the joint seven-variable model, BIO5's *individual unique slope* was unidentifiable for many species, because annual mean and maximum warm-season temperature are correlated, as are precipitation components.

The biological question is not whether BIO1 or BIO5 can "beat" the other under an extremely collinear 7-variable regression. It is whether **the shared level of warmth**, treated as a photo-outcome-blind environmental gradient, carries a repeated nine-hue PHOTO colour association once *precipitation, climate seasonality, and solar radiation* are accounted for.

## Fixed outcome-blind macroclimate coordinates

First standardize each of the authentic source seven WorldClim/photo climatic normals within each named species by its photographed SD. Define five variables BEFORE opening colour responses:

- thermal_level = (standardized BIO1 annual mean temperature + standardized BIO5 warmest-month maximum temperature) / 2
- precip_amount = (standardized BIO12 annual precipitation + standardized BIO14 driest-month precipitation) / 2
- BIO4_seasonality = standardized BIO4
- BIO15_seasonality = standardized BIO15
- SRADmean = standardized annual-average WorldClim solar radiation from the exact 12 monthly TIFF rasters.

The first two axes are simple equal-weight source-climate averages; **no colour-outcome PCA/feature selection** is permitted. Temperature seasonality BIO4 and precipitation seasonality BIO15 remain explicit nuisances, because they are not equivalent to the level of warmth or rainfall. Differences between BIO1 and BIO5, and between BIO12 and BIO14, are **not held fixed**; these new estimands are climate GROUP-level associations, not the unique effect of a specific original BIO predictor.

In each species, remove original photo-calendar-month and broad geographic Earth-centred X/Y/Z location effects from all five macro variables. Project EACH climatic macrovariable on the other FOUR. Require its unique remaining climate variance to exceed **5% of that species' original photo-climate macro-axis variance**. Report both numerator and retained fraction. When no climate-domain contrast is identifiable, leave a missing effect; do not infer a neutral/zero climatic hue response.

## Authentic photo source and precise inferential scope

- Source is the exact 447,500 original measured photographs of 4,475 named taxa; the source original classifiable 9-hue photo-ROI palette comprises **227,107** photographs, missing/failed source classes remain uninterpretable.
- Require ≥40 classified photo records with original coordinates/dates and complete WorldClim, ≥20 unique original observers and ≥5 distinct original raster climate cells per species. Repeat on one colour-blind-selected original photo per observer.
- Source original nine-palette channels: white, yellow, orange, red, pink, magenta, purple, blue, bronze. A nine-hue image vector is **photographic ROI colour**, not calibrated fresh petal reflectance.
- For each species/climate domain, compute the covariance between each photographic nine-palette fraction and the climate-domain residual. Subtract the **exact expectation** of permuting whole climate residuals within calendar month, then divide by the residual climate-domain variance. This ratio is a *descriptive partial-slope-like photo association*, **not a causal change or °C-normalized pigment response**.
- The exact PHOTO colour-label covariance identity partitions the original 4-state photographic label **between** and **within** contributions. The between contribution requires two colour labels with ≥5 source photos each. This is NOT a genetic discrete morph versus continuous plasticity distinction.
- A fixed unit white-versus-nonwhite photo-palette axis q=(8,-1,-1,-1,-1,-1,-1,-1,-1,-1)/sqrt(72) measures the mean signed direction; the complementary 8D chromatic-only subspace captures differences among nonwhite hues. A large percent white alignment means only that the **cross-species mean 9-hue covariance vector** follows white more than other hues.
- Use whole-genus cluster bootstrap with **999** draws for species-equal means, seven? Exactly **five macro climate axes**: simultaneous max-T family of five white projections, separately 5×8 chromatic-only hue comparisons and 5×9 source hue components, within each of two source photo panels and original photo-label component. Multiple families are exploratory and NOT controlled jointly.
- Report original five species-disjoint sampling cohort white slope signs as descriptive replication, not an untouched novel test; all photo cohorts share the iNaturalist photo source and measurement assay.
- Pair temperature and precipitation alignment within **exactly identical species and source-photo eligibility**, and cluster bootstrap their difference; do not compare white-axis percentages computed on different species supports or interpret them as fractions of molecular colour genes.
- The entire original 4,475-species ledger must remain in the output despite non-estimable climate/photo combinations. No original source photo, source classification, original FCP P1–P4 paper results or earlier climate-sensitivity results are altered.

## Pre-specified interpretation limits

**Possible positive result:** A same-signed photographic white/BIO1+BIO5 shared-temperature gradient persisting after moisture amount, seasonality and radiation controls would support a *climatically associated shared PHOTO colour geometry* more identifiable than a unique BIO5 effect in an overcollinear seven-variable model.

**Possible non-support:** Even the joint thermal domain may lack independent climate variation, or the photo hue-temperature direction may not replicate. This would bound the original source-photo observation, not prove absence of real flower biology.

**Nonclaims in all outcomes:** no independently botanist-certified open petals, genotype-linked flower colour polymorphism, native-versus-cultivated systematic disambiguation, photo background correction, selection or causal thermal pigment response. Earlier flower-minus-background validation was unsupported (mean partial rho −0.00328, p=0.82) and strict ≤5-km independent-observer local photo analyses were unsupported; neither is erased by this different photo/climate analysis. Long-term WorldClim site normals do not measure temperature at anthesis.

The story remains about **intraspecific differences in colour across species**, not mechanisms inferred only from a data-driven regression. Should a positive thermal-group finding survive, testing inherited flower colours and fitness remains a separate future task.

## Reproduction

Source code: scripts/analysis/fcp4475_shared_thermal_moisture_photo_hues_20261024.py.

CI with source archive SHA256 checks, exact synthetic source covariance, original photo and rasters, all 4,475 taxon ledgers and numerical summary: .github/workflows/fcp-4475-shared-thermal-moisture-photo-hues-20261024.yml.

Three output types: an original source ecological JSON with five-domain data, two taxon-level observer-thinned/full original support CSVs, and two source taxon-by-domain-by-photographic-colour component coefficient NPZs. Each excludes original individual photo identifiers, precise raw coordinates and individual observer IDs.
