# FCP: do continuous photo-colour and climate associations persist WITHIN geographic regions?

**24 October 2026 — independent spatial-conditional falsification.** This is a source-locked photographic study, not molecular colour-morph/selection inference.

## Why another test is essential

PR #265 measured a small +0.005278 signed white-palette/photo-BIO5 covariance under calendar-month and global *linear XYZ coordinate* controls on the same 1,578 multicolour-eligible taxa. About 75.8% of the signed photographed white-palette covariance arose between the four original PHOTO categories (not genetic morph states); the remaining 24.2% was within those labels. Original photo source 447,500 rows / 4,475 named taxa / 227,107 originally classifiable. All of these results are descriptive, with existing flower-minus-background validation having **failed** in a different original species-disjoint cohort (rho −0.00328, p=0.82).

A geographically varying environmental predictor can track species range composition or camera/photographer/photo-availability differences without producing any true within-region flower phenotype response. Global XYZ residualisation can miss regional jumps. **The next ecological discrimination is whether climatic association exists *inside* actual regional sampling blocks.**

## Test: region + month fixed-effects, then within-block-month matching

Reopen exact original data by SHA:
- original photo source compressed CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6;
- original WorldClim 2.1 BIO1/4/5/12/14/15 10-arcminute archived rasters SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf;
- original 12-month SRAD archived rasters SHA256 c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.

Maintain all 4,475 species in each audit ledger (no hidden zero imputation), with original 227,107 4-colour-classifiable nine-palette photo hues. In each named species separately:

1. Require ≥40 source classifiable photos, ≥20 distinct original observers, ≥5 occupied climate raster grid cells, original valid coordinates (reported accuracy at most 5 km), observed calendar month and complete seven-climate normal vector. Repeat all analyses using colour-blind one original photo per observer.
2. Define two a priori nested geographic locality scales as equal-*degree*, unequal-area latitude/longitude grids: **5° and 1°**. The latter has at most 36 WorldClim 10-arcminute cells per degree square away from boundaries. These are NOT 5 km/1 km sampling limits or true breeding-population boundaries.
3. Fit climate photo vectors on two sets of fixed effects **simultaneously**: *source regional block dummy* and *original calendar-month dummy*. Retain the climate residual in the original species climate-standardised units. Report its retained climate-variance fraction per predictor; require at least 5% rather than silently claiming identification under complete climate/geography collinearity.
4. Within each **region × calendar month** source-photo stratum, **permute the remaining climate residual vectors** (not colour labels) 199 times. Preserve the true photographed colour palettes, photographic labels, original sample size/strata, and geographic source locations. No different-calendar-month or different-region photo is mixed in the null. This cannot eliminate within-stratum confounding by camera, photographer or biological phenology.
5. Require at least **12 source photographs** and ≥2 source region/month strata in which there are **at least two distinct 10-arcminute raster cells with actual within-stratum climate variation** for the particular climatic axis. When the data lack local climate contrasts, label the species as *not identifiable*, not unaffected or constant-hued.
6. Exact 9-D photo hue basis for the total signed continuous-white-palette covariance, *between original coarse PHOTO labels* and *within these labels*. A between-label local effect must also show a sampled local region with ≥2 photographs from each of at least two source-photo labels, plus ≥5 original photos for two labels globally.
7. Report scalar **signed photographic white-palette covariance**, nine-palette covariance vector norm, and their 199-permutation regional/month null excess. Unit = within-source-species climate-standardised photo-fraction covariance, not a per-°C regression slope or causal effect.
8. Species equal in the pooled mean; taxonomic genus bootstrap (499 draws); Holm correction across seven climate variables within each signed/norm/label/scale/source-panel endpoint. Do not pretend this controls 2 scales × 2 observer panels × label/norm family multiplicity jointly.

**Primary question:** Does signed positive BIO5–white photo association survive these more conservative geographical and monthly constraints, particularly at 1° and after observer thinning? **Secondary:** BIO1, BIO4, BIO12/14/15, SRAD mean, and original photographed *within-label* hue variation.

## Interpretations determined by the data

- A strong broad 5° original-photo association but no comparable 1° estimate may mean there is too little within-region climate variation or photos, not that true biological effects are absent. Always report *n tested taxa* and *n informative photos*.
- A raw global photo-climate signal disappearing under 5°/1° controls is evidence that global geography or photo-source sampling can explain the earlier observational link (not proof of that single mechanism).
- A small, persistent positive regional BIO5 photo-colour covariance would be interesting as a within-region **photographic association** but still would not rescue prior negative floral background controls or demonstrate climate adaptation or genetic colour morph frequency.
- Source photo labels can vary due to petal hue, nectar guide/blotch, flower age, bud/fruit/sepal coloration, photography/background/cultivation, and taxonomy. This test does not independently reclassify each photographed organ.
- All seven climate predictors are long-term ecological normals rather than temperature at anthesis. Latitude/longitude grid blocks are not equal area; uneven density can cause residual dependence. Species/photographer geographic coverage remains observationally chosen.
- Earlier strict ≤5 km same-year/month distinct-observer photographic null (P=0.300 or 0.275) remains a separate **unsupported** result, and must not be overwritten by the different local-region/pooled-year-month null.

## Provenance and run contract

Code: scripts/analysis/fcp4475_seven_climate_within_region_falsification_20261024.py

CI: .github/workflows/fcp-4475-climate-regional-fixed-effects-20261024.yml

Product: four complete original **4,475-species** no-coordinate/no-photo-ID CSVs (5° and 1°, original and one-photo-per-observer), paired with a source-verified machine summary. No conclusion is authorised until genuine release checksums, synthetic block-confounding positive/negative fixtures, per-axis support checks and full repository source CI pass.

All previously source-verified results from PR #265 remain unchanged and are the benchmark, not redefined by this new post hoc test.
