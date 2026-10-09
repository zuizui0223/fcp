# Global 42,111 original flower photo classifiability: explicit precipitation negative control

## Scientific risk before interpreting the significant photographic colour signal

The already outcome-exposed global original flower-pigment visual atlas successfully classified **18,457 of 42,111 photographed nominal taxa**, while **23,654 taxa** have an original photo with no valid four-state flower-colour outcome. The preceding source-verified ecological results show a small within-genus BETWEEN-SPECIES rainfall (BIO12 annual precipitation, BIO15 precipitation seasonality) prediction increment across 13,307 evaluated congeneric original photo taxa, but no stable precipitation increment in 3,182 same-species photographed two-locality pairs.

**Alternative ascertainment explanation:** precipitation or climate could affect the likelihood that a citizen-science photo has classifiable focal floral pixels, rather than the actual biology of taxon's phenotype. Flower presence, flower orientation, camera/flower ROI segmentation, taxonomic and photographic opportunity all alter source missingness. This is a negative control for the observational measurement pipeline and cannot by itself correct bias.

## Full source denominator; source-verified input only

Use all 42,111 original source one-photo taxa, 18,457 classified and 23,654 unclassified. Retain missing source photo positions as missing (original public photo coordinates were source-verified for 42,038/42,111). Source geographic macroclimate and SoilGrids features were already sampled at the exact photograph site in [Actions 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446); do not rerun original image measurement, public photo metadata or raster acquisition.

The source-only colour classification outcome for THIS assay is **binary photo ROI classifiable vs unclassifiable**, not white/nonwhite, genetic heterogeneity or flower-colour categories. For accurate original-source denominators:
- 42,014 of 42,111 original photos have all four WorldClim predictors and elevation;
- 32,231 of 42,111 have those climate variables plus selected SoilGrids predictions.

Use fivefold geographic original photo-cell heldout GroupKFold for each fixed pool. Estimate training-only genus classifiability means and train-only within-genus centred geographic/environment covariate slopes; include only heldout photos from genera with >=2 different training species. All nested covariate models within a pool must evaluate exactly the SAME original photo IDs, original heldout folds and training-only genus baseline.

Fixed predictor families: genus baseline only, genus+geography/elevation, +two-temperature BIO1/BIO5, +two-moisture BIO12/BIO15, all four climate indicators, then five SoilGrids predictors only for soil-complete pool. Report original heldout binary Brier improvements and 999 original photo geographical-cell bootstrap intervals. These are **post-outcome conditional sensitivity tests**, not causal photography or plant biology mechanisms.

## Interpretation safeguards

- If rainfall/moisture predicts PHOTO CLASSIFIABILITY in addition to classified PHOTO COLOUR, the previous colour association could be selectively observed. The control does not prove it was entirely an artefact.
- If rainfall/moisture does not predict classifiability, this does not exclude other missing-not-at-random sources (plant flower state, photography, morphology, species taxonomic community, palette classifier biases).
- The original genomic/genetic state of the source species, breeding population, pigment chemistry, fitness and local adaptation are all unmeasured.
- The 23,654 missing source photograph labels are never set equal to white, yellow, a monomorphic colour or negative flower presence.

Results are source-verified only after the dedicated [GitHub Actions workflow](../.github/workflows/fcp-global42111-photo-classifiability-control-20261009.yml) succeeds and commits `results/fcp_global42111_photo_classifiability_moisture_control_20261009/result.json`.
