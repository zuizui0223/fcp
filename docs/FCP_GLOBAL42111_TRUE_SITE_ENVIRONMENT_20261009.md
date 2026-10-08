# FCP 42,111 actual flower-colour measurements and WorldClim / SoilGrids site extraction

## Origin and geographic confirmation (frozen, 2026-10-09)

This retrospective atlas uses **exactly the previously measured 42,111-species one-photo original breadth** (18,457 classifiable) plus the separately measured **85,337 species×equal-area-cell observations** (39,075 classifiable), never the smaller 1,499 high-depth photo-only source as the global denominator.

The original measured images themselves lacked longitude and latitude, but the original observation and photo IDs were intact. The full [validated original-ID recovery workflow](https://github.com/zuizui0223/fcp/actions/runs/37855900891) matched 100,543 unique old observation/photo IDs and recovered **100,348 exact source-consistent, publicly available point coordinates**, with **195 source photo sites unresolved**. Original breadth **42,038/42,111** and taxon-cell **85,161/85,337** received exact publicly permissible point coordinates, preserving original image measurement status and taxon identity. No new image pixels or colour states were generated.

The archived complete full-position source is the immutable workflow artifact `fcp-42111-verified-public-original-geocoordinates-20261009`, run **37855900891**, and the [durable original source position receipt](https://github.com/zuizui0223/fcp/blob/analysis/fcp-global42111-abiotic-atlas-20261009/results/fcp_global42111_all_public_original_coordinates_20261009/result.json).

## Actual abiotic covariate attachment

- **WorldClim 2.1** 10-minute `BIO1` mean annual temperature, `BIO5` maximum temperature of warmest month, `BIO12` annual precipitation and `BIO15` precipitation seasonality, plus `elev`. The official BIO zip must pass the previously pinned original SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`; elevation zip SHA is logged in each run. WorldClim data's stored scaling must be treated as native until checked; do not label raw BIO1 or BIO5 Celsius without conversion.
- **ISRIC SoilGrids** versioned/latest original 5-km aggregated raster `mean` product from three 0–30 cm depths (0–5,5–15,15–30 weighted 5:10:15): soil pH, soil organic carbon, nitrogen, clay and available-water proxy `wv0033 - wv1500`. Download only the specified 18 mean TIFF layers and validate each layer against the official provider's `checksum.sha256.txt`; these are modelled predictions, **not actual root-zone soil assays**.
- Do not replace missing or masked soil pixels with nearby soil, distance matching or cell centroids. Report field-by-field availability for the entire 100,543 unique original photos, all 42,111 species-equal old breadth images, and all 85,337 taxon-cell images. Every original flower-colour missing/unclassifiable image remains in the output, never white or monomorphic by implication.
- Geographic site point must match the identical old species taxon, observation ID and photo ID. All environmental samples with no recovered genuine point coordinate are `NaN`. No new images downloaded, no independent prospective 2,000+730 selected species inspected.

## Global scientific tests that follow in the SAME reproducible job

Use *original species-equal breadth* 42,111 taxa as the universe; 18,457 original colour-classifiable one-photo records form the original observable outcome sample. Within one **fixed, identical, soil+climate complete-case cohort**, compare four-state flower-colour prediction from (GEO) absolute latitude, longitude sine/cosine and elevation, (GEO+CLIMATE) these plus 4 WorldClim indicators, and (GEO+CLIMATE+SOIL) these plus 5 soil indicators. Report genus-held-out and 162-equal-area-cell-held-out fivefold `log_loss` on the identical tested species. Define soil incremental predictive gain as heldout log loss for GEO+CLIMATE minus GEO+CLIMATE+SOIL; a positive number is a predictive improvement, not causal soil selection or true polymorphism.

Also export median and interquartile native-scale climate/soil profiles for each of the four photo-colour classes among complete source-photos only. Separately export **all-original-42,111 species** regional classification and abiotic completeness, including `NO_EXACT_PUBLIC_GEO`.

This is exploratory, outcome-exposed observational biogeography. It is not a population genetics analysis, a new 42,111-species independent confirmation, a causal environmental sorting estimate or pigment-fitness study. Original high-depth 1,499-species H1/H2 and prospective reserved 2,000+730 species are completely unchanged.

## Execution boundaries

The dedicated workflow `.github/workflows/fcp-global42111-true-site-environment-20261009.yml` downloads only the completed source coordinates and the official geospatial raster inputs, runs pinned synthetic tests, verifies all denominators, then writes both full environmental source-joined CSVs and predictive ecology results as separate checked artifacts. It must not re-run the costly original iNaturalist source coordinate retrieval. **No environmental effect is claimed before the full extraction, blocked prediction and workflow validations finish.**
