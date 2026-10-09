# FCP global 42,111 flower-colour atlas: same-species moisture vs temperature photo contrast

## Reason for independent scale replication

The original species-equal climate model and the **genus-conditional BETWEEN-SPECIES** colour prediction found a small, spatially held-out climate association. [Verified within-genus block ablation](../results/fcp_global42111_within_genus_photo_environment_20261009/result.json) pointed mainly to the BIO12 annual precipitation + BIO15 precipitation seasonality block, not to unique BIO1/BIO5 temperature information, after holding genus means and source geographic position.

This does NOT entail that photographed flower colour shifts with annual precipitation *within a single plant species*. The original 13,416 species-unique observer-disjoint fixed across-cell pairs can test that separate scale.

## Exact historical same-species photographed comparison

- Keep original 13,416 source pairs; 3,196 classifiable at BOTH endpoints; 797 observed coarse-colour discordant, and 10,220 with at least one unclassifiable image. Never assume unclassified means same/monomorphic.
- Use the already verified *original photo ID / observation ID* joins into the 100,543-source-photo WorldClim and soil geographic data. Do not regenerate image pixels, source observation IDs, photo-colour labels, public coordinates or rasters.
- Use exactly the previously confirmed **3,182 original two-site climate-complete photographed same-species pairs**, of which **793** have discordant photo-colour states; do **NOT** select only 1,653 soil-complete pairs. Source taxon, photo identity, observer separation, region-cell and pair state remain unchanged.
- Outcome: original source `discordant` vs `same` photographed four-colour classes. This is NOT within-population genetic polymorphism or local adaptation.
- Fix fivefold genus-heldout and true original pair-midpoint-cell-heldout predictions of logistic discordance probabilities using identical source pairs/holdout folds and 4 predictor families: (A) geographic distance + absolute latitude difference + elevation difference, (B) + BIO1/BIO5 absolute temperature differences, (C) + BIO12/BIO15 absolute moisture differences, (D) + both climatic blocks.
- Infer unique temperature increment by comparing C vs D, unique precipitation/moisture increment by B vs D, and total climate increment by A vs D. The original BIO1/BIO5 and BIO12/BIO15 are CORRELATED, so unique block contributions do not have to add to total.
- Save original-group block-bootstrap percentile intervals (999 resamples) of already held-out **log-loss differences**, conditional on the fixed model and fold assignment. The test is POST-OUTCOME exploration, not independent preregistration. Do not choose other variables after inspecting the result.

## Interpretation ladder

If precipitation/moisture uniquely improves the within-genus BETWEEN-SPECIES result but not photographed within-species cross-region mismatch, this is evidence for **different predictive scales in the original photo data**, not proof of evolutionary species replacement as the sole mechanism. Differences in photograph classifiability, coarse flower colours, within-population variation, and the original single pair per species limit sensitivity.

A positive within-species moisture increment would instead challenge the prior conclusion that source photograph climate covariation does not transport within species, but require a separate confirmatory cohort and photographic/genetic reliability work.

Output: `results/fcp_global42111_same_species_moisture_blocks_20261009/result.json`. The unchanged original 1,499 species photo-intensive analysis and untouched prospective 2,000+730 taxa are not reused or relabeled.
