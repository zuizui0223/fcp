# FCP 42,111 global flower colour: same-species cross-cell environmental contrast

## Why the species-equal global result needs a second level

Source-verified [global model result](../results/fcp_global42111_colour_climate_soil_characterization_20261009/result.json): 42,111 original taxa; 18,457 one-photo four-state classifiable; 14,136 complete climate+soil cases. Mean four-state photo-colour cross-validation log loss improved slightly when four WorldClim variables were added to geographic/altitude covariates, but modelled SoilGrids topsoil added no stable out-of-group predictive gain:

| Group held out | Geography only | Geography+climate | Geography+climate+soil | Climate improvement | Soil incremental improvement |
|---|---:|---:|---:|---:|---:|
| genus | 1.178202 | 1.174389 | 1.174506 | +0.003813 | -0.000118 |
| source equal-area cell | 1.177283 | 1.173126 | 1.173943 | +0.004157 | -0.000817 |

These are **associations across different nominal species** based on one original photo per species. Climate's predictive increment may describe the geography of which species occur there, not adaptive change of flower colour within a lineage. Soil-null predictive increment is specific to these five topsoil indicators, logistic linear model, photo ROI classifiability, spatial folds and soil complete cases; it is not the absence of soil effects on petals.

## Nested within-species fixed-pair test

Original September photos comprise **13,416 one-per-species cross-cell pairs** drawn with observer disjointness; **3,196** have classifiable flower-colour states at both endpoints, **797** are coarse four-state mismatches. The remaining 10,220 are outcome missing, not matching/missing polymorphism.

Join each endpoint using the exact `inat_taxon_id, observation_id, photo_id` key to the already source-verified **100,543 unique original-photo environmental records**, never latitude-cell centroids. Do not drop original pairs from the missingness denominator. For both endpoints:

- actual geodesic distance (log1p transformed), unsigned absolute-latitude difference, unsigned elevation difference;
- unsigned differences of WorldClim BIO1, BIO5, BIO12 and BIO15;
- unsigned differences of SoilGrids pH, SOC, N, clay and available water proxy.

The conditional outcome is `discordant` vs `same` original photographed four-state colour, **not** genotype, population coexistence or local phenotype fitness. Only species with **both** original photos classifiable, geolocated and climate+soil complete can enter the predictive model. A minimum of 300 complete species-pairs, original 3,196 two-classifiable and 797 discordant denominator checks, and at least five independent fold groups apply. Other species remain in explicit coverage/missing records.

Models use identical complete species-pairs for comparisons:
- GEOGRAPHY;
- GEOGRAPHY+CLIMATE;
- GEOGRAPHY+CLIMATE+SOIL.

Compare genus-held-out and true geographic great-circle pair-midpoint-cell held-out GroupKFold predictions, log loss, source-discordance AUC and incremental prediction gains. Group bootstrap intervals resample the original grouping of **fixed out-of-fold losses** and do not pretend to integrate model refitting uncertainty; avoid elevating a small improvement to a confirmatory evolutionary discovery.

## Decision

A genuine within-species environmental contrast would be a useful new macroecological observation only if positive improvement appears under both genus- and geographic-blocking, with a compatible uncertainty interval and adequate coverage. The previously measured paired photos cannot establish that petals were genetically polymorphic, or that climate caused adaptive selection.

A below-300, unclassifiable-heavy or opposite-sign result is a **negative/uncertain study-specific predictive outcome** and should not be reframed as a universal biological absence or rescued by changing climate/soil feature sets after inspecting the target.

### Provenance

- Original 13,416 pair file at fixed September commit `2b5390ea74f8d196012499d2d35dd477f9795938`, SHA256 `26755316b2cbd07e1e2241eb4ddc1cddc85219fd7e96424df8c14df32575a4d2`.
- Original WorldClim and SoilGrids extraction, [Actions 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446), SHA-checked source files; 70,211 of 100,543 old photo records with climate+soil+altitude complete, and 14,136 classifiable one-photo species (not 14,136 genetic polymorphisms).
- Output: `results/fcp_global42111_within_species_abiotic_pairs_20261009/result.json` after source-verified Actions succeeds. Run has no image downloads, photo classification or new environment acquisition.
