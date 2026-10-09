# Global flower-colour and climate across 42,111 species: tested ecological scale

## Completed results (2026-10-09)

This manuscript-scale **photograph-based macroecological** synthesis keeps three nested inference domains that must never be mistaken for the same statistic:

1. **Between sampled species worldwide:** original 42,111 iNaturalist-discovered flower taxa, 18,457 with one four-state classifiable flower photo; 14,136 original species with one-photo colour, real-site geographic, WorldClim and selected SoilGrids values. Adding 4 climate features improves held-out multiclass log loss after coarse geography by ~0.0038–0.0042 (genus and photo-cell validation). Soil addition does not help robustly.
2. **Between different species WITHIN genus:** source climate-complete 18,413 nominal species, 13,307 source photographed species in 1,417 training-known genera testable under held-out original source geographic cells. Same real photo, no reclassification; training-only genus colour means. Unique full climate block improvement **+0.001465** species-equal Brier, original cell-block 95% interval **[+0.001005,+0.001994]**. Soil-complete sample 14,136 taxa, 9,498 evaluated among 1,122 genera: climate improvement **+0.001201**, interval **[+0.000623,+0.001843]**.
3. **Between photographed localities WITHIN species:** one fixed observer-disjoint cross-cell pair per 13,416 species, 3,196 both endpoints flower-colour-classifiable, 797 source mismatches; 3,182 climate-complete tested pairs, 793 mismatches. No stable independent improvement in photographed colour mismatch prediction when climate added to geographic distance, latitude and elevation.

## Thermal vs rainfall/moisture block attribution

BIO1 = mean annual temperature, BIO5 = maximum temperature of warmest month; BIO12 = annual precipitation, BIO15 = precipitation seasonality. Four environmental predictors are correlated; unique removal/addition scores are noncausal partial predictive associations.

| Scale and fixed original photo sample | Thermal unique gain beyond geography+moisture | Moisture unique gain beyond geography+thermal | Confidence status |
|---|---:|---:|---|
| Within genus, 13,307 evaluated different species, spatial photo-cell holdout | +0.000235 (95% original cell bootstrap −0.0000003, +0.000470) | **+0.001367** (+0.000882, +0.001865) | Moisture predictive increment positive under cell resampling; temperature cell CI includes 0 |
| Within genus, 9,498 evaluated different species from soil-complete source, spatial photo-cell holdout | +0.000200 (−0.000055, +0.000452) | **+0.000935** (+0.000503,+0.001404) | Moisture replicates as descriptive conditional association in second nested completeness set |
| Within one species, 3,182 original pairs (793 photo colour mismatches), genus-heldout | +0.000190 (−0.000780,+0.001197) | −0.000854 (−0.001965,+0.000253) | Neither thermal nor moisture gives robust mismatch prediction |
| Within one species, same 3,182 original pairs, original pair-midpoint geographic cell-heldout | −0.000448 (−0.001621,+0.000701) | −0.000155 (−0.001038,+0.000740) | Neither supports unique climatic mismatch prediction |

Within-genus models use **multiclass Brier improvement**, while within-species pairs predict photo-colour mismatch with **binary log-loss improvement**. Their numeric magnitudes are NOT directly comparable; only their respective within-design signs, intervals and source-population scopes can be contrasted. Confidence intervals are conditional resampling of existing out-of-fold errors, NOT independent evolutionary replication.

## Genus-size sensitivity

Species-equal within-genus climate gain persists if each evaluated genus is given equal aggregate weight: **+0.003291** for 1,417 genera (genus-cluster bootstrap [0.002355,0.004346]), and **+0.002250** for 1,122 soil-complete genera (genus-cluster bootstrap [0.001249,0.003335]). But approximately **half** the genera show a positive individual mean prediction gain: different lineages do not respond universally.

Soil increment remains unstable: original cell-block species-equal climate+soil Brier gain −0.000340, [−0.000990,+0.000242]; even though genus-equal mean is +0.000811, the **source-cell-block genus-equal** interval [−0.000082,+0.001677] spans zero.

## Strongest current source-limited ecological finding

*The sampled geographical composition of flower colours exhibits a small precipitation-associated signal among different species within genera, conditional on fixed source photographs, genus baselines and spatial holdouts. The same environmental signal does not transport to single photographed within-species geographical colour contrasts.* This is compatible with differences between ecological/phylogenetic taxonomic scales, but **does not prove that climate selected flower colours during speciation, or that same-species local adaptation is absent**.

An unmeasured critical alternative remains: geographic climate may also affect whether citizen-science photos are **flower-ROI-classifiable**. Original 23,654/42,111 source species lacked a classifiable four-state colour; this massive missingness must be explicitly audited before any evolutionary mechanism assertion.

## Verified machine results

- [Within-genus photo colours, genus-equal and moisture/heat source result](../results/fcp_global42111_within_genus_photo_environment_20261009/result.json), source verified by [Actions 37888110321](https://github.com/zuizui0223/fcp/actions/runs/37888110321).
- [Original 3,182 same-species photo-pair moisture and heat block falsification](../results/fcp_global42111_same_species_moisture_blocks_20261009/result.json), source verified by [Actions 37888370629](https://github.com/zuizui0223/fcp/actions/runs/37888370629).
- [Global 42,111 sampled species and original 85,337 geographic taxon-cell ecological atlas PR #138](https://github.com/zuizui0223/fcp/pull/138).

All analyses are retrospective and retain the old FCP photograph image classifications, source photo IDs and original environmental rasters. The original New Phytologist 1,499 high-depth H1/H2 analyses and unexposed preselected 2,000+730 study remain unchanged.


## Direct negative control: precipitation vs source PHOTO CLASSIFIABILITY

Completed original photo ascertainment diagnosis from [source-verified Actions 37888767069](https://github.com/zuizui0223/fcp/actions/runs/37888767069), 6 synthetic tests, and [machine result](../results/fcp_global42111_photo_classifiability_moisture_control_20261009/result.json). The binary response is **original flower ROI four-state CLASSIFIABLE or not**, NOT whether the biological species' flower is white/pigmented.

- Source-wide original photo denominator preserved: **42,111** nominal taxa, **18,457** photo-classifiable, **23,654** photo-unclassifiable.
- **42,014** source photos have climate+altitude available (18,413 classified, 23,601 unclassified); **33,810** train-genus-supported heldout photo opportunities. Unique BIO12/BIO15 precipitation increment beyond genus+geography+BIO1/BIO5: **+0.0000363 binary Brier**, 95% fixed geographic-cell bootstrap **[−0.0000488,+0.0001274]** — NOT supported. Unique BIO1/BIO5 thermal increment: **+0.0001728**, 95% **[+0.0000419,+0.0002999]** — very small positive classification-success predictability.
- **32,231** original source photos with full soil+climate+altitude available (14,136 classified, 18,095 unclassified); **24,731** in train-genus-supported spatially held-out photos. Unique rainfall-classifiability increment **+0.0000543**, CI **[−0.0000735,+0.0001729]** — NOT supported. Soil increment on classifiability also spans zero.

**Interpretation:** The source photograph-classifiability negative control does **not** show an independent rainfall classification-success signal analogous to the small rainfall/photo-colour association among classified *different congeneric species*. This weakens the simple explanation that wet places merely yield more classifiable flower photos. However, the same assay cannot test **colour-dependent unclassifiability** (true flower colour is unknowable for original 23,654 unsuccessful images), residual geography/taxonomy/photography confounding, gene flow, pigment physiology, or environmental selection. Binary classification-success Brier improvements are not numerically comparable to multiclass colour Brier improvements.

The supported ecological claim remains **retrospective, source-photographic and scale-conditional**, not experimental: a subtle precipitation-associated visible flower-colour signal among sampled different species within a genus, no robust precipitation predictor of single-pair within-species spatial photo-colour mismatch, and no simple independent precipitation prediction of overall ROI classification success. No claim that climatic selection occurred during speciation is warranted.
