# FCP 42,111-species world flower-colour × distribution × climate × soil
## Source-verified results and interpretation, 2026-10-09

### Main result in one sentence

Across old source-verified global iNaturalist flower photographs, **climate adds a small, portable predictive signal to between-species visible flower-colour composition beyond geography**, whereas **SoilGrids predictors do not improve the reported held-out species-level model**, and **within-species cross-region differences in photographed colour are not stably predicted by site climate or soil differences**. These describe original source-photograph association, not trait heredity, evolutionary fitness or climate adaptation. This is a *weak prediction* result, not a universal abiotic null or negative proof of local adaptation.

### Immutable source samples

| Estimand/source | Original photos, species or pairs | Classification/sampling | Important distinction |
|---|---:|---:|---|
| Complete world taxon universe | 42,111 named species | 18,457 four-state classifiable | one photo per species, species-equal descriptive sample |
| Regional photo atlas | 85,337 taxon × original equal-area cell | 39,075 four-state classifiable | repeated species across cells, never 85,337 independent species |
| Fixed within-species observer-disjoint cross-cell pairs | 13,416 one-per-species pairs | 3,196 classifiable on both sides; 797 discordant | 10,220 pairs' colour comparison undetermined |
| Original point-coordinate recovery | 100,543 distinct photographed observations across the first two tables | 100,348 real public points | 195 unavailable/changed/obscured or otherwise excluded, no centroid substitute |

Original photograph colours and IDs were already measured in September 2026 at original source commit `2b5390ea74f8d196012499d2d35dd477f9795938`. Location recovery was [GitHub Actions 37855900891](https://github.com/zuizui0223/fcp/actions/runs/37855900891), all 8 source metadata shards and final seal passed. None of these are fresh out-of-sample biological findings: they overlap source historical atlas and original FCP.

### Species-equal ecological predictors: geography / plus climate / plus soil

Completed actual WorldClim BIO1/BIO5/BIO12/BIO15/elevation and 0–30cm modelled SoilGrids pH/SOC/N/clay/available water were attached to all original geolocated photos in [Actions 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446).

- 42,038/42,111 original one-photo taxon species geolocated.
- 18,413 of 18,457 classifiable species have the four WorldClim climate values.
- **14,136 classifiable species** have complete selected climate+soil+elevation; **4,321 colour-classifiable taxa are outside this full complete case**.
- Within the ORIGINAL 42,111 source denominator, 23,654 photographs were unclassifiable and never reassigned to a flower colour.

Fivefold held-out multinomial log loss, lower is better; exactly the SAME 14,136 species in every compared model:

| Held-out groups | Geography/elevation | Add 4 climate covariates | Add 5 SoilGrids covariates | Improvement due to climate | Incremental improvement due to soil |
|---|---:|---:|---:|---:|---:|
| Genus | 1.1782016 | 1.1743888 | 1.1745063 | +0.0038128 | −0.0001176 |
| 162-cell source-photo site grouping | 1.1772831 | 1.1731262 | 1.1739430 | +0.0041569 | −0.0008168 |

The climate increment is small (approximately 0.3% of the log loss); soil does not robustly add incremental predictive value **for this particular four-colour logistic model / exact variable set / source complete-case population**. There is no causal inference and no claim all soil variables are biologically irrelevant.

The complete-case original source colours are: white 6,532; yellow/orange 4,788; red/pink 1,612; blue/purple 1,204. Soil missingness is observational selection; do not extrapolate the conditional complete-case fraction to the entire original 42,111 dataset.

### Within-species photographed cross-cell pairs: geographic distance vs climate vs soil gradients

Source-matched one-original-pair-per-species tests are [verified](https://github.com/zuizui0223/fcp/actions/runs/37863081465); [checked climate-only source JSON](../results/fcp_global42111_within_species_climate_only_pairs_20261009/result.json), and [soil-complete original-pair JSON](../results/fcp_global42111_within_species_abiotic_pairs_20261009/result.json).

| Fixed paired estimator | n both-colour classified and environment-eligible pairs | Discordant original pair colours | Genus-held-out climate gain | Spatial-midpoint-cell-held-out climate gain |
|---|---:|---:|---:|---:|
| Climate-only pair analysis | **3,182** | **793** | −0.000320 | −0.000574 |
| Climate+soil complete pairs | **1,653** | **419** | −0.001178 | −0.000736 |

For soil-complete pairs, adding SoilGrids soil-differences beyond geographic+climatic differences gave +0.000704 under genus-held-out but −0.002155 under spatial-held-out validation. The group-bootstrap confidence intervals of all these increments include 0, and point-signs are not consistent across the two independent grouping policies.

**Ecological interpretation:** the photo-based worldwide colour composition is weakly environmentally patterned across sampled taxa, but direct cross-place same-species pair evidence does not reproduce a similar stable abiotic predictive improvement. This contrast is compatible with species turnover and/or other uncontrolled sources of colour and detection differences; it does **not** isolate species turnover as the only cause, nor disprove substantial within-species genetic or ecological selection.

The climate-only 3,182 pair audit is crucial: the lack of signal in the 1,653 soil-complete pairs is not *solely* explained by soil masking. However, just one photo per endpoint, observation noise, coarse four-state classes and highly selective pair classification severely limit sensitivity. All 13,416 original pairs remain in the denominator, with 10,220 classifiability missing (not coded identical).

### World original 162-cell map with environmental descriptions

Map, all 85,337 original source records and global 4-colour composition are [source-verifiably generated](https://github.com/zuizui0223/fcp/actions/runs/37863503218) and stored in [the 162-cell atlas folder](../results/fcp_global42111_cell_colour_climate_soil_atlas_20261009/result.json).

- Source records 85,337 across **128 occupied of 162** original equal-area geographic cells; 39,075 source photos classifiable four-state.
- Classifiable original taxon-cell states: white **18,374**, yellow/orange **13,018**, red/pink **4,343**, blue/purple **3,340**.
- 85,161 source taxon-cell photos have real public positions; 85,073 have climate; 59,224 have the selected soil covariates; 59,222 have all selected climate+soil+elevation, while **27,032** of the original 39,075 colour-classifiable cases are complete in all selected abiotic fields.
- Regional original taxon-cell WHITE proportions among *classifiable photographed cell anchors*: |latitude| 0–30° **46.39%** (21,705 classified), 30–60° **47.00%** (12,500), 60–90° **49.90%** (4,870); red/pink conditional share declines 12.41→10.14→7.84%. These are **regional taxon-cell frequencies with recurrent taxa** and changing classification success (42.75→49.48→52.38%). They **cannot** be interpreted as evolution toward white with latitude without species-matched correction.
- [All 162 cell photo/env summary CSV](../results/fcp_global42111_cell_colour_climate_soil_atlas_20261009/all_162_cells_photo_colour_climate_soil.csv), [three-latitude-band summary](../results/fcp_global42111_cell_colour_climate_soil_atlas_20261009/three_abs_latitude_regions_colour_and_abiotic_coverage.csv), [original equal-area polygon GeoJSON](../results/fcp_global42111_cell_colour_climate_soil_atlas_20261009/original_162_equal_area_photo_colour_map.geojson). Cell geometry is for drawing only; **soil medians derive from actual original photo sites**, not from the cell centroid.

### Claim boundary / publication

Safe current macroecological claim: *The sampled global geography of visible photographed flower colours shows modest climatic predictability beyond coarse geography at the species-composition scale, while this signal does not transport cleanly to single-pair within-species colour discordance. Added modelled soil gradients show no robust independent improvement in these fixed comparisons.*

Do **not** claim globally genetically maintained flower-colour polymorphisms, phenotypic plasticity, a universal soil zero, causal climatic selection, balancing selection, net flower fitness, or a phylogenetically universal adaptive mechanism. Prior New Phytologist H1/H2 results for 1,499 deeply sampled taxa and the unopened future 2,000+730 confirmation roles remain unchanged.

The strongest potential next biological test needs multiple genetically validated colours per species/population with repeatable location, genotype, pollen success, pigment/fitness data; 42,111 one-photo taxa alone cannot supply that causal evidence.
