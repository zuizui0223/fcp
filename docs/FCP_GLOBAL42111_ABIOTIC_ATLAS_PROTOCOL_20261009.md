# FCP global 42,111 species: flower-colour biogeography, climate and soils (2026-10-09)

## New priority and separation from 1,499 deep-photo series

Target all **42,111 originally discovered taxa** from the completed 2026-09 flower-ROI measurement set. Do NOT reopen or reclassify image pixels. This is the GLOBAL research mainline; the retrospective 1,499-species high-depth dataset remains a higher-resolution **nested subset** used for within-species checks, not the global population denominator.

Three already-measured denominators cannot be conflated:

1. **Global breadth species**: 42,111 one-photo-per-species outcomes (18,457 classified four-state, 23,654 unclassifiable) for species-equal one-observation colour geography and source availability.
2. **Geographic taxon-cell anchors**: 85,337 distinct species×162 equal-area cell records, 39,075 classified, 128 occupied cells. One taxon can appear in multiple cells, so report species vs cell weights distinctly.
3. **Within-species cross-cell pairs**: 13,416 fixed observer-disjoint, species-unique cross-cell comparisons; 3,196 doubly classifiable, 797 coarse-colour discordant. These provide useful controls against species turnover but are severely missingness-limited; never call 10,220 unclassifiable pairs monochromatic.

## Ecological questions (ranked)

**A. Global species composition:** how does the *conditional colour composition of one classified original flower photograph per taxon* differ across continents, macroclimate and edaphic provinces? All 42,111 remain in the denominator; classification probability is an explicit outcome, and each classified species receives equal weight. A known photographic colour is not a genetically fixed species morph.

**B. Biogeographic maps:** build both observation-richness and species-balanced colour maps from the 85,337 taxon-cell anchors, with classifiability coverage and geographic effort. Report unclassified and uncertain observations; do not convert them to white or monomorphic. Coarse equal-area cell centres describe cells, not organism locations.

**C. Within-species environmental contrast:** paired same-species observer-disjoint cells could test whether larger within-species climate/soil differences are associated with photo-colour discordance *above* geographic separation. Use true observed per-photo coordinates where available, otherwise hold site-level abiotic analysis. Include geographic distance, genus/clade clustering and missing-pair sensitivity. With one fixed pair per species, results are associative only, not fitted local genotype-environment slopes.

**D. Climate and soils:** freeze the existing FCP source convention: WorldClim 2.1 BIO1–BIO19 and monthly mean radiation (plus elevation); SoilGrids 5 km 0–30 cm weighted pH, organic C, N, CEC, bulk density, coarse fragments, sand/clay, plant-available water proxy. Prefer a small interpretable hypothesis set (temperature extremes, precipitation seasonality, aridity proxies, soil pH/N/water retention) before high-dimensional PCA. SoilGrids has substantial land-mask missingness in the old 1,499-series (~68% soil complete; only 172/369, 173/363, 180/377 species in full gradients). Never interpret complete-case nulls as global nulls.

## Necessary first executed preflight: spatial source authenticity

The preflight script replays **exact SHA256-verified old measured source bytes**, inspects their real geographic column schema and valid latitude/longitude coverage, and reports the coordinate join opportunity separately for global breadth, taxon-cell and paired tables. No climate/soil download is justified until the source geometry is known. If the photographed coordinates are absent in the 85,337-cell table, a cell-centroid raster extraction is only a geographically coarse descriptive layer; a sea-centred cell is not a land soil sample.

### Decision sequence

1. PASS measured table and actual coordinate audit. If photo coords absent, locate their original observation-ID/metadata lineage (with matching SHA and IDs) before doing site-scale covariate extraction; do not invent from cell centres.
2. Reuse source-verifiable public WorldClim/SoilGrids layers already used in FCP, not new image/label classifications; extract covariates for the **true** sample location. Record source URL/version, units, depth, layer SHA and footprint; preserve per-feature missingness.
3. For global species-equal conditional colour and taxon-cell geography, fit effort/classifiability-aware associations with covariate coverage diagnostics, spatial blocking and genus-aware resampling. Compare colour structure against taxon turnover instead of presenting one planetary regression as ecological causation.
4. Within species, condition on geographic distance, include all fixed cross-cell pairs in the missingness denominator and report the 3,196 doubly classified pairs separately. Associations with climate/soil are exploratory unless clearly replicated out of sample.
5. Keep the 1,499 original H1/H2 original estimates and the untouched prospectively selected 2,000+730 taxon experiment separate. No redefinition of historical sample denominators or after-seeing-the-result threshold relaxation.

This design is **global ecological characterisation of visible flower colours**, not a test of heritable flower-colour polymorphism across all 42,111 or of adaptive pigment-cost tradeoffs.

### Primary source lineage

- Actual September measured original tables from `2b5390ea74f8d196012499d2d35dd477f9795938`.
- Existing integrated global phenotype/geography PR #135, branch `analysis/fcp-global42111-phenotype-geography-synthesis-20261008`.
- Earlier soil/climate code and source method in PR #126 `scripts/analysis/run_polymorphism_full_gradient_partition_20261007.py`. Its species-selection HOLD is an explicit warning and not a blocker for the broader *location/coverage audit*.
