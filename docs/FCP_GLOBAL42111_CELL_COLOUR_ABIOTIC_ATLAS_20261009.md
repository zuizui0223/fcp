# FCP 42,111 species / 85,337 original taxon-cell photos: full 162-cell colour/climate/soil atlas

The already completed historical photographic atlas measured 42,111 taxon-original breadth images (18,457 four-colour-classifiable) and separately 85,337 unique `species × equal-area cell` photo outcomes (39,075 classifiable, across 128 occupied cells). The source has **162 fixed equal-area world cells**, 18 × 9 bins in longitude × sine latitude.

The environmental attachment is already DONE and independently source-verified ([Actions 37861582446](https://github.com/zuizui0223/fcp/actions/runs/37861582446)): exact original public photo locations were recovered for 100,348 of 100,543 unique source image identities, including 85,161 of 85,337 taxon-cell observations. WorldClim 2.1 BIO1/BIO5/BIO12/BIO15+elevation and official SoilGrids 5km topsoil pH/SOC/N/clay/available-water covariates were sampled from these TRUE photo locations and not from arbitrary centroid coordinates.

## Global geographic atlas estimand

The independent geographic script `scripts/analysis/build_fcp_global42111_cell_colour_abiotic_atlas_20261009.py` summarizes the 85,337 already classified or unclassifiable original taxon-cell observations, retaining empty and undetermined cells. Export exactly:
- `all_162_cells_photo_colour_climate_soil.csv`: every original equal-area cell's sample count, classifiable/missing photo count, conditional counts and fractions for four coarse colours, original-source latitude/longitude cell-centre *display coordinate*, original true public point/soil/climate coverage, and median environmental covariates over the source photo locations in each cell (NOT cell centres);
- `three_abs_latitude_regions_colour_and_abiotic_coverage.csv`: original taxon-cell totals, photo-classification denominator, each colour fraction, soil/classification coverage in fixed |latitude| bands 0–30, 30–60, 60–90°;
- `original_162_equal_area_photo_colour_map.geojson`: original equal-area grid polygon geometry and descriptive photographed source fields, suitable for world map; polygon midpoint is **rendering only**, not an organism soil/geographic measurement;
- `result.json`: explicit original numerator and denominator receipts.

A photo equal-cell colour fraction is NOT a global equal-species prevalence (widespread species recur across geographic cells), a within-species colour frequency, a genetically proven natural polymorphism, a climate-causal effect or a reproductive-fitness contrast. A species photographed as red in one cell and white in another is an observational floral-colour pair, not an observed allele replacement.

## Results and interpretation boundaries

Separate from this mapping objective, the global one-photo-per-species species-equal logistic model has 14,136 fully geolocated/climate/soil/photo-classifiable old taxa. Adding WorldClim climate to geography slightly improved out-of-genus and out-of-cell four-colour prediction (heldout log-loss improvements +0.00381 and +0.00416); adding the five SoilGrids proxies did not (+-0.00012 / -0.00082). The species-level model is subject to species turnover and photographic sampling bias.

The completed **within-species cross-cell fixed pair** analysis found 3,196 doubly colour-classifiable pairs among 13,416 original species-pairs. Only 1,653 had climate+soil at both sites (419 mismatched photographs), and no climate/soil incremental out-of-group predictive advantage with interval excluding zero. Crucially a second soil-UNRESTRICTED climate-only check retained 3,182 pairs (793 photo-colour mismatches), and climate still showed no stable positive cross-genus or cross-region predictive gain (−0.00032 / −0.00057; both intervals cross zero). See the corresponding separate result receipts. These do **not** prove zero climate selection: one photographed pair per taxon and coarse colour categories have limited sensitivity.

This 162-cell atlas is descriptive geography: unlike site-level input prediction, it must retain all 85,337 observations and never report per-cell phenotype percentages without cell-specific classification fractions and denominators. It does not alter the original 1,499 high-depth cohort, New Phytologist H1/H2 analyses or the untouched independently preselected 2,000+730 species.

Reproducibility: [original full abiotic Actions](https://github.com/zuizui0223/fcp/actions/runs/37861582446), original fixed photographed pair study [Actions](https://github.com/zuizui0223/fcp/actions/runs/37863081465), and dedicated 162-cell atlas GitHub Actions in this PR.
