# Post hoc soil-neighbourhood coverage recovery diagnostic — 2026-10-08

## Why this route is necessary
The frozen complete-case soil analysis requires >=200 evaluable species per cohort but retains only discovery 172/369, validation 173/363, and third 180/377. Soil completeness is ~68% of D-eligible photo records while macroclimate and elevation completeness is >99.9%. The retained species are far less colour-diverse than excluded species (median D discovery 0.127 vs 0.271, validation 0.109 vs 0.254, third 0.117 vs 0.279). These effects prohibit interpreting the complete-case null as a general soil-null result.

ISRIC SoilGrids 5-km is a coarsened prediction masked in portions of the land surface, including urban, bare, water and ice classes. Missing pixels must not be replaced without explicitly labelling the spatial substitution.

## New diagnostic, defined after the first coverage HOLD
This **separate exploratory sensitivity** asks whether the same fixed soil feature set can be represented using nearby modelled soil predictions, without choosing features or sources in response to the colour findings.

1. Recover the same official, checksum-verified SoilGrids 5-km mean layers: 10 raw properties x three 0–30 cm depths.
2. Construct a raster-level common-valid-cell mask across all 30 inputs. A cell is valid only if all inputs are finite and non-nodata, and the derived thickness-weighted water content at 33 kPa is not smaller than at 1500 kPa.
3. Preserve each photo's true geographic coordinates for all flower-colour and geographic calculations.
4. For photos whose soil location does not lie in a common-valid cell, look for the nearest **common-valid prediction cell** with center at most 10.0 geodesic km from the real photo location, searching a fixed +/-8-pixel window.
5. Never extrapolate past 10.0 km. Never use flower colour, species identity, D, or IBD/IBE results to choose the replacement cell.
6. Where no such cell is available, keep soil missing.
7. Report how many species would meet >=40 valid photos and whether each cohort recovers the frozen >=200-species *coverage threshold*. Do not relax either cutoff.
8. Record the number of substitutions, their geodesic distances, and by-cohort D and absolute latitude for recovered/remaining-excluded species.

## Inference boundary
- This run is a **coverage feasibility** analysis only. It cannot support soil selection, local adaptation, or any IBE-like effect.
- Even if coverage is restored, soil at a neighbouring model grid cell is **a spatial proxy**, not the soil at the photographed plant.
- A subsequent, separately identified model must repeat the same five-component analysis on proxy-supported rows and compare the macroclimate-only benchmark on the **same set of species**.
- The original complete-case route and its coverage HOLD remain unchanged.
