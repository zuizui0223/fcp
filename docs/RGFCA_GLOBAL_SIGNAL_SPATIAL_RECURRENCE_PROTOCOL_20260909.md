# RGFCA global signal spatial recurrence — exploratory protocol

Date frozen: 2026-09-09 JST

Upstream recurrence result is retained in `RGFCA_GLOBAL_SIGNAL_RECURRENCE_RESULT_20260909.md`. This stage asks **where recurrent within-species flower–background signal modes repeatedly occur in the same direction**, without using island, climate, latitude bands, realms, pollinators, or other ecological labels to define the map.

## Modes

Map reference modes M1–M3 from the retained recurrence decomposition. These three modes are chosen because the upstream stage, completed before this protocol, showed cosine >=0.90 in all 200 recurrence realizations. M4–M5 remain available as secondary outputs but are not used to define the main spatial pattern because their individual-mode recurrence is weaker.

## Realizations

Reuse the exact upstream realization design:

- 200 realizations, seed `20260909`;
- eligible photos and 12-anchor flower-minus-background vectors unchanged;
- 300 km EPSG:6933 local-state grid with random x/y shift;
- local state >=3 eligible photos and >=2 observers;
- observer bootstrap with replacement plus independent within-observer photo bootstrap;
- species admitted only with >=2 retained local states separated by >=500 km;
- species-centred state vectors;
- each realization's first-five modes matched/sign-aligned to the retained reference modes.

## Fixed reporting lattice

After state scores are computed, assign each state centroid to a **fixed 500 km EPSG:6933 reporting cell** with origin (0,0). The reporting lattice is not shifted.

Within a realization and reporting cell, first average multiple states from the same species, then average those species values equally. A cell-mode value is admitted only when at least **10 distinct species** contribute in that realization.

The 500 km / 10-species rule was selected from coordinate-only support census before spatial mode scores were opened. In the zero-shift support census it leaves 49 reporting cells.

## Spatial recurrence quantities

For each reporting cell and mode report across supported realizations:

1. number and fraction of realizations meeting the >=10-species support rule;
2. median species-equal cell score;
3. 2.5% and 97.5% score quantiles;
4. `P_positive`, fraction of supported realizations with score > 0;
5. directional sign recurrence `R_sign = max(P_positive, 1-P_positive)`;
6. recurrent direction (`positive` if P_positive > 0.5, otherwise `negative`).

A main-map cell additionally requires support in >=80% of the 200 realizations. No score magnitude or sign threshold is used to choose cells.

## Species-conditioned spatial null

For every realization and admitted species, independently permute that species' state score vectors across its retained state centroids before fixed-reporting-cell aggregation. This preserves:

- the species' score values;
- number of states;
- spatial occupancy and support;
- the global mode basis;

while breaking the within-species association between signal state and geography.

Use one independently seeded permutation per species per realization from null seed `20260910`.

For every cell-mode compute the same null recurrence summaries. Report the observed-minus-null difference in `R_sign`. This null is exploratory calibration, not a family-wise-error-controlled significance test.

## Ecological hard stop

No island/mainland status, island area, isolation, climate, latitude, realm, pollinator, reproductive system, or taxonomic composition variable may be used to choose cells, modes, directions, or thresholds in this stage. Ecological overlays are a later stage only after this spatial recurrence result is retained.

## Claim boundary

The output may identify geographically recurrent **photo-derived within-species flower–background signal directions**. It does not establish a biological boundary, animal-perceived contrast, attraction, adaptation, or a causal environmental driver.
