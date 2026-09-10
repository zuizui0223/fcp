# RGFCA Step 7B — minimum global atlas rarefaction

Date: 2026-09-10 JST
Status: frozen before opening rarefaction outcomes.

## Question

How few species and raw candidate photos per species are needed to recover the repeatable global flower-colour pattern present in the full frozen RGFCA photo frames?

## Independent frames

Run discovery and reserve separately. A design is called globally robust only if it passes the same predeclared recovery rule in both species-disjoint frames.

## Sampling grid

Species counts: **25, 50, 100, 150, 250, 350**.
Raw photos/species: **10, 20, 40, 60, 80, 100**.
Realizations: **200** per design per tranche.

Each realization independently permutes species and the frozen raw photo rows within species. Raw rows are sampled before the `global_classifiable` gate; downstream estimators use only admitted four-state rows. Thus classification yield is part of the minimum-data problem rather than conditioned away.

Only species with at least 100 frozen raw rows enter the rarefaction universe.

## Full-data reference atlas

Use all eligible species and the first 100 frozen raw rows/species (equivalent to all 100 when exactly 100 exist).

Spatial grid: fixed 10-degree latitude x 10-degree longitude cells. Within each species x cell, calculate the four admitted morph frequencies. Then average those frequency vectors equally across species represented in the cell. Reference atlas cells require >=3 represented species.

This produces a species-equal global atlas rather than a photo-density-weighted atlas.

## Recovery metrics

For every realization/design:

1. **atlas cell coverage** = fraction of full-reference cells represented by >=1 sampled species with >=1 classifiable photo in that cell;
2. **atlas composition cosine** = cosine similarity of the flattened four-morph vectors over cells shared with the full reference, with cells weighted equally;
3. **D-distribution quantile MAE** = mean absolute difference between sampled and full-reference Simpson-D quantiles at q=0.1,...,0.9. Sampled D is calculated for species with >=5 classifiable photos in that realization.

Also report, but do not use for the primary minimum-data gate:

- median number of classifiable photos/species;
- fraction of sampled species with >=5 classifiable photos;
- mean/median D;
- recovered second-morph >=10% fraction.

## Predeclared pass rule

A realization passes when all are true:

- atlas cell coverage >= **0.80**;
- atlas composition cosine >= **0.90**;
- D-distribution quantile MAE <= **0.05**.

A design passes a tranche when >= **90% of 200 realizations** pass.

The **minimum robust design** is the passing design with the smallest nominal raw-photo budget `species_count * raw_photos_per_species`; ties are broken by fewer species, then fewer photos/species.

A cross-tranche minimum is declared only if the same nominal design passes independently in both discovery and reserve. If no design passes both, report `NO_COMMON_MINIMUM_WITHIN_GRID`; thresholds are not relaxed.

## Scope

Step 7B establishes data sufficiency and repeated atlas recovery. It does not test causes of C*/S* or old 34-species ecological hypotheses. Those move to Step 7C only after the minimum-data behavior is known.