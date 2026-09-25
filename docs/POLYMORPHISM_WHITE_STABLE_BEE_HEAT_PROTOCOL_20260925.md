# Stable local bee versus heat discriminator — 2026-09-25

Status: **joint-analysis design frozen before local bee predictors from the Noori et al. 2026 curated dataset are computed or inspected.**

This is a mechanism follow-up after the white outcome and the prospective BIO5 result are already known. It is therefore not an untouched confirmation of the white phenotype. Its purpose is to discriminate the existing broad-scale heat clue from a stable independent bee-community covariate.

## Sources

Flower outcome:
- frozen third-cohort iNaturalist measurement artifact, workflow run 35177668182.

Exposure control:
- frozen response-blind highlight technical seal, workflow run 35687421896.

Climate:
- WorldClim 2.1 10 arc-minute BIO5, using the same definition as the earlier white-environment test.

Bee interactions:
- Noori et al. 2026 curated GloBI v3.1, Zenodo record 18303036;
- exact file: `GloBI_Curated.csv`;
- expected SHA256 is recorded at execution;
- relevant columns are `bee_species`, `bee_genus`, `plant_species`, `x`, `y`, `coordinated`, and `plant_level`.

The curated source contains 981,982 bee-plant records and is versioned. A prior outcome-blind exact-name coverage audit matched 346/499 third-cohort species.

## Outcome-blind spatial preflight

Before reading morph or white/non-white outcome columns, use only flower `species`, `latitude`, and `longitude`.

Filter bee records to:
- exact plant-species matches;
- `plant_level == species_level`;
- coordinated/georeferenced records;
- finite coordinates;
- identified bee species.

Candidate local radii are fixed at 50, 100, 250 and 500 km.

At each radius and plant species, calculate the fraction of flower locations having >=3 same-plant bee-interaction records inside the radius.

A species is coverage-qualified at a radius when:
- it has >=20 flower coordinates; and
- >=50% of its flower coordinates have >=3 local bee-interaction records.

Choose the **smallest** candidate radius with >=100 coverage-qualified species. If no radius reaches 100 species, the joint test is `NOT_ESTIMABLE_STABLE_BEE_COVERAGE` and no biological model is opened.

This radius-selection rule uses no morph labels.

## Local bee predictors

After the radius is frozen by the preflight, calculate for each flower record:

1. `local_bee_richness`: number of distinct bee species among same-plant interaction records within the frozen radius.
2. `local_bombus_fraction`: fraction of distinct local bee species whose genus is `Bombus`.

Records with fewer than three local bee-interaction records have missing local predictors.

## Biological panel

After local predictors and radius are frozen, open morph outcomes.

Keep:
- globally classifiable rows in the four frozen biological states;
- rows with valid BIO5;
- rows with valid near-clip fraction;
- remove the exact response-blind high-clip set;
- rows with valid local bee predictors.

A species contributes to conditional inference only if it retains both white and non-white outcomes. The stronger species-level support audit requires >=5 white and >=5 non-white rows.

Minimum 60 informative species are required for the joint conditional models.

## Models

Within species, standardize BIO5, log1p(local bee richness), local Bombus fraction and near-clip fraction.

### Model B1 — bee richness versus heat

`white ~ z_BIO5 + z_log_local_bee_richness + z_near_clip | species`

Predictions:
- BIO5 coefficient > 0 if the broad heat clue survives independent local bee adjustment.
- bee-richness coefficient < 0 under the hypothesis that stronger/diverse bee interaction context disfavors white states.

### Model B2 — Bombus composition versus heat

`white ~ z_BIO5 + z_local_bombus_fraction + z_log_local_bee_richness + z_near_clip | species`

Predictions:
- BIO5 coefficient > 0;
- Bombus-fraction coefficient < 0.

The two pollinator-predictor p-values are Holm-adjusted. BIO5 is a previously positive discriminator and is reported as persistence/attenuation, not counted as a new discovery.

## Interpretation

- BIO5 stays positive while bee predictors are null: favors broad temperature sorting over the measured bee-community alternatives.
- bee predictor is supported and BIO5 strongly attenuates: consistent with bee-community turnover contributing to the broad heat association.
- both survive: compatible with joint abiotic and biotic sorting.
- neither survives: the earlier BIO5 result is likely scale/context dependent and should not be promoted mechanistically.

## Hard boundaries

This test does not establish:
- pigmented -> white evolutionary direction;
- a genetic pigment-loss pathway;
- causal selection;
- absence of unmeasured pollinator guilds;
- that interaction-record density equals pollination effectiveness.

The frozen New Phytologist manuscript remains unchanged.
