# Stable GloBI spatial preflight for the FCP white–pollinator mechanism — 2026-09-25

Status: **outcome-blind spatial preflight frozen after the live-API coverage audit and before any pollinator metric is joined to flower-colour outcomes.**

## Purpose

The live GloBI coverage audit passed all frozen identifiability gates, including 68/499 third-cohort plant species with a Sphingidae interaction and 295 species with substantial georeferenced pollinator coverage. The live API is not used for biological inference.

This preflight rebuilds the pollinator layer from one immutable, citable GloBI release and determines whether a within-species spatial hawkmoth test is identifiable.

## Frozen interaction release

Use exactly:

- Global Biotic Interactions Review Dataset Corpus, Zenodo record **18064921**, version v3, published 2025-12-26;
- file: `interactions.tsv.gz`;
- expected MD5: `7c12420410e0fea43608308c387fa89c`;
- file size reported by Zenodo: approximately 716 MB.

No live GloBI API result enters the biological test.

## Interaction extraction

Retain only these interaction types:

- `pollinates`;
- `visitsFlowersOf`;
- `pollinatedBy`;
- `flowersVisitedBy`.

Use exact canonical plant-name matching against the 499 third-cohort species.

Normalize orientation so every retained row is:

`plant species -> pollinator taxon -> interaction coordinates`.

Guild classification is the same taxonomic-path rule frozen in the live coverage protocol. Sphingidae is the focal hawkmoth class.

## Outcome firewall

This preflight may read from the third-cohort measurement table only:

- species;
- photo_id;
- latitude;
- longitude.

It must not read morph, colour fractions, D, H2 vectors, global_classifiable, exposure variables, or WorldClim values.

## Spatial descriptors

For every third-cohort photo coordinate, and separately within each plant species, compute great-circle distance to:

1. nearest georeferenced interaction of any retained pollinator;
2. nearest georeferenced bee interaction;
3. nearest georeferenced Lepidoptera interaction;
4. nearest georeferenced Sphingidae interaction.

Distances are descriptors of documented interaction geography, not pollinator abundance or true absence.

The planned focal covariate for a later biological test is:

`
relative_hawkmoth_distance
  = log1p(nearest_sphingidae_km) - log1p(nearest_any_pollinator_km)
`

This asks whether a location is disproportionately close to documented hawkmoth interactions relative to the plant's general pollinator-documentation geography.

## Frozen spatial gate

A plant species is `hawk_spatial_eligible` if the stable release contains:

- >=1 georeferenced Sphingidae interaction;
- >=5 georeferenced any-pollinator interaction records;
- >=2 distinct 1-degree any-pollinator cells;
- finite nearest-any and nearest-Sphingidae distance for its third-cohort photo coordinates.

The future hawkmoth mechanism test is authorized only if >=30 plant species are `hawk_spatial_eligible`.

A stronger multi-site sensitivity is estimable only if >=20 species also contain >=2 distinct 1-degree Sphingidae cells.

For a broader Lepidoptera mechanism, require >=50 species with >=1 georeferenced Lepidoptera interaction and the same >=5-record / >=2-cell general pollinator coverage.

## No scale tuning

The primary metric is continuous great-circle distance and therefore has no fitted search radius.

For descriptive coverage only, report the fraction of third-cohort photo coordinates within 50, 100, 250 and 500 km of an any-pollinator, Lepidoptera, and Sphingidae interaction. These thresholds cannot be selected later as alternative primary predictors.

## Hard nonclaims

Passing this preflight does not establish:
- hawkmoth preference for white;
- pollinator-driven flower-colour evolution;
- an evolutionary pigmented-to-white direction;
- true pollinator absence where Sphingidae is undocumented;
- independence from temperature or other environmental gradients.

If the spatial gate passes, the subsequent biological protocol must be frozen before morph outcomes are joined. It must include BIO5 and exposure control so that hawkmoth geography is tested against the heat-sorting candidate rather than in isolation.
