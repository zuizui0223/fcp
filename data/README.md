# Data layers for the frozen 34-species paper

The active empirical paper uses a layered evidence pipeline. Files in this directory should be interpreted by **stage and unit**, not as interchangeable candidate lists.

## Canonical evidence counts

| file / stage | unit | count | role |
|---|---|---:|---|
| `global_flower_colour_works.csv` + `global_flower_colour_qc.json` | retained works | 1,075 | original broad literature-discovery pool |
| `global_flower_colour_species_ranked.csv` | candidate species | 664 | high-recall species candidate layer |
| `global_flower_colour_review_queue.csv` | species | 72 | initial direct-evidence review queue |
| `resolved_inputs/global_flower_colour_review_queue_resolved.csv` | species | 111 | resolved queue after targeted follow-up/evidence aggregation |
| final manuscript dataset | species | 34 | frozen binary analysis: 20 within / 14 among, 25 families |

The final 34-species five-metric model dataset is generated/frozen by the manuscript analysis workflow. Supporting classification manifests and source audits are under `docs/supporting/`.

## Spatial states

- `within_population`: explicit evidence that multiple discrete natural flower-colour variants coexist within at least one population.
- `among_population`: geographic/among-population differentiation without retained evidence of local coexistence.
- `mixed`: both signals.
- `unclear`: evidence does not resolve the spatial state.

Mixed and unclear states are not forced into the binary manuscript response.

## Active versus historical data

The active manuscript evidence chain is the one documented above and in `docs/PIPELINE_34SPECIES.md`.

Files with names such as `candidate_screening_*`, `background_sampling_strata*`, `mass_screen_*`, and old seed/control templates record earlier development or exploratory screening. They are not final analysis inputs and are candidates for archive/removal once their scientifically necessary provenance has been captured in the canonical pipeline documentation.

A later, broader systematic-map search is documented under `literature/`. Its unreviewed expanded species sets are not primary manuscript data.
