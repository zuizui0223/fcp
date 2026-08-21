# Data layers for the frozen 34-species paper

The active paper uses three data layers. They are separated by **unit and purpose** so literature records, candidate species and final statistical rows are never conflated.

## 1. Durable statistical freeze

The only production input to the current 34-species climatic analysis is:

- `frozen/frozen_34species_five_metric_dataset.csv`
- SHA-256: `bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed`
- 34 species, 25 families
- 20 `within_population`, 14 `among_population`
- all rows `baseline_unambiguous`
- >=20 occupied climate cells per species
- five climatic-niche metrics
- canonical row order: `canonical_name` ascending

`frozen/freeze_manifest.json` records the recovery history and checksum. The production workflow reads this committed file directly; it no longer depends on a short-lived GitHub Actions artifact.

## 2. Original literature-evidence provenance

| file / stage | unit | verified count | role |
|---|---|---:|---|
| `global_flower_colour_works.csv` + `global_flower_colour_qc.json` | retained works | 1,075 | original broad literature-discovery pool |
| `global_flower_colour_species_ranked.csv` | candidate species | 664 | high-recall species candidate layer across 140 families |
| `global_flower_colour_review_queue.csv` | species | 72 | initial direct-evidence review queue |
| `resolved_inputs/global_flower_colour_review_queue_resolved.csv` | species | 111 | resolved queue after targeted follow-up/evidence aggregation |
| durable freeze above | species | 34 | final binary comparative analysis |

The 72- and 111-species files are evidence-screening stages, not inferential samples. Mixed, unclear or otherwise non-binary cases were not forced into the final response.

## 3. Broader search-completeness infrastructure

A later systematic-map search is documented under `literature/` and its acquisition code under `scripts/literature/`. It used 15 query blocks and 52 shards and produced 79,242 deduplicated bibliographic records. That later corpus is retained as search-completeness/provenance infrastructure; its unreviewed expanded species sets are not the current statistical dataset.

## Spatial states

- `within_population`: explicit evidence that multiple discrete natural flower-colour variants coexist within at least one population.
- `among_population`: geographic/among-population differentiation without retained evidence of local coexistence.
- `mixed`: both signals.
- `unclear`: evidence does not resolve the spatial state.

The current labels are source-traceable and rule-derived. Completed independent blinded human review is not claimed unless completed reviewer sheets are actually supplied.
