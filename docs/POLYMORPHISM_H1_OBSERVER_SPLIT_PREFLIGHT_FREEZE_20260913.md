# H1 observer-disjoint split preflight — frozen receipt

Date frozen: 2026-09-13 JST

Protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260913.md`

Workflow run: `34706579235`
Artifact: `10302000467` (`polymorphism-h1-observer-split-preflight-20260913`)
Artifact digest: `sha256:21e03e6413e920e0f4943a793ab48e1a375c909ded63f85293fa67825fe87c6d`

## Outcome firewall

The preflight read only:

- `species`
- `observer_id`
- `global_classifiable`

Receipt flags:

- `morph_column_opened = false`
- `palette_columns_opened = false`
- `D_computed = false`
- `D_unbiased_computed = false`
- `H1_association_or_reproducibility_computed = false`

Thus the observer split and per-half gate were fixed before any H1 outcome was opened.

## Opportunity result

Full-cohort gate: `n_classifiable >= 40`.

Discovery:

- full eligible species: 369
- species with at least two classifiable observers: 369
- min-half >=10: 369
- min-half >=15: 369
- min-half >=20: 369
- min-half count range: 20–43; median 28

Reserve:

- full eligible species: 363
- species with at least two classifiable observers: 363
- min-half >=10: 363
- min-half >=15: 363
- min-half >=20: 363
- min-half count range: 20–46; median 28

## Frozen primary gate

The preregistered automatic rule selects the largest threshold in `{20, 15, 10}` retaining at least 100 species in both cohorts.

Therefore the H1 primary gate is now fixed at:

**minimum 20 classifiable photographs in each observer-disjoint half.**

All 369 discovery and all 363 reserve species satisfy it.

## Exact preflight file hashes

- `observer_assignment_preoutcome.csv`: `sha256:f5e4333e596a69cdb032d840ec95eb02eadf1ec322658e2738e1b043fcde0cc1`
- `species_split_opportunity_preoutcome.csv`: `sha256:bb30b6f73c6ecac06022c7ddd275cd01ea71d4bc40f54dd8889814819d1f92eb`
- `result.json`: `sha256:e1a115887ec1f8453b8b971f242d78eabe04c30077ffa864e5eec48de0c2ab2f`

The H1 outcome script must deterministically reconstruct the observer assignment from the same frozen preflight code and verify these hashes before opening `morph`.
