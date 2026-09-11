# Polymorphism acquisition provenance audit

Status: **PASS**

## Proven lineage

- Historical source commit: `29584f3ad7ae0cd99a1d8f43459252af38f615da`.
- Frozen candidate pool: `data/frozen/global_monte_carlo_candidate_photos_v1.csv`; SHA-256 `f1319461d8883f3094cec8ac0e5fc247ff902b464f34146af275575b94edc9d2`.
- Candidate denominator: 1000 species × 100 photos = 100,000 rows.
- The reserve is the species-disjoint complement of the original fixed hash-ranked 500-taxon discovery budget within that same frozen candidate pool.

## Proven acquisition conditions

- iNaturalist v1 observation source; `quality_grade=research`; photographs and georeferences required; taxon `rank=species`.
- Flowering annotation gate: term 12, value 13.
- Positional accuracy ≤ 5000 m; `obscuration=none`.
- Allowed photo licences: cc0, cc-by, cc-by-sa, cc-by-nc, cc-by-nc-sa.
- Observer cap: 2 photos per species per observer; final selection: deterministic geographic maximin to exactly the inherited raw-photo target.
- Page/species selection was outcome-blind to flower colour and candidate pixels remained unopened at acquisition.

## Explicitly not established

- No native-range filter was applied by the frozen acquisition query.
- No explicit `captive=false` query parameter was applied.
- No explicit `wild=true` query parameter was applied.

Therefore Research Grade must not be rewritten as native-range-only or guaranteed wild-population sampling. The paper should interpret the spatial estimand as geographic organization in the observed community-photograph sample.
