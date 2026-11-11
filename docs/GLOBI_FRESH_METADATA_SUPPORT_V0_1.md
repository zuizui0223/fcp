# GloBI fresh metadata support — v0.1

## Decision

**PASS_GLOBI_FRESH_METADATA_SUPPORT_CANDIDATE**

Both source-semantics survivors pass the frozen metadata-only support gate.

| System | Focal taxa | Local units | Georeferenced units | Focal taxa in ≥3 units | Units with method metadata |
|---|---:|---:|---:|---:|---:|
| CarniDIET | 103 | 1,774 | 1,741 | 82 | 1,774 |
| CropPol | 48 | 3,448 | 3,076 | 44 | 2,534 |

The frozen minimums were 20 focal taxa, five georeferenced units, ten focal taxa repeated across at least three local units, and five units with method metadata.

No partner identity, partner abundance/weight, rewiring statistic, phylogenetic interaction-memory statistic or time–space coupling statistic was used.

## Meaning

The interaction arm is no longer blocked by raw sampling geometry.

The next gate is biological design eligibility:

1. define the exact local-network unit without looking at partner outcomes;
2. pin a branch-length phylogeny for the focal guild;
3. qualify context heterogeneity from independent environment/opportunity data;
4. qualify specialization from an external source or disjoint baseline;
5. qualify mobility/dispersal externally, or mark it NOT_EVALUABLE;
6. pass known-truth self-detectability.

Only afterward may partner columns be admitted.
