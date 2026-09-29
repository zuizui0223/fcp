# GloBI spatiotemporal dataset-index preflight — v0.1

## Decision

**HOLD_GLOBI_DATASET_INDEX_LACKS_GEOMETRY_TIME_METADATA**

The stable GloBI v0.11 archive was checked at the dataset-index level only.

Source identity:

- Zenodo record: 22691479
- version: 0.11
- published: 2026-09-10
- file: `datasets.tsv`
- SHA256: `6b5191180e9bbfdb7867ff2fdef3e2e305f78cfe6f6405209e455497a06e30c4`

## What was opened

Only the small dataset index was inspected.

It contains:

- 458 dataset namespace rows;
- one column: `namespace`.

No interaction row, taxon identity, edge weight, local coordinate, event date or biological turnover response was opened.

## Why the gate holds

The frozen question was whether the stable dataset index itself exposes enough dataset-level spatial or temporal metadata to identify candidate repeated-network systems before interaction records are touched.

It does not.

There are no geometry/time columns in the index.

Therefore the programme does **not** open the 2.5-GB interpreted interaction table merely to find a candidate after the gate failed.

## Next legitimate gate

GloBI remains the second and final interaction source family under the frozen finite-source rule.

A bounded continuation is allowed **within GloBI**:

1. use dataset namespace identities only;
2. freeze a finite original-source / dataset-cache metadata candidate set;
3. inspect only repository/source metadata needed to establish network/site/date/sampling semantics;
4. keep interaction records sealed;
5. stop if no source-complete paired time-space system qualifies.

This is not a source switch and does not authorize interaction-outcome opening.

## Boundary

This HOLD says only that the integrated GloBI dataset index is too thin for spatiotemporal eligibility screening by itself.

It is not evidence that GloBI contains no suitable ecological systems.

Current FCP manuscript science and frozen CHUN science are unchanged.
