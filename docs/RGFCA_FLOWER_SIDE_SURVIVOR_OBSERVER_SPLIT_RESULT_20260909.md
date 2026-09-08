# RGFCA surviving flower-side candidate observer-disjoint nuisance audit — result

Date: 2026-09-09 JST

Protocol was frozen before outcome computation at `269845cb338b004fdd16a4687e5efa4fa6f506b7`; bootstrap seeds were fixed pre-outcome at `897ac0c2a67f64f545c099d3c687345a3df285dd`.

## Fixed target and design

Only the surviving reporting cell `(-23,11)`, M2 positive was tested. Eligible observations were split by SHA-256 parity of exact `observer_id`, so all photographs from one observer remained in exactly one half. Each half independently retained the original local-state thresholds (`>=3` photos, `>=2` observers), `>=500 km` within-species spatial support, and the previously used `MINSP=5` half-support floor. No new colour axis was fitted: each half was projected onto the already frozen full-data aligned M2 vector for the corresponding realization.

## Result

| observer half | supported realizations | species median (range) | median M2 flower score | 2.5–97.5% | positive fraction | half gate |
|---|---:|---:|---:|---:|---:|---|
| A | 191 | 8 (5–14) | +0.009494 | -0.042578 to +0.056367 | 0.6178 | FAIL |
| B | 200 | 11 (5–17) | +0.014373 | -0.015158 to +0.048805 | 0.7800 | PASS |

The prespecified candidate-level rule required both observer-disjoint halves to pass (`>=100` supported realizations and positive fraction `>=0.70`). Therefore:

**`observer_disjoint_robust = false`.**

## Reproduction

A complete same-seed rerun was byte-identical for both retained outputs:

- `observer_split_realizations.csv`: `3a6e6f57bdc43b5807dcc38e3493ff3468c007b19a275aa5acd72a98440a1d89`
- `observer_split_summary.csv`: `5a38ed0bb0d827ab0169dc2ad8e11e41f07b07bf87f582a7385c5f2d4b2ab65a`

## Interpretation boundary

This failure does not prove that observer or camera differences generate the candidate. Half A has lower taxonomic support and a broad realization distribution, so low power and observer-population composition remain possible explanations. But the candidate does not satisfy the prospectively fixed requirement that its direction recur across completely disjoint observer sets. It therefore must **not** be promoted as a nuisance-robust biological transition in the current reserve.

The surviving-cell discovery remains a within-reserve exploratory observation that passed taxonomic and coarse season/year audits but failed the observer-disjoint gate. Independent measurement-valid ecological replication is required before biological interpretation.
