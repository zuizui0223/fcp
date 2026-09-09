# RGFCA observer-consensus geography — result

Date: 2026-09-09 JST
Protocol freeze: `8b866383c8afff24f7a14968ace306cd243bdb60`

## Result

Rebuilding the geographic map independently across 20 observer partitions and both halves removed the previous pooled-map Canadian Rockies candidate. Only three cell × mode combinations met the prespecified observer-consensus rule.

| reporting cell | approximate centre | mode | evaluable salts | positive agreements | dominant agreement fraction |
|---|---|---|---:|---:|---:|
| `(-16,6)` | 80.3°W, 26.4°N | M2+ | 20 | 20 | 1.00 |
| `(29,-9)` | 152.9°E, 35.5°S | M1+ | 20 | 20 | 1.00 |
| `(-23,8)` | 116.6°W, 35.5°N | M1+ | 20 | 16 | 0.80 |

The first cell is in the South Florida region; the second lies in the southeastern-Australia / NSW-coast reporting region; the third is in southern California / Mojave-region coordinates. These labels are geographic orientation only, not ecological interpretation.

The prior `(-23,11,M2+)` Canadian Rockies candidate did **not** survive observer-consensus reconstruction.

## Reproduction

A complete second execution was byte-identical for the four science-bearing outputs:

- `observer_consensus_geography_summary.csv`: `583f0426a0cfed742a1d81c02c270cce52afd8fa19dcc1f0c19c5e0de850bd0f`
- `observer_consensus_salt_agreement.csv`: `3d56df012069c72ddfc19595ada0590fe1eda039e5abdd58fc4dc72cb9984c94`
- `observer_consensus_half_cell_summary.csv`: `c4369ca72ef4665b0513674baba8187be22b6df15360fce7485f369cc3c5da30`
- `observer_consensus_cell_shift_scores.csv`: `c80b7d23f40beff068571b5779a892c4a1c0bfc0a0b226a8624ca1acc08a0234`

## Interpretation

The stable result is not a single universal boundary. It is a sparse set of local recurrent mode directions that survive observer-disjoint reconstruction. Before any biological interpretation, these three cells must be decomposed into flower-only and background-only contributions under the same observer-consensus framework.
