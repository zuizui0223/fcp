# RGFCA surviving flower-side candidate observer-disjoint diagnostic

Date: 2026-09-09 JST

This is an outcome-opened post-hoc nuisance diagnostic of the sole surviving local flower-side candidate `(-23,11,M2+)`. It is **not** a prospectively frozen confirmatory test.

## Construction

Observers were assigned wholly to deterministic SHA-256 parity halves. Within each half, local states were rebuilt from that half alone using the existing 300 km grid shifts, minimum 3 photos and 2 observers per state, >=500 km within-species state separation, and the frozen realization-specific M2 loading. Each supported half required at least 5 candidate species; this lower support threshold is diagnostic and reflects the strong support loss after splitting observers.

## Result

| half | supported realizations | species median (min–max) | score median | 2.5–97.5% | positive fraction | half pass |
|---|---:|---:|---:|---:|---:|---|
| A | 191 | 8 (5–14) | +0.008790 | -0.044024 to +0.048794 | 0.6230 | no |
| B | 200 | 11 (5–17) | +0.013510 | -0.019231 to +0.049968 | 0.7500 | yes |

Therefore `observer_disjoint_robust=false` for this arbitrary deterministic split.

## Support diagnostic

A metadata-only support census showed that requiring the original >=10 species simultaneously in both halves is often infeasible after splitting observers: with 3 photos and 2 observers per half-state, both halves reached >=10 candidate species in only 3/200 realizations. This is why this result must be interpreted as a stress test rather than a clean replication.

## Reproducibility

A complete second execution was byte-identical for all three retained outputs:

- `observer_split_summary.csv`: `d55612e23da924bd20d3080a8d5f906f59667b7284ca1ba489967e2bc2c887da`
- `observer_split_realizations.csv`: `d4c14e5c6872c8e8f8be5dc54961ca1a9e49e309bd76394ee8e8481365adaa33`
- `observer_split_species_scores.csv`: `056b768908441874ff4a840131f3226106cd2bec63f873e8b3668c2f6875315e`

## Interpretation

The season/year-matched candidate survives coarse temporal matching, but the signal is **not stable across one observer-disjoint partition**. This leaves observer/camera/geographic sampling composition as an unresolved nuisance source. The candidate should therefore remain exploratory and should not be promoted to a biological transition claim from the current reserve alone.
