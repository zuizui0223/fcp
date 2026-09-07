# Real-climate synthetic response-sharing qualification

Date: 2026-09-07 (Japan). **Synthetic labels only; observed flower colours remain closed.**

## Answer first

Overall frozen synthetic qualification gate: **FAIL**.

The same 21,418 jointly climate-complete retained photos, species split, four longitude holdouts, 500-km buffer and 20 photos/species were used for all three climate blocks.

| Block | amp=1 full shared, threshold SD=.5 | amp=2 full shared, threshold SD=.5 | worst nuisance rejection |
|---|---:|---:|---:|
| thermal_regime | 0/250 = 0.0% | 56/250 = 22.4% | 0/250 = 0.0% (no_structure) |
| water_balance | 8/250 = 3.2% | 179/250 = 71.6% | 0/250 = 0.0% (no_structure) |
| atmospheric_energy_dryness | 17/250 = 6.8% | 200/250 = 80.0% | 8/250 = 3.2% (geographic_independent) |

## Guardrails

- The primary threshold is one common threshold per model: the maximum 97.5th percentile across all 36 block-by-nuisance calibration distributions.
- No block-specific decision threshold, best block, best fold, replacement photo, or post-outcome axis search is used.
- Species-specific threshold heterogeneity and all richer nuisance worlds remain in the benchmark.
- Real CHELSA values are used only as predictor geometry for artificial labels; no biological flower-colour value is read.
- Passing this benchmark does not itself open empirical inference; a separate observed-colour contract is still required.

Workflow `34082752756` recomputed all 72 nuisance quantiles, both global thresholds and all 108 evaluation decisions after all 60 shards completed.
