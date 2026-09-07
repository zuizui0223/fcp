# Why real-climate synthetic qualification failed

Date: 2026-09-07 (Japan). **Post-failure synthetic diagnosis only; Stage 6 remains FAIL and observed flower colours remain closed.**

## Answer first

The diagnostic separates fixed-axis approximation from actual environmental leverage. The simulator true shared axis is privileged information and is never used for inference.

| Block / amplitude | parent detect | vector coverage median | projection corr median | train threshold exposure | eval threshold exposure | train colour-variable | eval colour-variable | true-axis augmentation Δscore median | true-axis posterior mass median |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| atmospheric_energy_dryness / 1 | 6.8% | 0.927 | 0.969 | 60.6% | 58.8% | 100.0% | 100.0% | 0.00024 | 0.188 |
| atmospheric_energy_dryness / 2 | 80.0% | 0.922 | 0.970 | 53.8% | 56.2% | 81.2% | 82.5% | 0.00183 | 0.495 |
| thermal_regime / 1 | 0.0% | 0.933 | 0.993 | 32.5% | 32.5% | 100.0% | 98.8% | -0.00002 | 0.053 |
| thermal_regime / 2 | 22.4% | 0.927 | 0.994 | 33.8% | 32.5% | 76.2% | 76.2% | -0.00012 | 0.067 |
| water_balance / 1 | 3.2% | 0.931 | 0.974 | 58.8% | 57.5% | 100.0% | 100.0% | 0.00013 | 0.154 |
| water_balance / 2 | 71.6% | 0.921 | 0.969 | 58.8% | 57.5% | 83.8% | 82.5% | 0.00148 | 0.496 |

## Guardrails

- Parent Stage 6 remains failed; no diagnostic value can reclassify it.
- Every parent synthetic label and fold score is reconstructed from the frozen seed/schedule before a diagnostic row is accepted.
- The true-axis augmentation has no inferential threshold and is not a proposed empirical axis search.
- No observed flower-colour value is read and empirical climate inference remains closed.

Workflow `34083329114` verified all 30 diagnostic shards and 3,000 parent-linked worlds.
