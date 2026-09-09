# RGFCA sharedness-specific predictive qualification

Date: 2026-09-07 (Japan). **RGFCA mainline only; six-species and 34-species analyses are not used in this stage. Observed flower colours remain closed.**

## Answer first

Overall frozen synthetic qualification gate: **FAIL**.

The null explicitly allows every species to possess its own strong spatial boundary. The statistic only asks whether a common boundary learned from 100 training species predicts 100 disjoint evaluation species better than that structured independent-boundary baseline.

Global decision threshold: `0.013962128`.

| Arm | detected / 250 | rate |
|---|---:|---:|
| partially_shared_boundaries, amp=0.5, shared=0.50, threshold_sd=0.25 | 0 | 0.0% |
| partially_shared_boundaries, amp=1.0, shared=0.10, threshold_sd=0.25 | 0 | 0.0% |
| partially_shared_boundaries, amp=1.0, shared=0.25, threshold_sd=0.25 | 0 | 0.0% |
| partially_shared_boundaries, amp=1.0, shared=0.50, threshold_sd=0.25 | 0 | 0.0% |
| partially_shared_boundaries, amp=1.0, shared=0.75, threshold_sd=0.25 | 1 | 0.4% |
| partially_shared_boundaries, amp=1.0, shared=1.00, threshold_sd=0.25 | 2 | 0.8% |
| partially_shared_boundaries, amp=2.0, shared=0.25, threshold_sd=0.25 | 0 | 0.0% |
| partially_shared_boundaries, amp=2.0, shared=0.50, threshold_sd=0.25 | 13 | 5.2% |
| partially_shared_boundaries, amp=2.0, shared=1.00, threshold_sd=0.25 | 178 | 71.2% |
| shared_common_offset, amp=1.0, shared=1.00, threshold_sd=0.00, common_offset=0.25 | 20 | 8.0% |

## Structured-null control

Worst evaluation nuisance rejection: **8/250 = 3.2%** in `directionally_clustered_independent` (amp=2.0).

## Guardrails

- `shared_fraction=0, amplitude>0` independent-boundary worlds are calibration nulls, not alternatives.
- Training and evaluation species are disjoint fixed pools (184/185); each replicate samples 100/100 species and 20 photos/species.
- Evaluation labels never update the training posterior.
- One threshold is the maximum 97.5th percentile across all nine structured-null calibration arms; there is no amplitude-specific or favourable-axis threshold.
- No six-species or 34-species result is used as evidence or input in this stage.
- No biological flower-colour value or real climate value is read.
- A PASS would qualify the synthetic statistic only; it would not reclassify G1 or automatically open observed-colour inference.

Workflow `34084872505` recomputed all 9 null quantiles, the global threshold, all 19 evaluation decision counts, and the species/photo inclusion censuses after all 20 shards completed.
