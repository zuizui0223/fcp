# Hypervolume response robustness to shifted thresholds and richer nulls

Date: 2026-09-07 (Japan). **Synthetic robustness complete; observed-colour and real-climate inference remain closed.**

## Answer first

The prespecified intercept-marginal response model was calibrated against all twelve structured nuisance arms. Its largest observed nuisance rejection was **6/250 = 2.4%** in `mixed_independent` (amplitude 2.0).

Positive recovery under species-specific threshold shifts is reported below. These are artificial effects on the actual coordinate support, not biological flower-colour evidence.

| Synthetic positive arm | threshold SD/common offset | detected / 250 | rate |
|---|---:|---:|---:|
| environmental_shared_shifted, amp=1.0, shared=0.25 | sd=0.5 | 3 | 1.2% |
| environmental_shared_shifted, amp=1.0, shared=0.50 | sd=0.5 | 25 | 10.0% |
| environmental_shared_shifted, amp=1.0, shared=1.00 | sd=0.5 | 164 | 65.6% |
| environmental_shared_shifted, amp=1.0, shared=1.00 | sd=1.0 | 53 | 21.2% |
| environmental_shared_shifted, amp=2.0, shared=1.00 | sd=0.5 | 250 | 100.0% |
| environmental_shared_common_offset, amp=1.0, shared=1.00 | common=0.5 | 144 | 57.6% |

## Guardrails

- Training and evaluation species remain disjoint and evaluation regions use the existing 500-km buffer.
- The candidate integrates a fixed five-point threshold grid; no threshold was selected after outcomes.
- Geographic, mixed, shifted-independent, support-aligned, directionally clustered, nonlinear and no-structure nuisance worlds all remain in calibration.
- No biological colour value or real climate raster was read, and no existing empirical support decision was changed.
- A favourable result here is still insufficient to open real-climate inference without a separately frozen empirical contract.

GitHub workflow `34080881485` finalized the 20 non-overlapping shards and recomputed every threshold and decision count.
