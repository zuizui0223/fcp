# RGFCA surviving flower-side candidate observer multi-partition stress test — result

Date: 2026-09-09 JST
Protocol freeze: `62e951c041eea7dd67efbf20969a9d076e4c6e53`

## Result

The sole local flower-side candidate `(-23,11,M2+)` did **not** survive the prospectively fixed family of 20 new observer-disjoint partitions.

- salts tested: **20/20**
- salts in which both halves passed: **2/20**
- salt-pass fraction: **0.10**
- fixed robustness requirement: **>=16/20**
- `observer_multipart_robust = false`
- median weaker-half positive fraction: **0.5614**
- weaker-half positive-fraction range across salts: **0.345–0.870**

Only salts 3 and 10 had positive fraction >=0.70 in both halves. Many other salts produced a strong positive direction in one observer half and a weak or opposite result in the other.

## Support

The metadata-only census was completed before these colour outcomes were opened. All 20 salts were retained; each had at least 159/200 realizations in which both halves contained >=5 candidate species under the fixed >=3-photo, >=2-observer half-state rule. No salt was removed based on the colour outcome.

## Reproduction

A complete second execution was byte-identical for all four retained outputs:

- `observer_multipart_global_summary.csv`: `6abdcb78d7741337fbde6239a1d6fd34484a7709f24205625b34ec4c7e7296f8`
- `observer_multipart_by_salt.csv`: `7f38c64a1571ea191bb4f647abc9788535f0b17fa8d555c6c1d9e6ff9f66f937`
- `observer_multipart_realizations.csv`: `64b34fb090d1e8f4ccf276f46d67bc9733d85be0761d65c1702fe2f9df06eed0`
- `observer_multipart_species_scores.csv`: `55bfc5ca0681222ad088160504de3a41f179f631ae7560c9456ada796926df9d`

## Interpretation

The candidate had previously survived species-disjoint, genus-disjoint, leave-one-genus, and coarse season/year matching diagnostics, but it is not stable to observer composition. Therefore the current reserve does **not** support treating this local M2+ pattern as a robust biological geographic transition.

This does not show that the underlying flower signal is absent. It shows that the present observational photo reserve cannot separate the candidate from observer/camera/geographic sampling composition strongly enough for a biological claim. The appropriate next validation target is the global mode basis itself: test whether M1–M3 remain stable when observers are split completely before any local geography is interpreted.
