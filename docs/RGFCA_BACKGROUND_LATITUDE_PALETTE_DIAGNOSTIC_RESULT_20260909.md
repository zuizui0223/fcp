# RGFCA background-latitude palette diagnostic — result

Date: 2026-09-09 JST
Protocol freeze: `48ae3c19f56dab4f162e0d09ad6f57244f403cef`.

This was a post-hoc descriptive decomposition of the already outcome-opened background latitude signal. Across the 12 matched-background palette fractions, the clearest recurrent latitude pattern was **green increasing with absolute latitude**: median species-equal slope `+0.025626` per 10 degrees, 2.5–97.5% `0.008147` to `0.044145`, p-like calibration `0.0398`. Brown moved in the opposite direction (median `-0.015109`, interval `-0.027890` to `-0.000621`) but its p-like calibration was `0.1144`.

All 12 flower-only palette slopes remained weak under the same species-conditioned null; none had p-like calibration below 0.52. The matched flower-minus-background green slope was negative (`-0.022721`, interval `-0.041002` to `-0.003148`, p-like `0.0647`), consistent with the background becoming greener while the flower palette itself remains comparatively stable.

A complete independent rerun reproduced all three CSV outputs byte-for-byte. Upstream mode/state reproduction differed only at floating-point roundoff (`4.44e-16`).

## Claim boundary

This diagnostic does not establish a causal vegetation mechanism or photographic bias. It shows that the strongest latitude structure in the retained matched signal is concentrated in the background palette, especially green, rather than in a corresponding flower-only latitude gradient.
