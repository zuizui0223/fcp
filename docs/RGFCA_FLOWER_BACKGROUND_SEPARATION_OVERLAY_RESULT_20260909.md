# RGFCA flower-background separation overlay — result

Date: 2026-09-09 JST
Protocol freeze: `188966c7a0710323d405a929a1c5a05b23437496`.

The predeclared total-variation (TV) and Hellinger distances between the 12-part flower and matched-background palette distributions showed **no global within-species latitude or island-share association** across 200 recurrence realizations.

| Metric | Predictor | Median slope | 2.5% | 97.5% | p-like |
|---|---|---:|---:|---:|---:|
| Hellinger | absolute latitude / 10° | +0.000801 | -0.008800 | +0.012643 | 0.9254 |
| Hellinger | island share | -0.001328 | -0.020931 | +0.018428 | 0.9353 |
| TV | absolute latitude / 10° | +0.002545 | -0.011269 | +0.017059 | 0.7761 |
| TV | island share | -0.000791 | -0.021076 | +0.020944 | 0.9602 |

Thus the recurrent geographic signal found earlier is not a simple global increase or decrease in overall flower/background palette separation. The evidence is more consistent with **regional reorganization of the composition/direction of contrast** than with a universal conspicuousness-magnitude gradient.

A complete independent rerun reproduced all three CSV outputs byte-for-byte and the upstream recurrence realization support was reproduced to floating-point roundoff (`4.44e-16`).

## Claim boundary

TV and Hellinger are photo-derived palette-separation measures. They are not pollinator receptor-space JNDs and do not establish attraction or adaptation.
