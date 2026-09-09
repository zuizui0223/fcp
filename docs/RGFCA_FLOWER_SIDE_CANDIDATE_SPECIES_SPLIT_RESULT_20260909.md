# RGFCA flower-side candidate species-split robustness — result

Date: 2026-09-09 JST
Protocol freeze: `d416e94939b4a2e1218b2ca8f442450b3c987c0c`.

The deterministic SHA-256 species split contained 253 species in A and 247 in B. Reusing the exact 200 recurrence realizations and the fixed 13 flower-side candidate cells, only **3/13 candidates passed the cross-split robustness gate** (>=100 supported realizations and >=0.70 recurrence of the original full-data sign in both halves).

Passing candidates:

- `(rx=-24, ry=11)`, M1, positive: A sign recurrence 0.765, B 0.795;
- `(rx=-23, ry=11)`, M2, positive: A 0.980, B 0.715;
- `(rx=-19, ry=8)`, M2, positive: A 0.7919, B 0.7473.

The previously highlighted western North America M3 adjacent sign-reversal pair did **not** pass. At `(rx=-24, ry=8)` the negative M3 signal recurred 0.99 in A but only 0.46 in B; at `(rx=-24, ry=9)` the positive M3 signal recurred 0.83 in A and 0.665 in B. It should therefore be downgraded from a robust boundary candidate to a taxonomically unstable exploratory pattern.

Several other full-data candidates were similarly dominated by one species half, including the Florida M1/M2 pair and Australian M2/M3 candidates.

A complete rerun reproduced all four output CSVs byte-for-byte; upstream recurrence support matched to floating-point roundoff (`4.44e-16`).

## Claim boundary

Passing this post-discovery split is a robustness check, not independent confirmation. The three surviving cells remain exploratory and require further lineage-disjoint and independent-data validation before any biological boundary or adaptation claim.
