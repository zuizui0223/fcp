# RGFCA recurrent-signal ecological overlay — exploratory protocol

Date frozen: 2026-09-09 JST

This is a downstream interpretation stage. The signal modes, recurrence, spatial reporting cells, and flower/background decomposition were all frozen before this protocol. No ecological variable was used to choose or rotate M1–M3.

## Fixed signal outcomes

Reuse the exact 200 signal realizations (`seed=20260909`) and the already frozen reference-aligned modes M1–M3.

For each admitted species × local state retain:

- `S_F`: species-centred flower-only score projected onto the frozen mode;
- `S_D`: species-centred flower-minus-background score;
- `S_B`: species-centred matched-background score;

with the exact identity `S_D = S_F - S_B`.

`S_F` is the primary ecological-response component because the completed decomposition showed that M1 and M3 matched signals are globally background-heavy. `S_D` is a companion matched-signal response and `S_B` is a diagnostic for background-driven associations.

## Fixed ecological predictors

Only two global predictors are tested in this stage:

1. **absolute latitude**, defined as the mean absolute photo latitude of the eligible photos forming a local state and scaled per 10 degrees;
2. **island share**, defined as the fraction of eligible photos in a local state that fall within the exact frozen GSHHG 2.3.7 island polygons used by the `island` project.

The predictor values are computed from all eligible photos in each state, not from the outcome bootstrap sample.

For every admitted species and realization, each predictor is centred by that species' mean across its admitted local states. Thus slopes represent **within-species geographic association**, not species turnover.

## Species-equal slope estimator

If species `s` has `k_s` admitted states, each state receives weight `1/k_s`. For each mode, component, predictor, and realization fit the weighted through-origin slope on species-centred predictor and species-centred signal score:

`beta = sum(w*x*y) / sum(w*x^2)`.

Species with zero within-species variation in a predictor contribute no denominator for that predictor. Report the number of informative species per realization.

## Species-conditioned null

Use null seed `20260911`. Within each admitted species and realization, independently permute each predictor across that species' retained state locations before fitting the same slope. This preserves the signal scores, predictor distribution, number of states, and species support while breaking their within-species correspondence.

For each mode × component × predictor report across 200 realizations:

- observed slope median and 2.5–97.5% quantiles;
- fraction of observed realization slopes >0;
- null slope median and 2.5–97.5% quantiles;
- observed-minus-null median slope;
- two-sided empirical calibration against the null using the absolute observed median slope as the scalar: `(1 + count(|null slope| >= |observed median slope|)) / 201`.

This p-like quantity is exploratory calibration, not a confirmatory family-wise-error-controlled test.

## Metadata-only support gate already checked

Before opening any ecological slope outcomes, 20 random grid shifts were inspected using coordinates and land metadata only.

- species with island-share range >0.2: 84–102 per shift;
- species with both island-dominant (>=0.8) and mainland-dominant (<=0.2) retained states: 43–54 per shift;
- species with retained states on >=2 distinct dominant islands: only 14–20 per shift.

Therefore absolute latitude and island share are admitted for this overlay. Island area, isolation, climate PCs, and analysis regime are **not** fitted here because the current reserve has insufficient clean within-species multi-island support for a stable multivariable interpretation.

## No candidate-cell cherry-picking

All admitted M1–M3 states enter the global slopes. The 13 post-hoc flower-side candidate cells and the western-North-America M3 transition candidate are not used to select observations or predictors.

## Claim boundary

Any association found here is exploratory and within-species in the current photo-derived reserve. It does not establish adaptation, pollinator attraction, causal island effects, or a universal latitudinal rule. Independent measurement validation and an independent ecological replication tranche remain necessary for stronger claims.
