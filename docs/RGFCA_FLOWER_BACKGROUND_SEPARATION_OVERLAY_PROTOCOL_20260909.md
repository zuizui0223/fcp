# RGFCA flower-background separation overlay — exploratory protocol

Date frozen: 2026-09-09 JST

This stage follows the recurrent-mode ecological overlay and the post-hoc palette diagnostic. It asks about the **magnitude** of photo-derived flower/background colour separation, not the direction of M1–M3 and not a prespecified island-whitening hypothesis.

## Fixed outcomes

For every admitted species × local state in each of the exact 200 RGFCA recurrence realizations, use the same observer/photo bootstrap to obtain the 12-part flower palette distribution `F` and matched-background distribution `B`. Compute two predeclared bounded separation measures from the state distributions:

1. total-variation separation: `TV = 0.5 * sum_i |F_i - B_i|`;
2. Hellinger separation: `H = sqrt(0.5 * sum_i (sqrt(F_i) - sqrt(B_i))^2)`.

Both are zero when the retained 12-colour distributions are identical and increase as they separate. They are photo-derived colour-separation measures, not pollinator receptor-space JNDs.

## Fixed predictors and estimator

Reuse exactly the ecological-overlay predictors and estimator frozen at `106ffd9986a5a2316a0909df14fb26f0dee28296`:

- absolute latitude per 10 degrees;
- island share from the exact frozen GSHHG 2.3.7 island classifications;
- predictors centred within species;
- state weight `1/k_s` so species are equal-weighted;
- through-origin slope `beta = sum(w*x*y) / sum(w*x^2)`.

The separation outcome is likewise centred within species before fitting.

## Null

Within each admitted species and realization, independently permute each predictor across its retained states. Use null seed `20260913`. Preserve the exact signal bootstrap, state support, and species support.

For TV and H × each predictor, report across 200 realizations: observed median and 2.5–97.5% quantiles, fraction positive, null median and quantiles, observed-minus-null median, and the same two-sided p-like calibration `(1 + count(|null| >= |observed median|)) / 201`.

## Claim boundary

These outcomes measure separation in the retained photo palette representation. They may be relevant to visual conspicuousness because flower and background are matched within photographs, but they do not establish pollinator-perceived contrast, attraction, adaptation, or causal island effects. Any biological interpretation requires independent receptor-space or behavioural validation.
