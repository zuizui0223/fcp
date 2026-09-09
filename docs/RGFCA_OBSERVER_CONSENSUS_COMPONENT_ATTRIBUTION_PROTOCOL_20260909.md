# RGFCA observer-consensus geography — flower/background component attribution protocol

Date frozen: 2026-09-09 JST

This is a fixed follow-up to the observer-consensus geography result. Exactly three already-open cell × mode combinations enter:

- `(-16,6), M2+`;
- `(29,-9), M1+`;
- `(-23,8), M1+`.

No new cells or modes are searched.

## Fixed decomposition

For every eligible photo retain the 12-colour flower composition `F`, matched-background composition `B`, and `D = F - B`.

Reuse exactly the same 20 observer partitions, both halves, first 20 frozen grid shifts, 300 km state construction, >=3 photos and >=2 observers per half-state, >=500 km within-species state separation, and 500 km reporting lattice used in the observer-consensus geography.

Within each half-state calculate equal-observer means separately for `F`, `B`, and `D`. Within every retained species, centre each component across that species' admitted states. Require the numerical identity `D_centered = F_centered - B_centered` to tolerance `1e-12`.

Project centred states onto the frozen full-data reference mode for the fixed cell (M1 or M2). Define background contribution as `-B_score`, so the exact score identity is:

`D_score = F_score + (-B_score)`.

## Observer-consensus component rule

For each fixed cell × mode and each component `F` and `-B`, use the same rule as the geographic discovery:

- half-cell-shift supported at >=5 species;
- half evaluable at >=10/20 supported shifts;
- half positive when >=70% of supported shift scores are >0, negative when <=30%, otherwise unstable;
- salt agreement requires both observer halves directional with the same sign;
- component is observer-consensus positive when >=16 salts are evaluable and >=80% of evaluable salts agree positive.

## Classification

For each of the three fixed D-positive geographic survivors:

- `flower-supported` if `F` is observer-consensus positive and `-B` is not;
- `background-supported` if `-B` is observer-consensus positive and `F` is not;
- `both` if both are observer-consensus positive;
- `unresolved` if neither is observer-consensus positive.

The original D observer-consensus result must also reproduce under this component runner.

## Claim boundary

This attribution distinguishes which photo-derived component carries the stable geographic direction. `flower-supported` still does not establish pollinator perception, adaptation, or validated petal colour because independent flower-region measurement validity remains a separate gate.
