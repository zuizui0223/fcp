# Polymorphism M3 Step 3 — frozen diagnosis protocol

Date: 2026-09-09 JST

## Question

What is the species-level M3 signal that appeared to behave differently from M1/M2: a flower-biological colour axis, a flower-internal nuisance/possible developmental signal, or photographic/background context?

This is a diagnosis step only. Step 1 already rejected Claim 2 as a direction-specific polymorphism result, so M3 cannot rescue that claim.

## Frozen axis

M3 is reconstructed by the exact zero-shift, observer-equal reference procedure in `scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py` using the reserve cohort. No discovery outcome is used to rotate M3.

The previously frozen M3 is dominated by brown/black/yellow coordinates in the flower-minus-background representation. This diagnosis asks which component carries any D association.

## Discovery test: flower-side source of the association

Use the 369-species discovery frame from Step 1. For each classifiable photo, calculate:

- `M3_full`: projection of the 12-anchor flower composition onto frozen M3;
- `nuisance_fraction`: flower palette fraction assigned to green+brown+black;
- `M3_bio9`: set green+brown+black to zero, renormalize the nine biological anchors to sum to one, and project onto the same frozen M3 without rotation.

Per species calculate mean and SD for `M3_full`, mean `M3_bio9`, and mean nuisance fraction.

Primary descriptive statistics are Spearman correlations of D with mean M3_full and SD(M3_full). The diagnosis statistics are:

1. `rho(D, mean M3_bio9)`;
2. partial Spearman `D ~ mean M3_full | mean nuisance_fraction`, by rank residualization;
3. `rho(D, mean nuisance_fraction)`.

For the first two M3 associations, use 20,000 label permutations with fixed seeds 20260909–20260911. No threshold tuning after outcome opening.

## Reserve component diagnosis

Using reserve rows with valid flower and matched-background 12-anchor counts and the reserve species with >=40 classifiable morph rows, calculate species means of:

- flower-only M3 score F;
- background-only M3 score B;
- flower-minus-background M3 score Dsig = F-B;
- flower nuisance fraction.

Report Spearman correlations of discrete morph diversity D with mean F, mean B, mean Dsig and flower nuisance fraction. Also report partial Spearman `D ~ mean F | mean B`.

These reserve component results are diagnostic rather than independent tests because M3 itself was learned from the reserve signal geometry.

## Frozen diagnosis rule

- `FLOWER_BIOLOGICAL_SUPPORTED` only if the discovery biological-only association remains positive and permutation p<0.05, the nuisance-controlled association remains positive and p<0.05, and reserve |rho(D,F)| > |rho(D,B)|.
- `BACKGROUND_CONTEXT_SUPPORTED` if the discovery full M3 association is present but either the biological-only or nuisance-controlled association fails, and reserve |rho(D,B)| >= |rho(D,F)|.
- otherwise `MIXED_OR_UNRESOLVED`.

A finding tied to flower nuisance fraction but not background may be called `flower-internal nuisance / possible developmental-stage signal`; it cannot be called flowering stage because the frozen tables contain no direct phenological-stage label.

No new images, manual stage scoring, axis rotation, species filtering, or post-result threshold selection are allowed in this step.
