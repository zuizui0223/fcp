# Polymorphism directionality Step 1 — frozen protocol

Date: 2026-09-09 JST

## Question

Does discrete four-state flower-colour diversity align preferentially with the already frozen recurrent M1–M2 subspace after removing the trivial effect of total continuous colour spread?

This is the first test after the FCP mainline pivot from shared global boundaries to species-level polymorphism and the geometry of within-species variation. No new image acquisition is part of this test.

## Frozen inputs

- Discovery cohort: `data/derived/global_monte_carlo_measured_photos_v1.csv`.
- Independent axis-learning cohort: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.
- M1–M2 are reconstructed exactly from the zero-shift, observer-equal reference procedure in `scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py`; the discovery cohort is not used to rotate the axes.
- Continuous representation: the 12-dimensional flower-minus-matched-background palette proportions already used by RGFCA.
- Discrete representation: `global_morph` among `white`, `yellow_orange`, `red_pink`, `blue_purple` from the frozen global classifier.

## Fingerprint gate for D

For species i, with discrete morph proportions p_ik,

`D_i = 1 - sum_k p_ik^2`.

Before inferential results are accepted, the reconstructed discovery summary must be compatible with the already reported fingerprint: approximately 369 species, maximum D approximately 0.71, second-morph frequency >=0.10 in approximately 46.6% of species, and >=0.20 in approximately 26.6%. Exact equality is not required because upstream row eligibility may make the raw measurement file slightly broader than the previously summarized table; any material mismatch must be reported and inference treated as diagnostic rather than confirmatory.

## Species-level continuous geometry

Using discovery photos that are globally classifiable and have valid flower and background 12-anchor measurements, form the within-species covariance matrix Sigma_i of the 12D flower-minus-background vectors.

Define:

- `T_i = trace(Sigma_i)`: total continuous within-species variance.
- `v1_i = m1' Sigma_i m1`.
- `v2_i = m2' Sigma_i m2`.
- `R1_i = v1_i / T_i`.
- `R2_i = v2_i / T_i`.
- `R12_i = (v1_i + v2_i) / T_i`.

`R12` is the primary total-variance-controlled estimator. Species require at least 3 usable discovery photos and positive finite total variance.

## Primary test

Primary statistic: Spearman rho between `D` and `R12`.

Uncertainty/stress tests:

1. label-permutation p-value for rho (`20,000` permutations; fixed seed 20260909);
2. species bootstrap 95% CI (`5,000` bootstraps; fixed seed 20260909);
3. partial Spearman robustness of `D` versus `R12` controlling `log(T)` by rank residualization;
4. individual `R1` and `R2` results as secondary diagnostics.

## Direction-specificity null

To distinguish a privileged M1–M2 direction from generic multivariate spread, compare the observed `D`–`R12` rho with `10,000` random two-dimensional orthonormal subspaces drawn inside the 11-dimensional zero-sum palette subspace. For each random subspace, recompute the same variance share and its Spearman rho with D.

One-sided rotation-null p-value:

`p_rotation = (1 + number(null_rho >= observed_rho)) / (B + 1)`.

## Frozen decision rule for Claim 2

Claim 2 ('polymorphism follows a reproducible small set of colour-space directions') survives Step 1 only if BOTH are true:

1. primary `rho(D, R12) > 0` with label-permutation `p < 0.05`; and
2. M1–M2 exceeds generic random 2D orientations with `p_rotation < 0.05`.

Interpretation if the rule fails:

- primary positive but rotation-null non-significant: discrete D tracks continuous spread/anisotropy, but M1–M2 are not direction-specific; current Claim 2 is rejected or substantially weakened;
- primary non-significant: current Claim 2 is rejected;
- both pass: Claim 2 survives as a direction-specific geometric result, with M1/M2 individual effects remaining secondary.

No post-hoc axis rotation, threshold tuning, species filtering, or new acquisition is allowed to rescue this decision.