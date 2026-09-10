# Polymorphism Step 6b — opportunity decomposition protocol

Date: 2026-09-10 JST

## Status

Step 6b is a post-outcome diagnostic. It does not change the frozen Step 6 verdict `POSITIVE_BUT_NOT_TWO_TRANCHE_HOLM_SUPPORTED` for the joint sampled-span + classifiable-n adjustment.

It is frozen before opening any span-only or n-only decomposition result.

## Why decompose the two controls

The nuisance variables in Step 6 have different causal/measurement status.

- **sampled geographic span** is a pre-colour geometry/opportunity variable and directly represents the nuisance hypothesis that broadly sampled species make both polymorphism and spatial structure easier to detect;
- **n_classifiable** is the number of rows that survive the colour-classification gate. The raw photo opportunity is fixed at 100 photos per species in both discovery and reserve frames. Therefore `n_classifiable` is a post-measurement yield variable, not raw sampling effort. It can be related to colour complexity or ambiguity and is not guaranteed to be a clean pre-outcome confounder.

Consequently, the threat-specific primary diagnostic in Step 6b is **span-only adjustment**. Conditioning on `n_classifiable` is retained as a deliberately conservative measurement-yield sensitivity, not promoted to the primary opportunity correction.

## Frozen inputs

Reuse only the completed Step 6 joined tables and their frozen source tables:

- `results/polymorphism_span_adjusted_robustness_step6_20260910/discovery_joined.csv`
- `results/polymorphism_span_adjusted_robustness_step6_20260910/reserve_joined.csv`
- `results/polymorphism_species_attributes_step4_preflight_20260910/covariate_panel_preoutcome.csv`
- `data/frozen/rgfca_reserve_geometry_audit_v1.csv`

No new data acquisition is allowed.

## A. Span-only continuous D–spatial diagnostic

For discovery and reserve separately:

1. compute the partial rank correlation between D and species spatial rho controlling only log sampled span;
2. test using the same 20,000 Freedman–Lane-type residual permutations as Step 6;
3. Holm-adjust across the two tranches.

Classify the measured sampled-span nuisance as **insufficient to explain the continuous D–spatial relationship** only if both span-only partial correlations remain positive and both Holm p<0.05.

## B. Finite-sample-bias sensitivity without conditioning on n

Use unbiased Simpson diversity

`D_unbiased = D * n_classifiable / (n_classifiable - 1)`

and repeat the span-only partial-rank tests. This corrects the finite-sample bias of Simpson diversity without treating classifiable yield as a covariate.

These two p-values are Holm-adjusted across discovery and reserve as a separate sensitivity family.

## C. Fixed >=10% threshold span-only diagnostic

For discovery and reserve separately, regress raw spatial rho on the fixed second-morph >=10% indicator plus rank(span), using Freedman–Lane residual permutations under the reduced span model. Holm-adjust across the two tranches.

## D. n-only decomposition

For raw D and D_unbiased, report partial-rank D–spatial associations controlling only `n_classifiable`. These are diagnostic and are not a new claim gate.

Also report Spearman correlations:

- D with n_classifiable;
- D_unbiased with n_classifiable.

This identifies whether the loss in the joint Step 6 adjustment is attributable mainly to post-measurement yield rather than sampled geographic span.

## E. Genus clustering decomposition

For discovery and reserve separately, apply the frozen equal-genus-weighted clustering statistic to:

1. raw D;
2. D_unbiased;
3. residual D after span-only adjustment;
4. residual D_unbiased after span-only adjustment;
5. residual D after n-only adjustment;
6. residual D after span+n adjustment (already represented in Step 6, recomputed only as a fingerprint).

No new taxonomic levels or thresholds are introduced.

## Interpretation ceiling

If span-only results survive but n-conditioned results weaken, the admissible conclusion is **not** that measurement bias is absent. It is that the specific sampled-geographic-span opportunity explanation is not sufficient, while sensitivity to post-classification yield remains and must be reported.

No causal, adaptive, genetic, or universal-boundary interpretation is allowed.
