# FCP selected environmental-sorting transport protocol — 2026-10-07

## Status and chronology

This protocol is frozen **after** the post hoc 500+500 environmental-sorting screen and **before** evaluating these selected statistics in the species-disjoint 499-species prospective-confirmation resource.

The third cohort has already been opened for H2 and other analyses, so this is **not untouched prospective confirmation**. It is a fixed, species-disjoint post hoc transport of two newly selected spatial-environment hypotheses.

No result here can alter frozen H1/H2/D-spatial/BIO5 decisions.

## Why exactly two targets are transported

The 500+500 screen tested one primary multivariate endpoint and three single-variable sensitivities (BIO5, BIO14, mean solar radiation) for each of two biological questions.

The joint multivariate adaptive-sorting hypothesis was not replicated.

Among the single-variable sensitivities, only two patterns were positive in both discovery and validation and remained supported after Holm correction across the three single-variable sensitivities within each endpoint and cohort:

### Target A — solar heterogeneity predicts species-wide D

Environmental-heterogeneity-versus-D screening:

- discovery SRAD: partial rho = 0.1330436, raw p = 0.004, Holm p = 0.012;
- validation SRAD: partial rho = 0.1524160, raw p = 0.001, Holm p = 0.003.

BIO5 and BIO14 did not satisfy this replicated criterion.

### Target B — BIO5 difference predicts flower-colour difference beyond geography

Within-species environment-colour sorting screening:

- discovery BIO5: equal-species mean partial rho = 0.0063943, vertex-null raw p = 0.025, Holm p = 0.050;
- validation BIO5: equal-species mean partial rho = 0.0081627, vertex-null raw p = 0.005, Holm p = 0.015.

BIO14 and solar did not satisfy this replicated criterion.

No other environmental variable is eligible for the third-cohort transport.

## Storage-schema implementation note

The first execution attempt failed before any target statistic was computed because the third-cohort table does not store the legacy pre-grouped `colour_*` vectors. It stores the nine frozen biological `flower_fraction_*` coordinates.

The target estimand is unchanged. For the four-state BIO5-sorting transport, the third-cohort script reconstructs the frozen biological groups deterministically:

- white = white;
- yellow/orange = yellow + orange + bronze;
- red/pink = red + pink + magenta;
- blue/purple = blue + purple.

The continuous-colour sensitivity uses the nine `flower_fraction_*` coordinates directly. This amendment is a storage-schema correction made before any third-cohort transport result was opened.

## Third-cohort target A — solar heterogeneity versus D

Use the frozen prospective-confirmation measured table from commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`.

Eligibility is unchanged:

- globally classifiable rows in the four biological states;
- at least 40 such rows per species.

For each species:

- D = Gini–Simpson diversity across white, yellow/orange, red/pink, blue/purple;
- solar heterogeneity = within-species SD of mean monthly WorldClim 2.1 solar radiation;
- covariates = log1p sampled geographic span, classifiable-row count, observer count.

The statistic is partial Spearman(D, solar heterogeneity | covariates), implemented by rank residualization exactly as in the screen.

A 999-permutation upper-tail residual permutation is used.

**Target-A support:** partial rho > 0 and p < 0.05.

## Third-cohort target B — BIO5 colour sorting beyond geography

For every eligible species, use:

- pairwise geographic great-circle distance;
- pairwise BIO5 absolute difference;
- pairwise flower-colour Jensen–Shannon dissimilarity from the frozen four-state colour vectors.

The species statistic is `partial Spearman(BIO5 distance, colour distance | geographic distance)`.

Inference uses 199 deterministic matched vertex permutations of complete colour vectors among fixed coordinates.

The cohort statistic is the equal-species mean partial rho.

**Target-B support:** observed mean > 0 and matched vertex upper-tail p < 0.05.

## Pre-frozen technical sensitivities for target B

These sensitivities cannot rescue a failed four-state primary transport.

If the frozen measured table contains the required columns, repeat the same BIO5-distance test using:

1. continuous flower palette-count Jensen–Shannon distance;
2. matched flower-minus-background distance, where the flower and background palette rows are permuted together as one observation.

These are intended to ask whether any transported BIO5 sorting is present in continuous flower colour and whether it exceeds same-image background colour structure.

## Multiplicity in the transport

Targets A and B are distinct biological endpoints selected before opening the third-cohort results. Raw p-values are reported for each.

For a joint statement that **both selected environmental signatures transport**, Holm correction across the two third-cohort target p-values must also leave both <= 0.05.

## Interpretation

If target B transports but the earlier white-minus-nonwhite directional BIO5 rule remains non-replicated, the allowed synthesis is:

> Warm-season temperature can repeatedly organize where flower colour changes within species without imposing one universal direction of colour change.

If target A transports, the allowed synthesis is:

> Species sampled across more heterogeneous solar environments tend to retain greater species-wide flower-colour-state diversity, conditional on sampled span and observation-depth covariates.

If both transport, they jointly support environment-associated structuring of FCP, but still do not establish local adaptation.

## Hard nonclaims

No outcome can establish:

- colour-dependent fitness;
- reciprocal local adaptation;
- genetic differentiation rather than plasticity;
- that BIO5 or solar radiation is the causal agent;
- that the same colour response direction occurs across species;
- universality across all flowering plants.
