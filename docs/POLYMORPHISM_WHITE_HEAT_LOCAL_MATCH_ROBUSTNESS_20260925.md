# White–heat local matching robustness — 2026-09-25

Status: **post-result robustness analysis**.

This analysis was designed after the prospective white-environment test showed support for BIO5. It is therefore not a new confirmatory test and cannot upgrade the original causal claim. Its purpose is to ask whether the BIO5 signal is merely a broad geographic contrast within species.

## Question

After matching white and non-white records of the same species at short geographic distances, do white records still occur at warmer WorldClim BIO5 locations?

## Inputs

Use exactly the same biological and technical sources as the prospective white-environment test:
- third-cohort biological artifact from workflow run 35177668182;
- frozen response-blind high-clip set from workflow run 35687421896;
- WorldClim 2.1 10 arc-minute BIO5.

Eligible rows must:
- pass the frozen global classifiability rule;
- belong to white, yellow/orange, red/pink, or blue/purple;
- have an available near-clip metric;
- not belong to the frozen high-clip set.

## Matching

Within each species, white and non-white records are paired one-to-one without replacement by minimum great-circle distance using optimal bipartite matching.

Four maximum pair-distance thresholds are evaluated:
- 25 km;
- 50 km;
- 100 km;
- 250 km.

A pair is retained only when its distance does not exceed the threshold.

For each threshold and species, compute the mean pair difference:

white BIO5 - non-white BIO5.

Species remain the primary replicate.

## Robustness readout

For each distance threshold report:
- number of matched pairs;
- number of species with at least one pair;
- median and mean species-level white-minus-nonwhite BIO5 difference;
- fraction of species with positive difference;
- two-sided Wilcoxon signed-rank p;
- one-sided sign-test p for a positive direction.

A paired conditional logistic model is also fit across all matched rows with:
- BIO5;
- near-clip fraction;
- pair identity as the stratum.

This row-level model is corroborative only.

## Interpretation

Persistence at short thresholds supports local thermal sorting rather than a purely broad-range geographic contrast.

Loss of the signal at short thresholds indicates that the prospective BIO5 result primarily reflects larger-scale spatial sorting and does not demonstrate local heat selection.

No threshold is used to rescue or redefine the prospective mechanism gate.
