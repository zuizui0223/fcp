# RGFCA within-species spatial organization omnibus

Date: 2026-09-07 (Japan). **RGFCA-only postoutcome exploratory exact randomization test; six-species and 34-species analyses are not used.**

## Answer first

Equal-species observed mean distance–colour Spearman rho = **0.02702**; exact within-species permutation upper-tail **p = 0.0010**.

Null mean = 0.00004; 95% null interval = [-0.00601, 0.00590].

Observed median species rho = 0.01640; fraction with rho > 0 = **60.4%**.

BH-detectable species at exploratory q<.05: **9/369**. This is explicitly a detectable-species count, not a prevalence estimate.

## Interpretation guardrails

- The primary null keeps each species, every coordinate, photo count and complete colour-vector multiset fixed; only colour vectors are reassigned among coordinates within species.
- The primary unit is species: 369 species receive equal weight regardless of photo or pair count.
- This tests a monotone geographic-distance versus colour-dissimilarity tendency. Non-monotone patchy structure can be missed.
- The analysis is postoutcome exploratory because the descriptive G3 direction was already known; it is not confirmatory.
- It is not a biological prevalence estimate, shared-boundary test, environmental mechanism test or causal test.
- G1 and all prior frozen decisions remain unchanged.

Workflow `34088925008` verified all 369 observed species statistics, all 999 synchronized global null statistics, and 12 direct SciPy Spearman spot checks.
