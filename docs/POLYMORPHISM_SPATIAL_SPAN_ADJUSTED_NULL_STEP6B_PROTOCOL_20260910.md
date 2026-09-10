# Polymorphism spatial organization — Step 6b geometry-preserving null propagation

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-spatial-span-adjusted-step6`

## Purpose

Step 6 showed that the positive D–spatial association remains directionally positive after rank adjustment for sampled span and `n_classifiable`, but its new species-label permutation P value was weaker in the independent reserve.

That species-label permutation changes the inferential null. The original Step 5 / Step 5b evidence used **within-species colour randomizations with every species' observed coordinates and sampling geometry fixed**. Those null arrays are the correct retained reference for asking whether geographic opportunity alone can generate the observed spatial association.

Step 6b therefore computes the same span-adjusted partial-rank statistic on every one of the already frozen 999 per-species spatial null realizations. No spatial randomization is regenerated.

## Immutable upstream artifacts

Discovery spatial run: `34088925008`.

Required files:
- 20 `global-rgfca-within-species-omnibus-shard-*/species_rho_permutations.csv` files;
- full-result species/global-null files for integrity checks.

Reserve spatial run: `34178957447`.

Required files:
- 20 `reserve-inference-shard-*/null.npz` arrays and receipts;
- full reserve result files for integrity checks.

## Statistic

For observed data and each frozen spatial-null realization, compute

`partial Spearman(D, species spatial rho | sampled log-span, n_classifiable)`.

Rank both focal variables and both controls, residualize each focal rank on intercept + ranked controls, then correlate residuals.

Use the exact same sampled-span definitions frozen in Step 6:

- discovery: Step-4 `log1p_span_primary`;
- reserve: maximum pairwise great-circle distance across all 100 fixed reserve photo coordinates, log1p-transformed before any morph/classifiability selection.

## Null P value

For each frame, the 999 retained spatial randomizations produce 999 adjusted partial correlations. The directional P value is

`(1 + number(null partial >= observed partial)) / 1000`.

Run for:

1. discovery primary spatial rho;
2. reserve primary spatial rho;
3. reserve matched flower-minus-background differential.

Repeat all three with `D_unbiased = D * n/(n-1)` as a predeclared sensitivity.

## Decision

This Step 6b result is the preferred span/opportunity robustness test because the null holds each species' sampling geometry fixed.

If discovery and reserve primary adjusted partial correlations remain positive and exceed their geometry-preserving nulls at P<0.05, the paper may state that the polymorphism–spatial association survives explicit sampled-span adjustment under the original within-species spatial randomization.

Reserve matched-background support is an additional flower-specific diagnostic; it is not required to define the primary spatial result.

If the adjusted geometry-preserving test fails in reserve, retain the original Step-5b replication but do not claim that the continuous D gradient is robust to explicit span adjustment.

## Claim boundary

No result here licenses adaptation, causal range-size effects, local selection, population-genetic differentiation, or a shared global boundary. The test only evaluates whether fixed sampled geographic opportunity is sufficient to explain the observed D–spatial association.
