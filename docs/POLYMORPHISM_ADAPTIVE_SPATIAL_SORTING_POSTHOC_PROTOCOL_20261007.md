# FCP post hoc adaptive spatial-sorting protocol — 2026-10-07

## Status

**Post-outcome ecological diagnostic.** This analysis is defined after the frozen H1/H2 and D–spatial results are known. It cannot change any frozen decision and cannot be described as prospective confirmation.

## Biological question

Does the replicated geographic organization of flower-colour variation contain evidence consistent with environmental sorting, rather than only geographic autocorrelation?

Two complementary predictions are tested in the frozen 500-species discovery and 500-species validation resources.

### Prediction A — environmental heterogeneity and species-wide colour diversity

If heterogeneous environments contribute to maintaining or generating species-wide flower-colour diversity, species sampled across more heterogeneous climates should have larger four-state D.

For each D-evaluable species:

- D is the frozen four-state species-wide sampled Gini–Simpson diversity;
- environmental heterogeneity is the median pairwise Euclidean distance in a three-variable climate space;
- climate space contains WorldClim 2.1 BIO5, BIO14 and mean monthly solar radiation;
- environmental variables are mean/SD scaled once using discovery classifiable rows and the same scaling is transported unchanged to validation;
- sampled geographic span, classifiable-row count and observer count are rank-residualized from both D and environmental heterogeneity.

The statistic is the partial Spearman correlation between D and multivariate environmental heterogeneity. A 999-permutation upper-tail test permutes the residualized environmental-heterogeneity values among species.

**Primary post hoc support:** positive partial correlation with p < 0.05 in both discovery and validation.

Single-variable BIO5, BIO14 and solar heterogeneity are sensitivities, not rescue routes.

### Prediction B — environmental distance and flower-colour distance within species

If colour variation is environmentally sorted, observations from more different environments should be more different in flower colour even after geographic distance is removed.

For each eligible species:

- flower-colour distance is pairwise Jensen–Shannon dissimilarity of the frozen four-state colour vectors;
- environmental distance is pairwise Euclidean distance in the same three-variable climate space;
- geographic distance is great-circle distance;
- the species statistic is partial Spearman(environmental distance, colour distance | geographic distance).

Inference uses matched vertex permutations of complete colour vectors among fixed coordinates. This preserves for each species:

- coordinates and geographic geometry;
- environmental values;
- the complete multiset of observed colour vectors;
- D and all univariate colour composition;
- sampling depth.

199 deterministic vertex permutations are generated per species. The cohort statistic is the equal-species mean partial rho, compared with the matched equal-species null distribution.

**Primary post hoc support:** positive observed equal-species mean with upper-tail p < 0.05 in both discovery and validation.

BIO5-, BIO14- and solar-specific pairwise distances are sensitivities, not rescue routes.

## Joint interpretation

- If A and B both replicate: environmental sorting is supported as a cross-species ecological pattern consistent with adaptive differentiation or phenotypic plasticity.
- If only B replicates: within-species environmental sorting is supported, but environmental breadth is not shown to explain between-species D differences.
- If only A replicates: environmental breadth covaries with D, but direct within-species environmental sorting is not replicated.
- If neither replicates: the existing D–spatial result remains a geographic structural pattern without added environmental-sorting support.

## Hard nonclaims

Even if both predictions are supported, this analysis does **not** establish:

- fitness differences among colour states;
- local adaptation;
- genetic differentiation rather than phenotypic plasticity;
- causal effects of BIO5, BIO14 or solar radiation;
- a universal environmental mechanism;
- independence from all unmeasured environmental, demographic or sampling processes.

A convincing local-adaptation claim would still require fitness or reciprocal/common-garden evidence, or an equivalent design connecting phenotype–environment matching to differential reproductive success.
