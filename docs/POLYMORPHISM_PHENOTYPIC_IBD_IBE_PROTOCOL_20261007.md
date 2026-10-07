# FCP phenotypic IBD versus IBE-like partition protocol — 2026-10-07

## Motivation

Geographic flower-colour differentiation can arise because populations are farther apart (IBD-like structure), because they occupy different environments (IBE-like environmental sorting), or both.

This analysis borrows the IBD/IBE logic from landscape genetics but applies it to flower-colour phenotype distances. It must therefore be described as **phenotypic IBE-like sorting**, not genetic isolation by environment.

## Data and representation

Use all three frozen FCP high-depth cohorts.

Species eligibility:
- globally classifiable;
- at least 40 retained rows;
- valid coordinates and BIO5.

Flower colour is represented continuously using the same nine biological palette coordinates: white, yellow, orange, red, pink, magenta, purple, blue, bronze.

Legacy discovery/validation use palette_count_*; the third cohort uses flower_fraction_*. Rows are normalized before pairwise Jensen–Shannon dissimilarity.

## Pairwise quantities within each species

For every unordered observation pair:

- Y = flower-colour Jensen–Shannon dissimilarity;
- G = great-circle geographic distance;
- E = absolute WorldClim BIO5 difference.

## Species statistics

### IBE-like unique association

rho_IBE = partial Spearman(Y, E | G)

### IBD-like unique association

rho_IBD = partial Spearman(Y, G | E)

### Relative balance

delta = rho_IBE - rho_IBD

Positive delta means the unique BIO5 association exceeds the unique geographic-distance association for that species.

## Matched null

For each species, permute complete colour vectors among fixed coordinates 199 times.

This preserves geography, BIO5 values, geography–BIO5 collinearity, the exact colour-vector multiset, sample size and overall colour distribution.

Both rho_IBE and rho_IBD are recalculated in every null world.

## Cohort-level inference

Species contribute equally.

For each cohort report:

- equal-species mean rho_IBE and upper-tail matched-null p;
- equal-species mean rho_IBD and upper-tail matched-null p;
- fraction of species with positive rho_IBE / rho_IBD;
- mean delta = rho_IBE - rho_IBD;
- sign-flip p for whether delta differs from zero;
- descriptive rank-distance commonality partition.

## Rank-distance commonality partition

Within each species regress rank(Y) on rank(G), rank(E), and both.

Calculate:

- R2_G;
- R2_E;
- R2_GE;
- unique_G = R2_GE - R2_E;
- unique_E = R2_GE - R2_G;
- shared = R2_G + R2_E - R2_GE.

These components are descriptive because pairwise rows are not independent. Inference is based on the matched vertex null for the partial-rank statistics.

## Interpretation ladder

1. rho_IBD supported only: geographic separation is sufficient to describe colour turnover under the tested environment.
2. rho_IBE supported only: colour turnover is more consistent with environmental sorting than simple distance.
3. both supported: both dispersal/history and environmental sorting contribute nonredundant structure.
4. rho_IBE supported but delta <= 0: environment adds information beyond geography, but is not stronger than the unique geographic component.
5. rho_IBE supported and delta > 0: under this BIO5 model, environmental sorting is stronger than unique geographic distance.

## Replication

The central IBE-like claim requires the same positive matched-null result in discovery and validation.

Third-cohort transport is reported as an additional species-disjoint test.

No result may be described as genetic IBE or local adaptation.

## Hard nonclaims

- colour phenotype is not a genetic distance;
- partial association does not establish causal selection;
- BIO5 can proxy correlated climate variables;
- phenotypic plasticity can generate the same pattern;
- drift/history and environmental sorting can coexist;
- local adaptation requires fitness or genotype-by-environment evidence.
