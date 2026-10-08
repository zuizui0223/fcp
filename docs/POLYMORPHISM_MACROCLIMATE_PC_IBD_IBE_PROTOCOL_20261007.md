# FCP macroclimate-PC IBD versus IBE-like protocol — 2026-10-07

## Motivation

The existing phenotypic IBD/IBE-like decomposition used BIO5 because warm-season temperature had already been singled out by an earlier signed white-versus-nonwhite analysis. That is too narrow for the broader evolutionary question: **does flower-colour turnover track the multivariate abiotic environment beyond geographic separation?**

This protocol therefore replaces the single BIO5 axis with a fixed, outcome-independent macroclimate space.

## Environmental variable set

Primary macroclimate space uses all reproducibly archived WorldClim 2.1 variables available from the frozen 10-minute resources:

- BIO1–BIO19;
- mean monthly solar radiation, averaged across the 12 archived SRAD layers.

These 20 variables jointly represent temperature level/seasonality/extremes, precipitation level/seasonality/extremes and incident solar radiation.

The analysis must be described as **macroclimate** or **macroabiotic environment**, not total environment. It does not include soil chemistry, UV-B directly, pollinator assemblages, herbivory or demographic/genetic history.

## PCA construction and transport

The environmental representation is frozen without using flower-colour outcomes.

1. Extract the 20 environmental variables for every classifiable discovery-cohort photo used by the spatial analyses.
2. Estimate variable means and standard deviations from discovery only.
3. Standardize discovery observations using those fixed moments.
4. Fit PCA to the standardized discovery environmental matrix by deterministic SVD.
5. Primary retained dimensionality is the smallest K for which cumulative explained variance is >=95%.
6. Freeze discovery means, SDs, loadings and K.
7. Apply the same transformation without refitting to validation and third-cohort observations.

Sensitivities use:
- the minimum K explaining >=80%;
- the minimum K explaining >=90%;
- full standardized 20-dimensional Euclidean environmental distance without PCA truncation.

No dimensionality is selected after seeing colour associations.

## Pairwise quantities within species

For every eligible species and every unordered photo pair:

- Y = continuous nine-colour flower Jensen–Shannon dissimilarity;
- G = great-circle geographic distance;
- E_PC = Euclidean distance in the frozen retained environmental PC space.

Species must retain at least 40 classifiable rows with valid coordinates, complete macroclimate values and valid continuous colour vectors.

## Symmetric partial-rank statistics

### Macroclimate IBE-like component

rho_ENV = partial Spearman(Y, E_PC | G)

### Geography / IBD-like component conditional on macroclimate

rho_GEO = partial Spearman(Y, G | E_PC)

### Relative balance

delta = rho_ENV - rho_GEO

The cohort statistic is the equal-species mean for each quantity.

## Matched vertex null

Within each species, permute complete continuous colour vectors among fixed observation vertices 199 times.

This preserves:
- coordinates;
- all macroclimate values;
- geography–environment covariance;
- the exact continuous colour-vector multiset;
- sample size and species-wide colour distribution.

Recalculate rho_ENV and rho_GEO for every null world.

Primary support for each component requires:
- positive observed equal-species mean;
- matched-null upper-tail p < 0.05.

Whether macroclimate exceeds geography is tested from delta using both:
- the matched permutation distribution of mean delta;
- a deterministic 9,999-replicate species-level sign-flip test.

## Comparison with the earlier BIO5-only analysis

Report side-by-side, without changing the older result:

- BIO5-only IBE-like mean partial rho;
- macroclimate-PC IBE-like mean partial rho;
- geography conditional on BIO5;
- geography conditional on macroclimate PC.

This comparison is descriptive. It asks whether a multivariate environmental representation captures more nonredundant flower-colour turnover than the previously selected single temperature axis.

## Replication rule

The macroclimate-PC IBE-like pattern is called replicated only if:
- discovery is positive with matched-null p < 0.05;
- validation is independently positive with matched-null p < 0.05.

Third cohort is a species-disjoint transport and is reported separately.

A third-cohort result cannot rescue failed discovery-validation replication.

## Interpretation ladder

1. **ENV unsupported, GEO supported**:
   flower-colour turnover is adequately described as spatial/history/dispersal structure under the measured macroclimate set.

2. **ENV supported, GEO supported, ENV < GEO**:
   macroclimate adds nonredundant colour sorting beyond distance, but geography/history remains the stronger component.

3. **ENV supported and ENV > GEO**:
   multivariate macroclimate explains a larger unique component than geographic distance under this phenotype-distance analysis.

4. **Macroclimate ENV > BIO5 ENV**:
   environmental sorting is better described as multivariate rather than a single warm-season temperature axis.

5. **BIO5 >= macroclimate ENV**:
   the previous BIO5 association is not strengthened by broadening to the full macroclimate space.

## Hard nonclaims

Even a replicated macroclimate-PC result does not establish:
- genetic isolation by environment;
- causal climate selection;
- local adaptation;
- fitness differences;
- genetic differentiation rather than plasticity;
- absence of unmeasured soil, pollinator, UV, demographic or historical drivers.

The strongest allowed wording is **phenotypic macroclimate-associated IBE-like sorting**.
