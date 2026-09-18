# Methods classification for the New Phytologist polymorphism paper — 2026-09-18

Status: reporting-only methods audit. No biological analysis is changed.

## Purpose

Separate standard statistical components from study-specific design choices in the active flower-colour polymorphism paper.

## Standard or widely used components

The paper uses several established methods without claiming methodological novelty for them individually:

- Gini–Simpson diversity: `D = 1 - sum p_k^2`;
- Spearman rank correlation;
- Lin-style concordance correlation coefficient;
- Spearman–Brown reliability projection;
- Hellinger transformation of compositional data;
- k-means/two-means clustering;
- Jensen–Shannon divergence for compositional colour dissimilarity;
- great-circle geographic distance;
- permutation/randomization tests with plus-one Monte Carlo p-values;
- partial-rank residualization;
- Blomberg's K;
- Pagel's lambda;
- phylogenetic generalized least-squares sensitivity analyses.

These methods are not presented as new.

## Study-specific design choices

The methodological contribution lies in the way established components are assembled around opportunistic image data and species-level trait distributions:

1. fixed high-depth photo budgets per species rather than effort-proportional weighting;
2. observer contribution caps and geographic maximin photo selection;
3. location-blind image measurement;
4. explicit terminal distinction among classifiable colour states, clear technical failures and ambiguous palette compositions;
5. observer-disjoint reliability testing of a species-level distributional phenotype before ecological interpretation;
6. separation of target discovery from a new species/photo-disjoint prospective confirmation cohort;
7. a construction-preserving H2 null that keeps every species × coarse-morph row count while permuting continuous palette rows within coarse states;
8. reuse of matched, species-conditioned spatial nulls that preserve observed sampling geometry;
9. matched flower-minus-background spatial diagnostics;
10. one-shot/no-rerun and no-replacement prospective execution rules.

## What is generalizable

The generalizable workflow is:

[
	ext{opportunistic individual observations}
ightarrow
	ext{validated species-level trait distribution}
ightarrow
	ext{distributional geometry}
ightarrow
	ext{species-specific spatial organization}
ightarrow
	ext{cross-species comparison}.
]

This framework is not intrinsically restricted to flower colour. In principle it applies to georeferenced within-species phenotypes for which repeated individual-level observations and a defensible measurement representation exist.

## What is flower-colour specific

The following parts are application-specific:

- the four biological colour states;
- the nine-colour palette;
- the fixed white-versus-nonwhite contrast `q_white`;
- ROI-v4/EfficientSAM flower localization;
- interpretation of flower/background colour structure.

## Claim boundary

The paper should not claim a wholly new statistical method. Its methodological contribution is a rigorous validation and inference architecture for species-level intraspecific trait distributions from heterogeneous community-science images, assembled from established statistics plus study-specific sampling, null-model and prospective-confirmation rules.
