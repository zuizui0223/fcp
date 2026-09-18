# disttrait v0.1 package design — 2026-09-18

Status: reusable methods package extracted from the FCP/New Phytologist analysis architecture.

## Goal

`disttrait` generalizes the species-level distributional-trait inference used by the active flower-colour polymorphism paper without encoding flower-colour biology.

The package lives at:

`packages/disttrait/`

and is intentionally isolated from the historical `fcp_pipeline` package.

## v0.1 public API

### Diversity

- `gini_simpson(values)`
- `categorical_diversity(labels, categories=...)`

These convert repeated categorical observations into a continuous species-level diversity phenotype.

### Observer-disjoint reliability

- `observer_disjoint_reliability(...)`

This reproduces the FCP design principle of keeping each observer intact, balancing complete observer groups between halves without reading outcomes, calculating the species-level diversity phenotype independently in each half, and summarizing split-half Spearman/CCC/Spearman-Brown diagnostics.

The function is generic in:

- species column;
- observer column;
- categorical state column;
- allowed state set;
- minimum full/half information gates;
- number of deterministic partitions and seed.

### Distributional geometry

- `hellinger_rows(...)`
- `two_mode_axis(...)`
- `one_vs_rest_contrast(...)`
- `alignment_statistic(...)`
- `permute_rows_within_strata(...)`

This layer separates continuous trait geometry from categorical diversity. It supports the general version of the H2 logic: infer a species-specific displacement direction, compare sign-invariant axes with a fixed contrast, and preserve a construction through row permutations within fixed strata.

The FCP-specific nine-colour palette and white/non-white contrast are not hard-coded.

### Spatial organization

- `great_circle_pairwise_km(...)`
- `jensen_shannon_pairwise(...)`
- `spatial_rho(...)`
- `spatial_permutation_null(...)`
- `matched_difference_spatial_rho(...)`

This implements the species-specific spatial statistic used upstream in RGFCA/FCP:

`rho_i = Spearman(pairwise geographic distance, pairwise trait JSD)`.

The matched-background function uses the exact form:

`Spearman(distance, focal JSD - matched-background JSD)`

rather than subtraction of two separate correlation coefficients.

### Cross-species association

- `distribution_spatial_association(...)`
- `partial_rank_correlation(...)`
- `partial_distribution_spatial_association(...)`

These implement the threshold-free D-spatial logic and its matched spatial-null / rank-residualized controls.

## What is general versus application-specific

### General package layer

- repeated individual observations;
- categorical diversity;
- observer-disjoint validation;
- compositional geometry;
- fixed-contrast alignment;
- construction-preserving row permutations;
- species-specific spatial organization;
- matched spatial controls;
- across-species distribution-spatial association.

### Remains in FCP application code

- iNaturalist acquisition;
- flowering-annotation filters;
- ROI-v4 and EfficientSAM;
- four named flower-colour states;
- nine named flower-colour palette coordinates;
- the frozen white/non-white `q_white`;
- H3 phylogenetic analyses;
- prospective workflow authorization / no-rerun infrastructure.

## Methodological claim

The package does not claim invention of Gini-Simpson diversity, Spearman correlation, Hellinger transforms, k-means, Jensen-Shannon divergence or permutation tests.

The reusable contribution is the **inference architecture**:

[
	ext{repeated observations}
ightarrow
	ext{validated species-level distribution}
ightarrow
	ext{distributional geometry}
ightarrow
	ext{species-specific spatial organization}
ightarrow
	ext{cross-species inference}.
]

## Release path

v0.1 is an in-repository standalone package.

A later standalone repository/PyPI release should only occur after:

1. exact equivalence tests against frozen FCP H1/H2/spatial fixtures;
2. API naming review;
3. versioned documentation/examples;
4. package license selection;
5. semantic versioning/release automation;
6. at least one non-flower synthetic or empirical demonstration.

## Verification

Dedicated CI:

`.github/workflows/disttrait-package.yml`

The v0.1 test suite checks:

- Gini-Simpson calculations;
- deterministic Hellinger/two-mode geometry;
- one-vs-rest contrast recovery;
- stratum-preserving row permutation;
- observer-disjoint split-half reliability;
- pairwise geographic/JSD geometry;
- deterministic spatial nulls;
- matched focal-minus-background definition;
- distribution-spatial and partial-rank association.
