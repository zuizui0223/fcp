# disttrait v0.2 package design and FCP equivalence — 2026-09-18

Status: active in-repository methods package.

Package:

`packages/disttrait/`

Current version:

`0.2.0`

## 1. Purpose

`disttrait` extracts the reusable inference architecture underlying the active FCP/New Phytologist flower-colour polymorphism study.

It does not encode flower-colour biology. Its target is the more general problem:

> How can repeated, heterogeneous individual observations be converted into validated species-level trait distributions, their geometry and their species-specific spatial organization?

## 2. Public methodological layers

### Layer A — distribution amount

- `gini_simpson`
- `categorical_diversity`

A repeated categorical phenotype is converted into a continuous species-level diversity value.

### Layer B — measurement/reliability validation

- `observer_disjoint_reliability`

Observers remain intact and are assigned outcome-blind to two halves. The species phenotype is recalculated independently in both halves and compared across species using split-half Spearman correlation, CCC and Spearman-Brown diagnostics.

v0.2 matches the frozen FCP H1 implementation in:

- blank/whitespace observer exclusion;
- optional explicit classifiability column;
- species/observer string normalization;
- full-row observer balancing;
- exact SHA256 tie-break string `seed | species | observer`;
- category-count construction;
- split seed schedule.

### Layer C — distributional geometry

- `hellinger_rows`
- `two_mode_axis`
- `one_vs_rest_contrast`
- `alignment_statistic`
- `permute_rows_within_strata`

v0.2 matches the frozen FCP H2 deterministic two-means construction in:

- row normalization;
- Hellinger transform;
- grand-mean farthest-point initialization;
- second center chosen farthest from the first;
- deterministic empty-cluster repair;
- maximum 200 iterations;
- major/minor cluster orientation;
- lexicographic tie orientation;
- separation-ratio definition.

The generic package does not hard-code the FCP nine-colour palette or white/non-white target.

### Layer D — species-specific spatial organization

- `great_circle_pairwise_km`
- `jensen_shannon_pairwise`
- `spatial_rho`
- `spatial_permutation_null`

The core statistic is:

`rho_i = Spearman(pairwise geographic distance, pairwise trait Jensen-Shannon divergence)`.

The vertex-permutation null preserves coordinates, species identity and the complete set of individual trait rows while breaking their assignment to positions.

### Layer E — matched-background / matched-control spatial structure

- `matched_difference_spatial_rho`
- `matched_difference_spatial_permutation_null`

For paired focal and control measurements from the same observations:

`rho_i = Spearman(distance, focal dissimilarity - matched-control dissimilarity)`.

The null applies one common vertex permutation to the paired focal/control observations. This matches the joint same-photo reserve construction used by FCP/RGFCA.

### Layer F — cross-species association

- `distribution_spatial_association`
- `partial_rank_correlation`
- `partial_distribution_spatial_association`

These implement the general form of the FCP D-spatial analysis, including reuse of matched species-level spatial null realizations and rank-residualized controls.

## 3. Exact-equivalence evidence

Manifest:

`packages/disttrait/fixtures/fcp_equivalence_manifest.json`

Always-on tests:

`packages/disttrait/tests/test_fcp_equivalence.py`

They compare the generic package against compact fixtures generated from the frozen legacy FCP/RGFCA algorithms for:

1. observer-disjoint reliability;
2. deterministic Hellinger two-means;
3. the generic one-vs-rest contrast versus the frozen q_white loadings;
4. the species-specific vertex-permutation spatial null;
5. the matched focal-minus-background joint-permutation null.

The manifest records exact source branches, file paths and Git blob SHAs for the legacy reference implementations.

## 4. Full raw-data replay boundary

The clean submission main branch intentionally does not contain the original 100,000-row discovery/reserve measured-photo tables.

Therefore the always-on equivalence suite establishes **algorithmic equivalence** of the portable methods layer, not a full raw-data replay of every biological result.

When the original artifacts are available, full H1 replay is supported by:

`packages/disttrait/scripts/replay_fcp_h1.py`

Observed H2 fixed-axis alignment can be replayed from the legacy delta-vector artifacts with:

`packages/disttrait/scripts/replay_fcp_h2_observed_w.py`

Example:

```bash
python packages/disttrait/scripts/replay_fcp_h1.py \
  --discovery /path/to/global_monte_carlo_measured_photos_v1.csv \
  --reserve /path/to/rgfca_reserve_replication_measured_photos_v1.csv \
  --expected results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json
```

This recalculates the 200 observer-disjoint partitions from raw rows and checks the main frozen H1 summary statistics. The H2 replay checks the four primary/strict discovery/reserve species denominators and observed W values against the frozen targeted-result JSON. Full H2 structured-null replay still requires the row-level palette artifacts.

## 5. What remains application-specific

The following remain outside `disttrait`:

- iNaturalist acquisition/API logic;
- RGFCA equal-area global discovery;
- observer/photo availability census;
- flower localization and segmentation;
- the four named flower-colour states;
- the named nine-colour palette;
- the frozen q_white biological target;
- phylogeny-specific H3 analyses;
- workflow-level no-rerun authorization and chain-of-custody controls.

Those can be supplied by an application layer around the generic inference package.

## 6. Methodological claim ceiling

`disttrait` does not introduce a new diversity index, rank correlation, clustering algorithm, dissimilarity measure or permutation principle.

The defensible methodological contribution is:

> a reusable inference architecture that separates species-level distribution construction, measurement validity, trait geometry, species-specific spatial organization, matched-control nulls and cross-species association.

The architecture is intended for repeated individual-level traits beyond flower colour, provided the application supplies a defensible measurement representation and observation provenance.

## 7. Current verification layers

1. **unit tests** — basic mathematical behavior;
2. **synthetic integration tests** — reliability and spatial associations;
3. **frozen FCP algorithm-equivalence tests** — exact compact fixtures;
4. **optional raw-artifact replay** — full H1 result replay when source measurement tables are supplied;
5. **dedicated package CI** — `.github/workflows/disttrait-package.yml`.

## 8. Non-flower generalization demonstration

Executable example:

`packages/disttrait/examples/nonflower_binary_trait.py`

Frozen result:

`results/disttrait_nonflower_synthetic_v0_2_20260918/result.json`

The example contains no flower-colour representation.

### Induced-signal world

Across 80 synthetic species, an abstract binary-trait spatial gradient is constructed so that species-level diversity and latent spatial strength covary. The generic package recovers:

- rho(diversity, species-specific spatial organization) = **0.8115**;
- rho(latent spatial strength, recovered spatial organization) = **0.8382**.

### Between-species geographic-confounding world

Across 40 species, baseline trait frequency changes strongly among geographic species centres, but trait state is spatially independent **within each species** by construction.

- naive pooled distance–trait correlation = **0.5778**;
- species-conditioned mean rho = **0.0053**;
- species-conditioned median rho = **0.0000**.

This demonstrates the estimand distinction that motivated RGFCA/disttrait: between-species geographic turnover should not be interpreted as within-species spatial organization.

The demonstration does not establish calibrated superiority across arbitrary data-generating processes.

## 9. Next methods-paper gate

Before claiming an independent methods paper/package release, add:

1. full artifact-backed equivalence for H1 and selected H2/spatial fixtures;
2. simulation benchmarks against simpler alternatives;
3. at least one non-flower demonstration;
4. calibrated false-positive/power experiments under observation imbalance and missingness;
5. standalone repository/package licence and release automation.

At v0.2, the package is a verified reusable implementation layer, not yet a separately validated general methods publication.
