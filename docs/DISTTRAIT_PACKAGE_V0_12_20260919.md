# disttrait v0.12 package status — 2026-09-19

Status: validated candidate with multivariate continuous-trait orientation benchmark.

Package: `packages/disttrait/`

Candidate version: `0.12.0`

## 1. New in v0.12

v0.12 extends the continuous-trait spatial layer from scalar traits to complete
multivariate trait vectors.

New public API:

- `euclidean_pairwise`
- `multivariate_spatial_rho`
- `multivariate_spatial_permutation_null`

Rows are individual observations and columns are continuous trait dimensions.
The default dissimilarity is Euclidean distance. The permutation unit is the
**complete multivariate row**, so covariance/dependence among trait dimensions
is preserved within every observation.

The package deliberately does not rescale trait dimensions automatically.
Applications must standardize or otherwise define the metric upstream when
dimensions are not commensurate.

## 2. Benchmark design

Canonical benchmark:

`packages/disttrait/benchmarks/multivariate_orientation_surface.py`

Frozen outputs:

- `results/disttrait_multivariate_orientation_v0_12_20260919/cells.csv`
- `results/disttrait_multivariate_orientation_v0_12_20260919/result.json`

The benchmark uses two-dimensional continuous traits and crosses:

- effect size: **0, 0.4, 0.8, 1.2**;
- response-orientation spread: **0, 0.5, 1.0** of the full trait-space circle;
- MCAR observation loss: **0%, 50%**;
- 40 replicate worlds per cell;
- 20 species per world;
- 20 observations per species before missingness;
- 39 matched spatial-null permutations.

There are 24 cells.

Every species also has a baseline trait vector linked to its geographic centre,
creating deliberate between-species geographic confounding.

## 3. Compared estimands

### Naive pooled pairwise analysis

All species are pooled before relating geographic distance to Euclidean
multivariate trait distance.

### Equal-species matched-null trait-space organization

Within each species:

`rho_i = Spearman(pairwise geographic distance, pairwise Euclidean trait distance)`.

Species contribute equally and complete multivariate trait rows are permuted
among observed positions within species.

This estimand asks whether individuals become more different in multivariate
trait space with geographic separation, without requiring species to share one
signed response vector.

### Common signed multivariate response vector

A species fixed-intercept multivariate linear model estimates one common
two-dimensional signed response vector.

This estimand is efficient when species share orientation. As response
directions spread around trait space, the common vector can cancel even when
every species is strongly organized.

## 4. Null behavior

Across all six effect-zero cells:

- minimum naive pooled rejection = **1.00**;
- maximum equal-species matched-null rejection = **0.05**;
- maximum common-vector rejection = **0.05**.

These are finite 40-world simulation frequencies rather than universal
calibration guarantees.

## 5. Shared-orientation result

At effect size 0.4 and orientation spread 0:

- common-vector detection = **1.00**;
- matched-null detection = **0.975–1.00**.

When the shared response vector is the scientific target, the common-vector
model is fully competitive.

## 6. Full orientation-spread result

At effect size 0.4 and orientation spread 1.0:

- matched-null detection = **0.95–1.00**;
- common-vector detection = **0.075**.

At effect size 0.8 and orientation spread 1.0:

- matched-null detection = **1.00**;
- common-vector detection = **0.10–0.25**.

The common response-vector norm approaches zero because species directions span
the full trait-space circle and cancel in the cross-species mean, despite strong
within-species organization.

## 7. Interpretation

v0.12 closes the multivariate continuous-trait benchmark gap.

The supported interpretation is:

> Complete-row multivariate matched-null inference can quantify
> within-species trait-space organization without assuming a shared signed
> response orientation across species. A common multivariate response vector is
> efficient when orientations are shared, but answers a different question and
> can cancel when equally strong species responses point in different trait-space
> directions.

This is an estimand-alignment result rather than a universal method ranking.

## 8. Computational implementation

The multivariate vertex null uses an exact rank-matrix optimization.

A complete-row vertex permutation preserves the full multiset of pairwise
Euclidean trait distances. Pairwise trait-distance ranks are therefore computed
once and re-indexed under each vertex permutation rather than re-ranked from
scratch.

A dedicated equivalence test confirms that the optimized null is numerically
identical to direct re-ranking on a fixed fixture.

## 9. Claim boundary

v0.12 does not establish that:

- Euclidean distance is optimal for every multivariate phenotype;
- trait dimensions can be combined without biologically justified scaling;
- a common multivariate response model is inappropriate when species share
  orientation;
- matched-null inference identifies which trait dimension drives organization;
- the finite null rejection rates are universal type-I error guarantees;
- the benchmark covers arbitrary high-dimensional, nonlinear or joint-MNAR
  multivariate processes.

## 10. Validation

Dedicated workflow:

`.github/workflows/disttrait-multivariate-orientation.yml`

Exploratory frozen-result run:

- PR #53;
- head: `0595fe35254d0b1e7336a6611ec690b17f07b388`;
- run: `35412417908`;
- job: `105814495686`;
- artifact: `10573828639`;
- conclusion: success.

The final workflow regenerates all 24 cells and verifies them against the frozen
v0.12 receipt.

## 11. Methods-paper state after v0.12

The validation stack now includes:

1. exact FCP/RGFCA implementation-equivalence fixtures;
2. pooled-confounding benchmarks;
3. effect/imbalance/MCAR performance surfaces;
4. nonparametric and model-based comparators;
5. signed direction heterogeneity;
6. permutation-calibrated species-slope meta-analysis;
7. two external non-flower empirical transports;
8. explicit MNAR observation-process limits;
9. nonlinear curvature / shape-heterogeneity benchmarks;
10. **multivariate continuous-trait orientation heterogeneity**.

The previously listed nonlinear and multivariate scientific benchmark gaps are
now closed for the tested families.

The main remaining pre-release task is standalone package infrastructure:
repository separation, licence, citation metadata and versioned release
automation.
