# disttrait

`disttrait` is a small Python package for **validated species-level distributional trait inference** from repeated individual observations.

It extracts the general inferential core used by the FCP flower-colour polymorphism study without hard-coding flower colours, iNaturalist, or the RGFCA world-map application.

## v0.2 scope

The package provides reusable components for:

- species-level categorical diversity using Gini-Simpson diversity;
- observer-disjoint split-half reliability;
- Hellinger-transformed two-mode trait geometry;
- fixed contrast / axis-alignment statistics;
- construction-preserving row permutation within fixed strata and full fixed-contrast structured-null refitting;
- pairwise great-circle distances;
- pairwise Jensen-Shannon trait dissimilarity;
- species-specific spatial organization `rho_i`;
- matched trait-minus-background spatial organization and joint same-observation permutation nulls;
- species-level distribution-versus-spatial association;
- partial-rank adjustment with matched spatial nulls.

The package deliberately does **not** claim a new standalone statistic. Most components are established statistics. The reusable contribution is an inference architecture that separates measurement validity, species-level trait distributions, trait geometry, spatial organization, matched nulls, and prospective confirmation.


## FCP equivalence status

v0.2 adds compact exact-equivalence fixtures against the frozen FCP/RGFCA implementations for:

- observer-disjoint partitioning, classifiability and split-half reliability;
- deterministic Hellinger two-means and mode orientation;
- the generic one-vs-rest contrast versus the frozen nine-dimensional `q_white`;
- species-specific distance-trait vertex-permutation nulls;
- matched focal-minus-background joint-permutation nulls.

Reference source Git blobs and result contracts are frozen in:

`packages/disttrait/fixtures/fcp_equivalence_manifest.json`

The clean submission branch does not include the original 100,000-row discovery/reserve measurement tables. Therefore the always-on tests establish **algorithmic equivalence**, not a complete raw-data replay of every biological result.

For a full H1 replay when the original measured tables are available:

```bash
python packages/disttrait/scripts/replay_fcp_h1.py \
  --discovery /path/to/global_monte_carlo_measured_photos_v1.csv \
  --reserve /path/to/rgfca_reserve_replication_measured_photos_v1.csv \
  --expected results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json
```

For an H2 observed-W replay when the four legacy delta-vector artifacts are available:

```bash
python packages/disttrait/scripts/replay_fcp_h2_observed_w.py \
  --primary-discovery /path/to/primary_0_10_discovery_delta_vectors.csv \
  --primary-reserve /path/to/primary_0_10_reserve_delta_vectors.csv \
  --strict-discovery /path/to/strict_0_20_discovery_delta_vectors.csv \
  --strict-reserve /path/to/strict_0_20_reserve_delta_vectors.csv \
  --expected results/polymorphism_white_axis_targeted_test_20260912/result.json
```

The H2 replay checks species denominators and the generic fixed-contrast alignment statistic. It does not reconstruct the 999 construction-preserving null worlds unless the underlying row-level palette artifacts are also supplied.

## Install from this repository

```bash
python -m pip install -e packages/disttrait
```

## Small API example

```python
import numpy as np
from disttrait import gini_simpson, spatial_rho

D = gini_simpson([80, 20])

rho = spatial_rho(
    latitude=[35.0, 35.5, 36.0, 36.5],
    longitude=[135.0, 135.5, 136.0, 136.5],
    traits=np.array([
        [1.0, 0.0],
        [0.8, 0.2],
        [0.2, 0.8],
        [0.0, 1.0],
    ]),
)
```

## Non-flower demonstration

`packages/disttrait/examples/nonflower_binary_trait.py` uses an abstract binary trait with no flower-colour assumptions.

Two deterministic examples are frozen:

- **signal recovery:** induced species-level diversity–spatial association, rho = **0.812**;
- **pooled confounding:** naive pooled geographic analysis gives rho = **0.578** even though within-species spatial dependence is absent by construction, while the species-conditioned mean rho is **0.005**.

Frozen receipt:

`results/disttrait_nonflower_synthetic_v0_2_20260918/result.json`

This is a synthetic demonstration, not a general superiority benchmark.

A repeated 40-null / 40-signal benchmark targeting between-species geographic confounding is also frozen in:

`results/disttrait_species_conditioning_benchmark_v0_2_20260918/result.json`

Under that specific failure mode:

- naive pooled false-positive fraction: **1.00**;
- species-conditioned matched-null false-positive fraction: **0.00**;
- species-conditioned signal detection fraction: **1.00**.

The benchmark is intentionally narrow and does not establish universal superiority.

## Validation gate for v0.2

The release-candidate gate is the dedicated `disttrait package` workflow. It installs the package from `packages/disttrait`, compiles the public modules, and runs the full package test suite including compact frozen-FCP equivalence fixtures, non-flower generalization tests, and the repeated species-conditioning benchmark.

## Relationship to FCP

The active New Phytologist manuscript remains an ecological application. `disttrait` is the reusable methods layer.

Flower-colour-specific objects such as the four frozen colour states, ROI-v4/EfficientSAM measurement, and the white-versus-nonwhite `q_white` target remain outside the generic package.
