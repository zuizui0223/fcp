# disttrait v0.6 candidate — continuous traits and external empirical transport

Date: 2026-09-19 JST

Status: validation candidate; empirical numerical result not yet frozen.

## 1. Purpose

Extend the reusable spatial-inference layer beyond compositional traits and Jensen-Shannon divergence.

v0.6 candidate adds:

- scalar continuous-trait spatial organization using absolute pairwise trait difference;
- arbitrary symmetric observation-by-observation trait-distance matrices;
- matched vertex-permutation nulls for those generic distances;
- an external empirical transport workflow using a public Ricinus communis individual-trait dataset.

## 2. New generic APIs

- `pairwise_absolute_difference(values)`
- `spatial_rho_from_distance_matrix(latitude, longitude, trait_distance_matrix)`
- `spatial_permutation_null_from_distance_matrix(...)`
- `continuous_spatial_rho(latitude, longitude, values)`
- `continuous_spatial_permutation_null(...)`

For a scalar trait y, the default continuous-trait statistic is:

`rho = Spearman(pairwise great-circle distance, |y_i - y_j|)`.

The matched null keeps the observed coordinates and complete trait values fixed, then permutes individual trait assignments among geographic observations.

## 3. External empirical transport source

The transport example uses an external dataset from:

- repository: `The-Frederickson-Lab/castor_latitudinal_gradient`
- pinned commit: `ff8cb26f83cd271c516e83510cd06c4b5563db22`
- source file: `morph_natint_20May24.csv`
- source Git blob: `c4dcca178d7c885585edbe0cdbeee89eef2a484a`

The original study sampled Ricinus communis individuals across a broad latitudinal range and records individual geographic coordinates and multiple morphological/ecological traits.

The source data are **not copied into fcp**. The example reads the commit-pinned public file at execution time.

Canonical transport script:

`packages/disttrait/examples/castor_empirical_transport.py`

Dedicated networked validation workflow:

`.github/workflows/disttrait-empirical-transport.yml`

## 4. Empirical traits

The transport script evaluates three continuous individual-level traits:

- leaf-area proxy: `areaavg`;
- herbivory fraction: `pherb`;
- leaf-blade EFN density: `EFNLB`.

For each trait it reports:

- complete georeferenced record count;
- observed species-level spatial rho;
- 999-permutation upper-tail p;
- null mean, median and 95% interval.

These are transport diagnostics rather than predeclared biological hypotheses.

## 5. Claim boundary

A successful run would support only that the generic continuous-trait API executes reproducibly on an independent empirical dataset with a different phenotype representation.

It would not establish:

- cross-species comparative validity, because this dataset is one focal species;
- causal ecological interpretation of any spatial association;
- general calibration under the external study's sampling process;
- superiority over the original study's statistical models;
- independence from all community-science or geographic sampling biases.

## 6. Promotion gate

Do not promote to v0.6.0 until both are green:

1. the hermetic `disttrait package` workflow, including continuous-trait unit tests;
2. the networked `disttrait empirical transport` workflow against the pinned external commit.

After the first successful empirical workflow, freeze its JSON receipt in the repository and convert the workflow into a regression check against that receipt.
