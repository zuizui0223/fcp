# disttrait

`disttrait` is a small Python package for **validated species-level distributional trait inference** from repeated individual observations.

It extracts the general inferential core used by the FCP flower-colour polymorphism study without hard-coding flower colours, iNaturalist, or the RGFCA world-map application.

## v0.12 scope

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
- partial-rank adjustment with matched spatial nulls;
- scalar continuous-trait spatial inference using absolute pairwise differences.
- multivariate continuous-trait spatial inference using Euclidean pairwise distances with complete-row permutations.

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

### Full-data numerical-identity boundary

The compact fixtures establish equivalence for the tested contracts, but `disttrait` should **not** be described as a bitwise numerical reproducer of every frozen FCP analysis. A post-confirmatory audit on the immutable third-cohort H2 artifact found:

- frozen study-specific W = **0.5172457461**;
- `disttrait.two_mode_axis` over the same frozen 158 species = **0.5183899565**.

The difference is small enough that the scientific decision is unchanged, but it is real. The generic implementation defensively normalizes input rows; the frozen FCP loader had already normalized them. For nearly tied farthest-point initializations, the extra floating-point normalization can change the deterministic argmax and lead two-means to a different local partition. The largest audited discrepancy was *Cirsium vulgare* (absolute axis cosine **0.6333**; frozen cluster sizes 5/45 versus renormalized 11/39); *Vicia benghalensis* showed a smaller near-tie shift (|cosine| **0.9830**; 5/44 versus 6/43).

Accordingly, the frozen study-specific pipeline and immutable biological result files control manuscript numbers. `disttrait` is a later general-purpose implementation of the same estimand/architecture, not the numerical provenance layer for this paper.

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

## License

`disttrait` is licensed under the **MIT License**. The package-local licence text is in `LICENSE`.

## Public-release metadata workflow

The monorepo intentionally does not guess software ownership metadata. A complete
template is provided at:

`release_metadata.example.toml`

For a public standalone release:

1. copy it to `release_metadata.toml`;
2. keep the frozen licence identifier `MIT`, and fill the ordered authors/maintainers and standalone
   repository URL;
3. place the selected full licence text in `LICENSE`;
4. set `release.release_ready = true`;
5. run `scripts/release/finalize_disttrait_release_metadata.py` from the FCP
   source repository;
6. run the public-release gate before tagging.

The finalizer renders `CITATION.cff`, author/maintainer/license metadata in
`pyproject.toml`, final repository URLs, and `RELEASE_METADATA.json`. It
refuses placeholders and does not choose a licence or authorship automatically.

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

v0.3 adds a broader 24-cell performance surface:

`results/disttrait_performance_surface_v0_3_20260919/`

It varies within-species effect size (0, 0.8, 1.6, 2.4), species-level observation imbalance (1× vs 4×) and MCAR observation loss (0%, 25%, 50%), with 40 replicate worlds per cell. Across the six null cells the maximum species-conditioned false-positive fraction is **0.025**, while naive pooled inference rejects in **100%** of null worlds under the deliberately strong between-species geographic-confounding design. Conditioned detection rises from **0.225–0.675** at effect 0.8 to **0.75–1.00** at effect 1.6 and **1.00** throughout at effect 2.4.

v0.4 adds an alternative-method comparator surface on the same synthetic worlds:

`results/disttrait_comparator_surface_v0_4_20260919/`

It compares naive pooling, equal-species matched-null inference, pair-count-weighted matched-null inference, and a one-sided one-sample test on species-specific rho values. Under the six null cells, maximum false-positive fractions are **0.025**, **0.025**, and **0.05** for the three species-conditioned summaries respectively, while naive pooling rejects in **100%** of null worlds. Under weak effect (0.8) plus 4× observation imbalance, equal-species matched-null detection spans **0.35–0.675**, compared with **0.325–0.575** for pair weighting and **0.35–0.70** for the simple species-rho t-test.

The interpretation is therefore bounded: species-conditioning is the major protection against between-species geographic confounding; equal species weighting can help under imbalance but is not uniquely optimal, and simple species-level tests can be competitive when their assumptions are adequate.

v0.5 adds a model-based species-conditioned comparator on the same 24 cells:

`results/disttrait_model_comparator_surface_v0_5_20260919/`

The added model uses species-specific intercepts with one common nonnegative within-species logistic slope and a one-sided profile-likelihood ratio test. Under the correctly specified binary-logistic data-generating process, weak-effect detection is **0.975–1.00**, compared with **0.225–0.675** for the equal-species matched-null omnibus. Across null cells, the fixed-effect logistic mean false-positive fraction is **0.05**, but its worst cell reaches **0.10**; the matched-null omnibus remains at or below **0.025** in every null cell.

The interpretation is a trade-off rather than a winner: a correctly specified response model can gain substantial power, while the matched-null approach is more assumption-light and has more stable worst-cell null behavior in this benchmark.

These benchmarks are targeted demonstrations, not evidence of universal superiority or robustness to missing-not-at-random trait observation.

## External non-flower empirical transport

v0.6 adds a public-domain empirical transport using a fixed snapshot of the San Francisco street-tree inventory, independent of FCP/RGFCA.

- source snapshot: TidyTuesday `2020-01-28/sf_trees.csv`;
- upstream source: San Francisco Public Works Street Tree List;
- source Git blob: `bdc06c1297b7dd88bea0df77de4007eaab30198e`;
- fixed eligible set: 52 taxon labels with >=500 usable DBH+coordinate records;
- outcome-blind deterministic selection: 20 taxon labels × 80 trees = 1,600 observations;
- trait: `log(DBH)`;
- within-taxon dissimilarity: absolute pairwise difference in `log(DBH)`;
- vertex-permutation nulls: 99 per taxon.

Frozen result:

`results/disttrait_sf_street_tree_empirical_v0_6_20260919/result.json`

The generic continuous-trait layer gives an equal-taxon mean spatial rho of **0.09219** with matched-null **p = 0.01**. Across taxa, the association between the SD of log(DBH) and spatial rho is not supported (**rho = -0.186, p = 0.83**). A naive pooled all-pairs analysis gives a small rho (**0.02783**) but an extremely small nominal p because it treats the very large number of pairs as ordinary independent evidence.

This empirical transport is deliberately non-causal. Urban planting, management, age structure, cultivar usage and inventory processes can all contribute to the observed DBH geography. Its role is to show that the package operates unchanged on an external, non-flower continuous trait dataset.

## Second external empirical transport: ShareTrait Gammarus

v0.9 adds a second external continuous-trait transport using the CC-BY-4.0 ShareTraitDatabase and the `Gammarus insensibilis` metabolic-rate dataset linked to Shokri et al. (2022; DOI `10.1242/jeb.244842`).

Frozen fixture:

`packages/disttrait/fixtures/sharetrait_gammarus_metabolic_v1.csv`

The frame contains **375 individuals** from three Adriatic populations with metabolic rate standardized to **Joule/day**. Individuals inherit their population-level site coordinates, so this is a between-population transport rather than fine-scale within-site spatial inference.

Primary analysis:

- trait = `log(metabolic_rate)`;
- dissimilarity = absolute pairwise log-rate difference;
- 999 within-species vertex permutations.

Frozen result:

`results/disttrait_gammarus_empirical_v0_9_20260919/result.json`

The result is intentionally **non-supporting**:

- log-rate rho = **-0.00683**, matched-null **p = 0.68**;
- raw-rate sensitivity rho = **0.01093**, matched-null **p = 0.198**.

This second transport is useful precisely because it does not reproduce the positive street-tree result. The package is being tested for executable generality across external datasets, not used to select only significant biological examples.

## Direction-heterogeneity / estimand-alignment benchmark

v0.7 adds a continuous-trait benchmark in which the within-species spatial effect may reverse sign among species.

Canonical script:

`packages/disttrait/benchmarks/direction_heterogeneity_surface.py`

Frozen result:

`results/disttrait_direction_heterogeneity_v0_7_20260919/result.json`

The design crosses:

- effect size: 0, 0.4, 0.8, 1.2;
- fraction of species with reversed signed effect: 0, 0.25, 0.5;
- MCAR missingness: 0%, 50%;
- 40 replicate worlds per cell.

When all species share the same direction, the common-slope fixed-effect linear model is more powerful at the weak effect: detection is **1.00** at effect 0.4, compared with **0.45–0.80** for the equal-species matched-null omnibus.

When half the species reverse direction, the common signed slope cancels. At effect 0.4 its detection falls to **0.05–0.10**, whereas the direction-invariant matched-null omnibus retains **0.475–0.90** detection. At effect 0.8 with 50% reversal, matched-null detection is **1.00** in both missingness cells while common-slope detection is only **0.05–0.125**.

This is an **estimand-alignment result**, not proof that the common-slope model is defective. A common signed slope is appropriate when a shared direction is the scientific target. Pairwise dissimilarity is appropriate when the target is spatial organization whose direction may differ among species.

## Permutation-calibrated species-slope meta-analysis

v0.8 adds species-specific signed slopes and a random-effects summary with a matched within-species permutation calibration.

Frozen benchmark:

`results/disttrait_random_slope_meta_v0_8_20260919/`

The first analytic random-effects preflight was deliberately **not** promoted: under 50% MCAR missingness, asymptotic heterogeneity and omnibus false-positive fractions reached **0.30**. Re-estimating species slopes and variances inside 99 matched within-species trait-permutation worlds reduced the calibrated combined omnibus maximum null rejection to **0.05**.

Under weak effect (0.4) with 50% directional reversal:

- common signed-slope detection: **0.025–0.05**;
- direction-invariant matched-null detection: **0.40–0.90**;
- calibrated slope-heterogeneity detection: **0.75–1.00**;
- calibrated slope-meta omnibus detection: **0.525–1.00**.

The interpretation is again estimand-specific. Common slopes target a shared signed response; random-effects slope meta-analysis separates average signed response from directional heterogeneity; pairwise distance-dissimilarity targets spatial organization without requiring a common sign.

## MNAR observation-process stress test

v0.10 adds a latent-null observation-process benchmark:

`results/disttrait_mnar_observation_v0_10_20260919/`

The latent biological process contains **no within-species trait-position
association**. Before analysis, rows are selected under four mechanisms:

- MCAR;
- trait-only selection;
- position-only selection;
- joint trait-by-position selection.

The benchmark uses 20 species, 100 latent observations/species, 40 worlds per
cell and two selection strengths.

All 8 cells retain **40/40 evaluable worlds**.

Under MCAR, trait-only and position-only selection:

- maximum equal-species matched-null rejection = **0.05**;
- observed species-level spatial rho remains near zero.

Under joint trait-by-position selection:

- matched-null rejection = **0.80–1.00**;
- species fixed-effect logistic rejection = **1.00**;
- mean absolute observed state-position rho rises to **0.21–0.32**.

This is not interpreted as an ordinary type-I calibration failure conditional
on the observed sample. Joint selection has changed the observed sample so that
trait and position are genuinely associated among retained rows.

The supported boundary is therefore explicit:

> matched within-species randomization calibrates inference conditional on the
> observed sampling geometry and observed trait multiset; it cannot, by itself,
> distinguish biological spatial organization from unmeasured observation
> processes that depend jointly on trait and location.

This is why v0.10 does **not** claim robustness to arbitrary MNAR sampling.

## Nonlinear within-species benchmark

v0.11 adds a continuous-trait curvature benchmark:

`results/disttrait_nonlinear_curvature_v0_11_20260919/`

The data-generating process uses a centered quadratic local-position effect
within species. Curvature sign is shared by all species, reversed in 25%, or
reversed in 50% of species.

Compared methods:

- naive pooled pairwise analysis;
- equal-species matched-null distance-dissimilarity;
- species fixed-intercept common linear slope;
- species fixed-intercept common quadratic curvature.

Under a weak shared curvature effect (0.4):

- correctly specified quadratic detection = **1.00**;
- matched-null detection = **0.10–0.375**;
- common linear detection = **0–0.025**.

When 50% of species reverse curvature sign:

- weak-effect common quadratic detection falls to **0.025–0.05**;
- matched-null detection remains **0.15–0.35**;
- at effect 0.8, matched-null detection is **0.85–1.00** while common
  quadratic detection is **0–0.075**.

The result is an estimand-alignment result rather than a method ranking.
A correctly specified common nonlinear model is highly efficient when species
share curvature direction. Distance-dissimilarity remains informative when
species are spatially organized but the sign of nonlinear curvature differs.

Across the six null cells, the largest observed rejection fractions in 40
worlds are **0.075** for matched-null, **0.10** for the common linear model and
**0.075** for the common quadratic model. These finite simulation frequencies
are not treated as universal calibration guarantees.

## Multivariate continuous-trait orientation benchmark

v0.12 adds complete-row multivariate spatial inference and a 24-cell
orientation-heterogeneity benchmark:

`results/disttrait_multivariate_orientation_v0_12_20260919/`

The benchmark uses a two-dimensional continuous trait with strong
between-species geographic turnover and a within-species response vector whose
orientation spans 0, half, or all of the trait-space circle.

Across the six effect-zero cells:

- maximum matched-null rejection = **0.05**;
- maximum common-vector rejection = **0.05**;
- naive pooled rejection = **1.00** in every null cell.

With effect size 0.4 and a shared response orientation:

- matched-null detection = **0.975–1.00**;
- common-vector detection = **1.00**.

With the same effect size but response orientations spread across the full
trait-space circle:

- matched-null detection = **0.95–1.00**;
- common-vector detection = **0.075**.

At effect size 0.8 under full orientation spread:

- matched-null detection = **1.00**;
- common-vector detection = **0.10–0.25**.

The interpretation is estimand-specific. A common signed multivariate vector is
appropriate when species share response orientation. Euclidean
distance-dissimilarity asks whether individuals become more different in trait
space with geographic separation and therefore remains informative when
species-specific response vectors cancel in the cross-species mean.

Complete multivariate rows are the permutation unit; covariance among trait
dimensions is preserved within observations. The package does not rescale
dimensions automatically, so biologically justified standardization remains an
application responsibility.

## Validation gate

The release-candidate gate is the dedicated `disttrait package` workflow. It installs the package from `packages/disttrait`, compiles the public modules, and runs the full package test suite including compact frozen-FCP equivalence fixtures, non-flower generalization tests, and the repeated species-conditioning benchmark.

## Release status

The validated in-repository package version is **0.12.0**.

Release engineering files:

- `CHANGELOG.md` — version history;
- `RELEASING.md` — external-release checklist and unresolved ownership choices;
- `.github/workflows/disttrait-release-check.yml` — sdist/wheel build, Twine
  validation and fresh-environment install smoke test.

The package licence, final author/maintainer metadata and standalone repository
are intentionally not guessed here and remain explicit release decisions.

## Relationship to FCP

The active New Phytologist manuscript remains an ecological application. `disttrait` is the reusable methods layer.

Flower-colour-specific objects such as the four frozen colour states, ROI-v4/EfficientSAM measurement, and the white-versus-nonwhite `q_white` target remain outside the generic package.
