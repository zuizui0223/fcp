# Hypervolume transfer: what separates environmental sharing from individual structure?

Date: 2026-09-07 (Japan). Status: **synthetic benchmark complete; empirical transfer inference not opened**.

## Answer first

Colour-conditioned Gaussian hypervolumes become informative about sharing when their partitions transfer to different species. Total volume is colour-blind; within-species separation alone does not distinguish shared from independent species-specific structure.

Under the deliberately constructed, moderate-signal synthetic example, environmental-only recovery was **472/500 = 94.4%** (Wilson 95% interval 92.03–96.10%). This is a method-development result, not evidence that real flowers share a climatic boundary.

The specification was committed before simulation at `bb23f879f6d84a27d0c7fdbde430b9d4f4b4f7ad`. The branch was forked from parent snapshot `fe336db8fe9344608886a421ac4a0e6fc6aa38d2`. No frozen G1, environmental mean-effect, interaction, or species-disjoint support decision was changed.

## 1. What was run

The benchmark contains 40,500 simulated worlds: 24,000 calibration worlds and 16,500 evaluation worlds. Each world contains 60 species with 20 synthetic photos each. The first 30 species train the colour-conditioned Gaussian models; the other 30 evaluate transfer. There are 500 evaluation replicates per scenario. The 900 dependent training/evaluation species pairs within one world are not treated as 900 independent replicates.

The two representations are a three-dimensional Cartesian embedding of a sphere and a synthetic two-dimensional environmental field. The latter is a smooth, even fourth-harmonic function of spherical coordinates, intentionally constructed to have repeating environmental conditions that are not recoverable from the same low-order geographic Gaussian summaries. It is **not** temperature, VPD, CHELSA, or a fitted ecological environment.

The worlds include shared geographic structure, shared environmental structure, independent geographic species structure, independent environmental species structure, a mixture of independent structures, and unstructured noise. A separate confounded control sets environmental coordinates equal to two geographic coordinates and plants their common signal.

For each species and synthetic colour, covariance is estimated with fixed 20% isotropic shrinkage. Models predict partitions in held-out species using equal-prior Gaussian discriminant scores. Adjusted Rand index (ARI) measures agreement without treating colour-name reversal as disagreement. No held-out colour label fits or tunes a training model.

Species with fewer than four photos in either synthetic colour contribute zero instead of disappearing from the denominator. Constant predicted partitions also carry zero information. Geographic and environmental scores are compared with separately calibrated thresholds: the maximum 97.5th-percentile threshold across all four nuisance worlds at that signal amplitude. All four outcomes—geography only, environment only, both, neither—are retained.

## 2. Recovery and failure are both material results

| Planted scenario / recorded outcome | Amplitude 0.5 | Amplitude 1.0 | Amplitude 2.0 |
|---|---:|---:|---:|
| Geographic sharing / geography only | 45.4% | 100.0% | 100.0% |
| Environmental sharing / environment only | 13.2% | 94.4% | 100.0% |
| Confounded sharing / both | 19.8% | 98.6% | 100.0% |

Amplitude is the artificial latent signal relative to noise SD = 1. It has no established conversion to the empirical RGFCA rho, beta, or effect size.

The weak environmental signal was missed in 86.8% of evaluation worlds. At amplitude 0.5, environmental recovery was 66/500 with interval 10.51–16.45%. Thus the prototype cannot treat non-detection as absence. At amplitude 2, 500/500 detections have a Wilson lower bound of 99.24%, not proof of literally certain detection.

The highest observed nuisance false-sharing rate was **21/500 = 4.2%**, in the confounded mapping with weak independent environmental structure. Its interval was **2.76–6.34%**. This does not certify type-I error below 5%: intervals omit uncertainty in the estimated calibration thresholds, the nuisance family is finite, and real-data nuisance strength is unknown. All nuisance and four-way outcomes are retained in the full decision table.

## 3. Why overlap alone is insufficient

For amplitude 1 in environmental coordinates:

| Statistic | Shared environmental response | Independent species responses |
|---|---:|---:|
| Mean within-species separation, 1 minus Gaussian affinity | 0.248787 | 0.247905 |
| Mean species-disjoint transfer ARI | 0.107159 | 0.059873 |

The separation scores are almost the same, despite very different sharing assumptions. Species-specific responses can strongly separate their own colour groups while using different environmental directions. Transfer adds a different observable: whether a partition estimated for one species predicts a partition for another.

Even the raw transfer score is positive under independent species structure, so a zero benchmark is inadequate. Calibration against structured nuisance worlds is essential.

The affinity is the analytic Bhattacharyya coefficient of Gaussian densities. It is **not** the geometric intersection fraction of 95% ellipsoids, nor a biological niche-overlap probability. Total Gaussian ellipsoid volume, computed without colour, was identical under label replacement on paired coordinates: maximum difference **0.0**. This is an expected invariance, not a discovered biological relationship.

## 4. Identifiability ceiling

If environmental coordinates satisfy E = f(G), then a common environmental response h(E) is also the geographic function h(f(G)). No observational representation comparison alone can generally distinguish the two as exclusive causes.

In the confounded moderate-signal control, 98.6% of worlds were positive in both representations. We retain that ambiguity instead of declaring a winner. At weak signal, however, this control was often undetected or positive in only one representation. Detecting confounding is therefore not guaranteed either.

The harmonic example tests deliberately representation-dependent recoverability with a low-order Gaussian model. A more flexible geographic model could learn additional environmental preimages. It is not evidence that geographic commonness is genuinely absent whenever only the environmental Gaussian succeeds. The current holdout separates species, **not geographic regions**. All species in the synthetic benchmark draw from uniform spherical support.

Consequently this result does not establish causal selection, adaptation, shared real-world environmental thresholds, a measured fraction of sharing species, actual-RGFCA power, or temporal hypervolume turnover.

## 5. Actual photo geometry: audited, not substituted for the toy

A checksum-verified export from the frozen measured table contains only photo ID, species, coordinates, and the pre-existing classifiability mask. It contains no biological colour values.

| Stage | Species | Photo rows |
|---|---:|---:|
| Measured coordinates | 500 | 50,000 |
| Classifiable, before the species eligibility rule | 500 possible species | 25,377 |
| At least 40 classifiable photos per species | 369 | 21,424 |

The species eligibility loss is 131. Among the same 369 eligible species, the median retained photo count is 56 versus 100 in their complete coordinate pools. Complete-versus-retained geometric volume comparisons change sample size as well as locations and cannot identify a pure missingness effect. Geographic covariance effective rank describes the spherical Cartesian embedding, not ecological niche dimensionality.

The synthetic benchmark has **not** yet been calibrated on this actual geometry or on real climate. That is the next substantive extension: use the fixed photo support and outcome-blind environmental representation, preserve species-specific structure in the null, and test transfer across held-out species and geographic regions before opening a new observed-colour claim. Partial sharing and time-varying environments also remain outside this benchmark.

## Reproduction and verification

The reference run used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3 and scikit-learn 1.8.0. Fifteen unit tests passed locally, including a slow independent SciPy Gaussian-log-density calculation with sklearn ARI agreeing with the vectorized transfer statistic to 12 decimal places.

One test-only numerical correction preceded simulation: the even-coordinate identity used a tolerance of 1e-14 rather than bitwise equality after machine-rounding discrepancies up to 2.22e-16. No generator, scientific threshold, scenario or seed changed.

From the repository root:

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pandas==2.2.3 scikit-learn==1.8.0
python -m unittest discover -s tests -p test_hypervolume_transfer_benchmark.py -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/analysis/run_hypervolume_transfer_suite.py --output-dir /tmp/hypervolume-transfer-reproduction
```

Use a fresh output directory for a fresh reproduction. Existing arm outputs are retained by the suite and are not overwritten by the arm runner. Every arm contains its stage, mapping, world, amplitude and replicate identifiers. Summary generation checks complete calibration/evaluation counts and colour-invariant total volume.

### Source files

- [Frozen contract](supporting/hypervolume_transfer_benchmark_contract_v1.json)
- [Recorded result and uncertainty](supporting/hypervolume_transfer_benchmark_result_v1.json)
- [All 33 four-way decision distributions](../data/derived/hypervolume_transfer_decisions_v1.csv)
- [Actual coordinate/mask audit](supporting/hypervolume_actual_geometry_audit_v1.json)
- [Benchmark runner](../scripts/analysis/run_hypervolume_transfer_benchmark.py)
- [Unit tests](../tests/test_hypervolume_transfer_benchmark.py)
