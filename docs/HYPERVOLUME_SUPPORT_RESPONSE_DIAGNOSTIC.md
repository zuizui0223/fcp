# Hypervolume transfer: support overlap and response estimation

Date: 2026-09-07 (Japan). **Completed synthetic diagnostic; observed-colour inference remains closed.** Specification: `ee3f603f75f29e9fe60cce59d9d4d47552f8365e`. Parent snapshot: `6c12203d8c635fd62c64d6da3c28f36f1cf7b8dc`.

## Main result

There is usable signal in some synthetic worlds on the actual photo geometry, but neither overlap weighting, colour-prior correction, nor a simple conditional logistic model makes moderate shared responses reliably detectable. A diagnostic using the simulator's true source-species response detects moderate full sharing much more often than the fitted models. This implicates response estimation/model approximation as an important limitation; it does not identify sample size alone as the cause.

## Fixed experiment

We retained the existing 369-species frame, fixed training/evaluation species assignment, all four longitude-sector holdouts, 500-km buffer, 20 training and 20 evaluation species per fold, and 20 photos per species. Only coordinate support and the existing measurement mask were read. The environment remains the parent's artificial two-dimensional fourth-harmonic field, not temperature, dryness, CHELSA or other real climate.

A fresh seed generated 2,500 calibration and 3,500 evaluation worlds, 250 per arm, for 6,000 unique worlds. The same world and sampling schedule feed all four models and both weighting schemes. Four fold scores are averaged equally. Neither folds, source-target pairs, nor the eight alternative designs multiply biological sample size.

The four predictors are: the original equal-prior Gaussian/QDA, the same fitted class densities with empirical training colour priors, a fixed-penalty conditional logistic regression, and an oracle using the simulator's true response for each SOURCE species. The oracle never uses a target-species response parameter, but it still has unrealistic privileged knowledge and is not a deployable model or a formal power upper bound.

The overlap weight is the fraction of target environmental points inside a colour-blind 95% Gaussian ellipsoid fitted to a source's 20 points. Source weights are normalized separately for each target; every target retains equal outer weight. Zero-overlap targets and invalid colour-group predictions remain in the denominator. This is a Gaussian support proxy, not exact hypervolume intersection; scores still evaluate full target samples rather than restricting every evaluated photo to common support.

## Detection rates when all species share the response

All entries use the fixed UNKNOWN-strength calibration: for each model/weight, take the maximum empirical 97.5th percentile over ten nuisance arms, including independent species responses of strengths 0.5, 1 and 2. No favourable amplitude-specific threshold is substituted.

| Model | Moderate signal, equal weights | Moderate signal, overlap weights | Strong signal, equal weights | Strong signal, overlap weights |
|---|---:|---:|---:|---:|
| Original QDA | 1/250 (0.4%) | 0/250 (0%) | 71/250 (28.4%) | 122/250 (48.8%) |
| QDA with empirical colour priors | 1/250 (0.4%) | 0/250 (0%) | 79/250 (31.6%) | 125/250 (50.0%) |
| Conditional logistic model | 0/250 (0%) | 0/250 (0%) | 115/250 (46.0%) | 124/250 (49.6%) |
| Oracle source response | 171/250 (68.4%) | 56/250 (22.4%) | 232/250 (92.8%) | 233/250 (93.2%) |

Moderate and strong mean artificial amplitudes 1 and 2, not measured biological effect sizes. The oracle's 171/250 interval is 62.40–73.85%. Its success shows that the synthetic coordinates and labels are not uniformly devoid of exploitable signal, but the oracle is not an attainable estimator.

For original QDA at strong full sharing, overlap weighting adds 51 detections net: 61 worlds gain detection and 10 lose it. The paired difference is +20.4 percentage points (Monte Carlo standard error 3.12 points). This is conditional on separately estimated diagnostic thresholds. It is not a biological effect.

Overlap weighting is not universally beneficial: it lowers the oracle's moderate full-sharing detection from 68.4% to 22.4%. Scores, their nuisance distributions and the source weighting all change together. We do not select the better weighting separately after seeing each scenario.

At moderate signal with only 25% or 50% of species sharing, all three fitted models detected zero of 250 worlds for both weighting schemes. The oracle detects 14/250 at 50% sharing with equal weights and 3/250 with overlap weights; it also detects zero at 25%. Weak partial sharing remains particularly poorly identified.

The largest observed primary nuisance rejection rate across all eight designs and ten nuisance arms is 9/250 = 3.6%, with Wilson interval 1.91–6.70%. This does not certify error below 5%: calibration-threshold uncertainty, a finite nuisance family and joint comparisons across designs remain limitations.

## What the upstream geometry audit adds

| Fixed fold | Mean target-photo fraction inside a source ellipsoid | Weight-based effective source count out of 20 |
|---|---:|---:|
| 0 | 17.35% | 4.09 |
| 1 | 22.13% | 5.60 |
| 2 | 34.93% | 8.31 |
| 3 | 15.18% | 4.07 |

These are averages over the 500 fresh colour-blind schedules, not estimated real ecological niche overlap. Effective source count is (sum w)^2/sum(w^2), not a count of independent species. Most target samples have at least one overlapping source, yet their weight can be concentrated on few sources.

Thus a fixed species/photo eligibility rule can be well-defined while leaving limited information for response transfer. This does not authorize outcome-driven revision of the existing species frame. Future designs should distinguish sample-frame definition from within-species response learnability and cross-species environmental comparability.

## Strength conditioning is diagnostic, not a shortcut

Using only strength-1/noise nuisance arms for original equal-weight QDA yields 34/250 detections (13.6%) for moderate full sharing. But applying that cutoff to independent strength-2 responses yields 140/250 false sharing calls (56%). This oracle-known-strength calibration is retained as a secondary diagnostic only. Actual response strength is unknown, so this is not an acceptable replacement for the primary pooled-strength cutoff.

## Interpretation of the model comparison

For an exactly specified generative classifier, class-density log odds plus log(n1/n0) recover the posterior log odds corresponding to the empirical training prior. Holding priors equal can move a decision boundary when the sampling distributions differ. Nevertheless, Gaussian misspecification and environmental range differences remain; the prior correction tested here barely improves moderate-signal performance. Changing to a simple conditional logit also does not solve moderate-signal recovery. These are methodological comparisons, not evidence for real selection or an exclusive environmental rather than geographic cause.

## Verification and reproduction

The simulation finished once with exit code 0. Forty-seven local tests passed with the actual geometry enabled. Independent verification checked all 2,000 scheduled folds, 48,000 score rows, 80 calibration quantiles, 24 thresholds and 336 detection-rate/interval rows. There were 432 independent Gaussian/logistic prediction-and-ARI spot checks; maximum saved-fold-score discrepancy was 1.11e-16. All fitted logistic gradients met the convergence rule.

A verification-only JSON export initially failed on a NumPy boolean. Converting check flags to native booleans fixed serialization; no simulation output, seed, model, threshold or scenario changed. All original empirical decisions remain unchanged.

Runner: `scripts/analysis/run_hypervolume_support_response_diagnostic.py`.
Tests: `tests/test_hypervolume_support_response_diagnostic.py`.
Contract and result: `docs/supporting/hypervolume_support_response_diagnostic_{contract,result}_v1.json`.

With the checksum-verified coordinate-only export and its audit.json in the same input directory:

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pandas==2.2.3 scikit-learn==1.8.0
python -m unittest discover -s tests -p test_hypervolume_support_response_diagnostic.py -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/analysis/run_hypervolume_support_response_diagnostic.py --geometry /path/to/geometry_only.csv --output /tmp/hypervolume-support-response-fresh
```

Use a fresh output directory; existing outputs are not overwritten. The portable package also contains the input export, complete scores, schedules, all rates and intervals, the independent verifier and test logs. The repository validation workflow checks code/tests and recorded guards, not a new full simulation run.
