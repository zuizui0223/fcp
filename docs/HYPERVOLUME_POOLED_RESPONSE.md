# From colour hypervolumes to uncertainty-aware shared-response prediction

Date: 2026-09-07 JST. **Synthetic method benchmark complete; real-colour inference remains closed.** Specification was committed at `d4ba51f5cb5dabf9080a2d58358f49191a673cf5`, after the earlier method failures but before the new outcomes. Existing G1, thermal, dryness, interaction and other empirical decisions are unchanged.

## Answer first

The fixed actual-coordinate design can recover moderate full sharing without receiving true species response parameters. The predeclared linear-logit joint model detected **248/250 (99.2%; Wilson 95% 97.13–99.78%)**, while original QDA on the same fresh worlds detected **2/250 (0.8%)**. At 50% sharing the new joint model detected **93/250 (37.2%; 31.44–43.34%)**. At 25% sharing it detected only **8/250 (3.2%)**: sparse partial sharing remains poorly identified.

This is not a better estimator of hypervolume intersection or evidence for real climatic colour boundaries. It connects the environmental representation to a model of shared partitions, retains response uncertainty, and evaluates whether a learned distribution improves prediction for different species in different regions.

## Fixed design and methods

The source is the same coordinate-only export: 369 eligible species and 21,424 retained photos, with the pre-existing measurement mask. No biological colour values or real environmental rasters were read. The environment is the earlier two-dimensional artificial fourth-harmonic field.

For each fold: 20 training species, 20 different evaluation species, 20 photos per species, the same four longitude-sector holdouts and 500 km buffer. Every global synthetic world has one response and sign per species and one noise draw per original photo, reused across its folds and models. Four fold scores are averaged; folds and species pairs are not independent simulation replicates.

There were 2,500 calibration worlds and 3,500 evaluation worlds, 250 per arm, for 6,000 unique worlds. All five methods use the identical fresh worlds and schedules. Every method has one threshold, the maximum empirical 97.5th percentile across all ten independent-structure nuisance arms and strengths. Cutoffs are not matched to the unknown true signal strength after the fact.

Two response families are used. The saturating probit knows the generator's link shape and origin but not any true species axis, sign or strength. The linear logit deliberately differs from the generating link and is the predeclared practical-family candidate. Both retain a strong zero-intercept/common-origin assumption.

For each species, likelihood is integrated over 48 unoriented axes, a fixed strength grid and one whole-species colour reversal. This sign is not independently redrawn per photograph. The species likelihood ratio at an axis is its likelihood divided by its own uniform-axis marginal likelihood. The independent baseline therefore allows strong species-specific structure; it is not a flat-colour or label-shuffle baseline.

The joint model learns an axis-and-sharing-fraction posterior using training species only. Its sharing grid includes zero, 0.25, 0.5, 0.75 and one. The target score is the marginal predictive log-density improvement over independently distributed species axes, averaged per target and photograph. Evaluation labels never update the training posterior. The fitted sharing fraction is a modelling quantity, not a validated estimator of population prevalence.

The individual-source comparator retains the same species response uncertainty but averages individual source-axis posteriors without fitting the shared-fraction model. Original equal-prior QDA/ARI is retained as a separate comparator. Unlike QDA, the likelihood methods retain constant-colour samples; changes in likelihood, priors, score and information use prevent attribution of the entire improvement to pooling alone.

## All positive outcomes

Each entry is detections out of 250 evaluation worlds.

| Method | Moderate 25% sharing | Moderate 50% sharing | Moderate full sharing | Strong full sharing |
|---|---:|---:|---:|---:|
| Original QDA | 0 | 0 | 2 | 99 |
| Linear logit, individual-source uncertainty | 0 | 12 | 246 | 250 |
| Linear logit, joint training | 8 | 93 | 248 | 250 |
| Saturating probit, individual-source uncertainty | 0 | 12 | 248 | 250 |
| Saturating probit, joint training | 6 | 88 | 249 | 250 |

At full moderate sharing, retaining uncertainty and changing the predictive score already gives 246/250 for individual-source logit; joint training adds only two net detections. At half sharing, joint training adds 81 detections with none lost relative to that comparator: **+32.4 percentage points**, paired Monte Carlo SE 2.97 points. This contrast still combines a different aggregation with a shared-fraction model; it is not a one-factor causal attribution.

These percentages cannot be compared directly to the earlier oracle ARI as a power upper bound: the oracle used a different score and validity rule, and was never an optimal-power bound.

## Failure controls and limits

The highest nuisance rejection was **12/250 (4.8%; Wilson 95% 2.77–8.20%)** for strong independent geographic responses, in both logit designs. Do not claim guaranteed error below 5%. Intervals condition on estimated thresholds, exclude their calibration uncertainty, and do not provide a joint guarantee across five designs. The entire 70-row evaluation family is retained in the result record and portable CSV.

The result is conditional on a two-dimensional common environment, zero-centred boundaries, finite priors and an isotropic independent-axis nuisance family. It has not tested unknown or species-specific thresholds, broader link misspecification, nonuniform independent direction distributions, phylogenetic dependence, real climate, temporal hypervolumes or multicolour responses. An environment-only predictive result is not identification of an exclusive environmental cause.

In particular, 25% moderate sharing is still essentially not recoverable here, and 50% sharing is missed in 62.8% of worlds. Non-detection is not absence. The successful full-sharing arm does not open observed-colour inference by itself. The next methodological prerequisite is robustness to unknown/shifted response thresholds and richer independent-structure nulls, before a separately specified real-climate analysis.

## Verification and execution record

67 local tests passed with actual geometry; all 20 new tests also passed on GitHub. Independent verification checked all 2,000 scheduled folds, 30,000 score rows, 50 calibration quantiles, five thresholds and 70 counts/intervals. Separately implemented probability calculations agreed in 192 fold/model checks (maximum error 1.41e-16); eight complete QDA fold calculations using SciPy densities and sklearn ARI agreed within 9.37e-17.

The first monolithic process hit the execution time limit during calibration, before evaluation or exported scientific rows. The same hashed functions, seeds and scenarios were recovered in nonoverlapping deterministic shards. A distinct one-line source-upload transcription was restored on GitHub to the exact executed SHA256. The downloaded CI source was byte-compared with the local scientific runner and tests. Neither incident changed any scientific computation or decision.

Runner SHA256: `e87f8343fdd84f64c7245dff0399fb6b33ec8040f5290ffe5a8fdf1a1158581d`.

## Reproduction

Use the same checksum-verified `geometry_only.csv` and `audit.json` from the geometry-only export. The portable package contains these inputs, the parent dependencies, all score rows, schedules, rates, independent verifier, test logs and technical recovery record.

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pandas==2.2.3 scikit-learn==1.8.0
HYPERVOLUME_GEOMETRY=inputs/geometry_only.csv OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_hypervolume*.py' -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/analysis/run_hypervolume_pooled_response.py --geometry inputs/geometry_only.csv --output /tmp/fcp-pooled-response-fresh
```

The output directory must be new. The standard runner uses the identical functions as the recovered shards. Verification is not a favourable re-run of the simulation.

Sources: the frozen specification and `docs/supporting/hypervolume_pooled_response_result_v1.json`; the portable detailed `decision_rates.csv`, `paired_detections.csv` and verification records. Stable finite-mixture marginalization is implemented in log space, consistent with the Stan User's Guide chapter on finite mixtures. No MCMC, LOO approximation, or empirical-climate validation is claimed.
