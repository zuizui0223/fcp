# Actual-support transfer: the current estimator misses planted sharing

Date: 2026-09-07 (Japan). Status: **simulation completed and independently checked; no empirical transfer claim opened**.

## Decision first

Do not carry the original uniform-sphere success rate into the real-photo analysis. On the actual coordinates, with the newly frozen nuisance calibration and sampling design, the retained-photo, species-and-region-disjoint estimator detected a planted, fully shared environmental response in **0/250 worlds**. It also detected the fully shared geographic response in **0/250**. These are failures of recovery under the tested methodological conditions, **not evidence that real flowers lack common responses**.

The specification was committed before synthetic outcomes at `db9e34175a22c8f21bde8f018db69b26dde7cc58`. The implementation commit is `0bae260b559a42a9d2506bb7512016094c5916b3`. No existing G1, H1-H6, environmental mean, interaction or species-disjoint biological support decision was changed. No real biological colour values or real environmental rasters were read.

## 1. What was executed

The new run contains **6,500 distinct synthetic worlds**: 2,500 calibration and 4,000 evaluation worlds. Each world is evaluated under four paired designs and four geographic folds, giving 104,000 fold-world-design evaluations. Those evaluations are not 104,000 independent worlds, and neither folds nor species pairs increase the biological sample size.

Only the checksum-verified five-column coordinate export was used: photo identity, species, latitude, longitude and the pre-existing classifiability mask. The source frame contains 369 eligible species, 36,900 complete-coordinate rows and 21,424 classifiable rows. The mask is measurement-derived; it is not assumed to be random or independent of real colour. The same 369 species define both coordinate-pool comparisons; the 131 excluded species are not restored.

A fixed SHA256 species split allocates 184 species exclusively to training and 185 exclusively to evaluation. Each synthetic world samples 20 eligible training species and 20 different evaluation species per fold, with 20 photos per species. The target photos lie inside one of four fixed 90-degree longitude sectors. In the blocked design, training photos lie outside that sector and at least 500 km from every complete-coordinate reference photo in it. This is a prespecified stress-test distance, not an estimated ecological independence range.

| Fold longitude interval | Eligible training species | Eligible evaluation species | Smallest training-to-sector-photo distance |
|---|---:|---:|---:|
| -180 to -90 degrees | 133 | 60 | 500.31 km |
| -90 to 0 degrees | 135 | 62 | 500.79 km |
| 0 to 90 degrees | 141 | 59 | 502.30 km |
| 90 to 180 degrees | 165 | 26 | 509.73 km |

All four folds passed the frozen capacity rule. Their equal-weight mean is the test statistic; no favourable fold was selected.

The four designs cross retained versus complete-coordinate pools with species-only versus species-and-region separation. The species-only comparator uses the same species identities and the same target-sector photos, but permits training photos from all sectors. This is not a within-species fit or a random-photo holdout. Random photo priorities, sampled species and the synthetic label assigned to each original photo are paired across designs. The complete-coordinate comparison therefore holds draw size at 20, rather than comparing 56 versus 100 observations directly.

## 2. The full recovery result

The positive worlds plant an amplitude-1 response shared by 25%, 50% or 100% of the 369 species; nonsharing species retain independent same-domain response directions. All nuisance worlds and all four possible classifications are retained in the accompanying decision table.

For **100% sharing**:

| Coordinate pool / training restriction | Geographic world: geographic detection | Environmental world: environmental detection |
|---|---:|---:|
| Retained / species only | 0/250 | 0/250 |
| Retained / species and region | 0/250 | 0/250 |
| Complete, same species / species only | 0/250 | 1/250 |
| Complete, same species / species and region | 0/250 | 1/250 |

These columns count positivity in the planted representation. Wrong-representation classifications are not counted as recovery. For example, the complete/blocked environmental world produced three geography-only and one environment-only classifications, rather than four correct environmental recoveries.

For zero detections out of 250, the Wilson 95% interval is approximately **0-1.51%**, conditional on the estimated fixed thresholds and this simulated-world family. The primary retained/blocked partial-sharing environmental arms detected the environmental representation in 0/250 at 25% sharing and 1/250 at 50% sharing. This small nonmonotonicity is Monte Carlo variation, not a claim that half-sharing is easier than full sharing. The geographic partial-sharing arms were 0/250.

The largest recorded nuisance false-sharing rate across the four designs was **11/250 = 4.4%**, interval **2.47-7.71%**. It occurred for strong independent environmental responses in the retained/species-only design. This does not certify a 5% upper error bound; threshold estimation uncertainty and nuisance-family generality remain outside these intervals.

## 3. What the diagnostics do and do not identify

### Restoring the coordinate pool did not recover the method

For the fully shared environmental world under blocked evaluation, the paired mean difference in transfer score, complete minus retained, was **-0.000298**, with Monte Carlo standard error 0.000225 and normal Monte Carlo interval **[-0.000739, 0.000142]**. No clear improvement was recovered by restoring the complete coordinate pool while holding species and draw size fixed.

This does not establish that measurement missingness is harmless: the 131 excluded species remain absent, the retained mask can be colour-dependent, and the artificial labels are not unmeasured real colours. It addresses a conditional support/sampling contrast only.

### Spatial exclusion reduces transfer, but is not the whole failure

For the same retained fully shared environmental worlds, the species-only minus blocked mean score was **0.001185**, with Monte Carlo interval **[0.000955, 0.001415]**. Nevertheless the species-only design also failed to recover environmental sharing under its fixed threshold. Spatial exclusion alone therefore does not explain this benchmark's failure.

Only **25.20725%** of the potential training/evaluation species pairs contributed an informative environmental prediction in the primary fully shared world. This requires both species to have at least four synthetic photos per colour and the predicted target partition to be nonconstant. The other pairs remain in the denominator at zero; they were not silently discarded. This diagnoses a bottleneck in usable partition predictions, not a separately proven cause of all non-detection.

### Strong private structure can outscore weaker shared structure

Primary retained/blocked environmental mean transfer scores were:

| Planted world | Signal amplitude | Mean transfer ARI |
|---|---:|---:|
| Independent species responses | 1 | 0.009450 |
| Fully shared response | 1 | 0.011721 |
| Independent species responses | 2 | 0.015873 |

The score is not ordered by sharing fraction independently of signal strength. The fixed environmental threshold was **0.026227**: the maximum 97.5th percentile over ten nuisance calibration arms, including strong independent responses. The amplitude-1 independent-environment arm alone had a lower calibration quantile, about 0.016221, but it was not substituted after outcomes. No alternative-threshold discovery is reported.

This is more restrictive than the original benchmark, whose thresholds were calibrated separately at each known artificial amplitude. The change was frozen deliberately because true nuisance strength would be unknown in an empirical application.

## 4. How this changes the earlier 94.4% result

The previous 94.4% result remains a result of its original uniform-sphere, moderate-signal construction. The current benchmark differs in actual sampling support, numbers of sampled species, restricted target regions and calibration across nuisance strengths. It is **not a one-factor experiment proving that coordinates alone caused a drop from 94.4% to zero**.

The scientifically useful conclusion is narrower: success under the original constructed geometry did not qualify the method for empirical application. Under the present prespecified actual-support stress test, fully shared signals were not reliably recovered.

The environment is still the parent's synthetic fourth-harmonic transformation of spherical coordinates, not temperature, VPD, CHELSA or SoilGrids. ARI identifies a partition up to colour-label reversal: it does not require all species to share the same directional colour response. No causal geographic-versus-environmental conclusion follows.

## 5. Verification and next gate

Thirty-four tests passed locally (15 parent tests and 19 new tests). These include species disjointness, all photo identities, mask handling, geographic exclusion, spherical distances, class-denominator rules and independent SciPy/sklearn score calculations. All 52,000 saved score rows, all 80 calibration quantiles and all 448 decision-count rows were checked. Across 480 actual-support direct-Gaussian spot checks, the maximum transfer discrepancy was **9.97e-17**; no simulation was rerun to obtain a more favourable result.

The new scientific gate stays closed: do not open an observed flower-colour transfer claim with this estimator. The next methodological target is transfer within overlapping environmental support and calibration conditional on species-specific signal strength. It requires a separate frozen design; it must not reclassify the present failure. Physical climate, empirical transfer, causal identification and temporal hypervolumes remain untested here.

## Reproduction

The portable package includes the original coordinate-only input and audit, both runners, all test and verification scripts, frozen protocol, sampling schedules, all per-world/per-fold statistics and checksums. It does not include biological colour values.

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pandas==2.2.3 scikit-learn==1.8.0
HYPERVOLUME_GEOMETRY=inputs/geometry_only.csv OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m unittest discover -s tests -p 'test_hypervolume*.py' -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/analysis/run_hypervolume_actual_support_blocked.py \
  --geometry inputs/geometry_only.csv --output-dir reproduction_new
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/analysis/verify_hypervolume_actual_support_blocked.py \
  --geometry inputs/geometry_only.csv --results reproduction_new \
  --output reports/reproduction_verification.json
```

Use a fresh reproduction directory; the runner refuses to overwrite completed arm outputs. The package's saved reference run is under `results/`.

### Files and methodological references

- Frozen specification: `supporting/hypervolume_actual_support_blocked_contract_v1.json`.
- Checked compact summary: `supporting/hypervolume_actual_support_blocked_summary_v1.json`.
- All 64 four-way decision distributions: `../data/derived/hypervolume_actual_support_decisions_v1.csv`.
- Valavi et al. (2019), *blockCV*, Methods in Ecology and Evolution, DOI 10.1111/2041-210X.13107, motivates evaluating structured data with separated folds. The present fixed-sector protocol is our own stress test, not an implementation or validation of every blockCV design.
- scikit-learn `adjusted_rand_score` documentation describes chance adjustment and label-permutation invariance. The present constant-partition-zero and minimum-colour-count rules are additional explicit information guards.
