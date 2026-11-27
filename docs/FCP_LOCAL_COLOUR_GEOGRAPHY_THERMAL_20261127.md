# Local colour contrasts, geographic structure and thermal associations (2026-11-27)

## Scope and data

Exploratory secondary analysis of the existing source-aligned Table S10, at manuscript PR #274 head cbe5301d4a3c57d91f001655a455ecefa0acb550. Source Git blob 887991a51009461a640f4cd2ab4b13583f53516b. Retrieved directly from GitHub, decompressed with verified gzip CRC32 715323063 and uncompressed size 291125 bytes. 4475 distinct taxon IDs, 2400 with repeated photographic colour categories. No new photo acquisitions, image classification, or fitted thermal models. The direct numerical calculation used JavaScript with 3999 whole-genus bootstrap resamples (LCG seed 20261127); the included equivalent Python standard-library script is provided for reproduction but has not yet run in CI.

## Question

Do local contrasts in photographic colour necessarily preclude broad geographic structure or a mean thermal-white association? For all taxa with repeated colour categories, select those with an opportunity to detect a contrast at <=5 km, same calendar month and year, distinct observers. A local contrast is a detected pair with different colour categories; nondetection is not proof of local monomorphism. Restrict the principal joint description to taxa with both the 50 km geographic estimate and the shared thermal-white estimate. Retain all 415 local-opportunity taxa in the descriptive counts. Positive geographic contrast = range-wide minus local mismatch; positive thermal slope = increasing photographed white fraction on the existing standardized thermal index.

## Direct results

| Panel | Local contrast detected | Not detected |
|---|---:|---:|
| All repeated-colour taxa with local opportunity | 228 | 187 |
| Mean geography, all local-opportunity taxa (percentage points) | +1.00677 | +2.91644 |
| Common geography and thermal estimates | 108 | 121 |
| Common-panel genera | 91 | 104 |
| Positive geography, common panel | 59 | 66 |
| Positive thermal, common panel | 69 | 58 |
| Both positive, common panel | 39 | 28 |

Common panel: 229 taxa, 176 genera. Joint whole-genus bootstrap draws (same sampled clusters for both strata); 3999 valid replicates for each quantity. All results below in original proportion units, not percent.

| Statistic | Estimate | Marginal 95% interval |
|---|---:|---:|
| Geographic, detected | 0.011308695 | [-0.000394822,0.023314212] |
| Geographic, not detected | 0.022828605 | [0.007731330,0.038577193] |
| Geographic, detected minus not detected | -0.011519910 | [-0.031013899,0.007538275] |
| Thermal, detected | 0.027237922 | [0.004760129,0.048992960] |
| Thermal, not detected | 0.005020325 | [-0.013334427,0.022736385] |
| Thermal, detected minus not detected | 0.022217597 | [-0.005887631,0.051410087] |
| Positive geography fraction, detected | 0.546296296 | [0.447901786,0.642857143] |
| Both positive fraction, detected | 0.361111111 | [0.267322465,0.455445544] |

## Interpretation and gate

Local photo-colour contrasts and a positive geographic point estimate can occur in the same taxon. Neither group contrast is supported by its unadjusted marginal interval; do not claim that local coexistence increases thermal sensitivity or reduces spatial sorting. The geographic estimator and local detection share source photos; the association may partly reflect mathematical overlap. Nondetection depends on the number of local pair opportunities and composition, which Table S10 does not contain. Positive point estimates are not species-specific statistical significance. The thermal estimate is based on historical location climate normals, not experimental temperature response, and the image pipeline does not establish homologous-petal pigmentation, heritability, fitness or population genetic differentiation. The groups differ in thermal estimability (108/228 vs 121/187); avoid extrapolating the common-panel contrast. Intervals are post hoc marginal, not multiplicity-controlled or full phylogenetic models.

Next biological test: on the original photo source, construct independent observer partitions for local detection, geographic sorting and thermal estimates, retain per-taxon local pair opportunity and composition, and evaluate whether local colour contrast persists alongside geographic sorting beyond a species-conditioned matched null. This is a new estimand and must not be represented as already tested. Keep the current manuscript and its frozen results unchanged.

Reproduction from repository root:
```sh
python scripts/analysis/fcp_local_colour_geography_thermal_20261127.py
```
