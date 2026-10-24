# FCP: source-authenticated shared thermal axis predicts the photographic white–chromatic direction

**Original source computation PASSED, 24 October 2026.** GitHub Actions https://github.com/zuizui0223/fcp/actions/runs/40066224301. This is a post hoc analysis of photographic floral-region hue mixtures, **not** validated genetic flower-colour morphs, pigment costs, experimentally supported climate adaptation or calibrated botanical reflectance.

## Main biological interpretation, conditional on photographic scope

The preceding #268 source study found that holding all seven climate variables separately fixed greatly reduced the amount of unique BIO5 climate information. This was a parameterization problem for *climate GROUP-level questions*: high correlation of BIO1 with BIO5 need not invalidate a biologically meaningful shared broad warmth gradient. We therefore constructed exactly five fixed outcome-blind macroclimate axes from the same seven WorldClim normals and adjusted the five axes mutually:

- **Shared thermal level** = (zBIO1 annual temperature + zBIO5 warmest-month maximum temperature) / 2.
- **Shared precipitation amount** = (zBIO12 annual precipitation + zBIO14 driest-month precipitation) / 2.
- Temperature seasonality BIO4, precipitation seasonality BIO15, and mean solar radiation SRAD as separately controlled axes.

Each source photo hue/climate model additionally controls broad 3D XYZ geography and the original observation calendar month and retains ≥5% unique variance. The mean signed photographed white-axis response to **shared thermal level** remains positive after controlling precipitation amount, climate seasonality and radiation. The precipitation-group response is smaller and lacks a reliable shared white direction. This is evidence for differing SOURCE PHOTO colour–climate **geometry**, not proof of two independent physiological selection pathways.

## Source opportunity and authentic data provenance

Original genuine photographed dataset: **447,500** independent photograph IDs across **4,475** named species (100 original photographs per species). **227,107** retained nine-colour classified source photo palette fractions; all 220,393 original unclassified photographs remained unassigned to hue/white phenotypes.

The source photo archive checksum SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6. Archived original WorldClim 2.1 BIO and monthly solar checksum SHA256 values: 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf and c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.

All **4,475 species** remain in the per-taxon output, including those failing the ≥40 classifiable photos, ≥20 original observers, ≥5 occupied source climate raster cells, original photo-month availability, or ≥5% original macro-axis unique climate variance eligibility filters. Complete original source photo and deterministic one-photo-per-observer separate panels are included.

## Result A: the common temperature-domain photo-white colour association persists

| Original photo panel | Thermal-domain n identifiable | Species-equal white-photo component partial-slope-like coefficient | Squared white-axis alignment of cross-species mean original nine-hue vector | White projection genus *five-domain simultaneous* 95% CI |
|---|---:|---:|---:|---:|
| Original full-source photos | **2,837** | **+0.008209** | **85.4%** | **[+0.004237,+0.013177]** |
| One photo per observer | **2,771** | **+0.008002** | **87.1%** | **[+0.004085,+0.012890]** |

The signed photo-white coefficient is an observed cross-taxon mean **partial covariance divided by remaining within-species climate-macro-axis variance**. The last column instead quantifies the pre-fixed UNIT q-projection and is in a different scaled coordinate. Do not call the result a causal per-1°C effect or a frequency of genetically white flowers.

Only **54.2%** (full) and **54.0%** (observer-thinned) of eligible species have a positive photo-white partial association, so many named taxa differ in direction. The mean shared thermal-white slope is positive in **all five pre-existing species-disjoint original photographic sampling cohorts** (discovery, validation, prospective, wave4 and new_U100), albeit each uses the same iNaturalist source and colour measurement assay.

| Source original cohort | Original thermal-eligible taxa | Mean white-photo partial slope |
|---|---:|---:|
| Discovery | 331 | +0.012570 |
| Validation | 323 | +0.012985 |
| Prospective | 322 | +0.007827 |
| Fourth photo cohort | 297 | +0.003559 |
| Fifth photo cohort | 1,564 | +0.007261 |

These are *descriptive source-cohort means*, not a claim that every cohort individually passed a preregistered significance test or that photo sources were independently calibrated.

## Result B: matched-species temperature and precipitation photo-hue directions differ

The source-correct comparison uses exactly the same taxa in both macroclimate axes, avoiding the misleading comparison of photo-hue orientations measured on different taxon sets.

| Same-source eligible original photos | Matched taxa | Thermal group white-axis alignment | Precipitation group white-axis alignment | Paired genus-bootstrap alignment difference, descriptive 95% |
|---|---:|---:|---:|---:|
| All eligible source photos | **2,689** | 83.7% | 24.3% | 59.4 percentage points [6.0, 86.8] |
| One photo per observer | **2,625** | 84.8% | 12.7% | 72.0 points [10.0, 88.3] |

These are the squared alignment of **mean climate-associated nine-hue photo palette covariance vectors**, not of individual source photo phenotypic variance or any molecular/heritability explanation. The precipitation-group average photo-hue vector is relatively weak, so its orientation is uncertain even after matching taxa. The difference was examined post hoc and CI is exploratory.

Within the full eligible original photo set, the precipitation-group mean photo-white partial-slope coefficient is +0.002230 (3,022 taxa) with a five-domain q-projection CI that crosses zero [−0.002041,+0.006771]; one per observer gives +0.001779 (2,957 taxa), CI [−0.002095,+0.005868]. Chromatic-only yellow/bronze source-photo contrasts are not established after five-domain source familywise intervals. Do not infer specific causal water-linked pigment allocation.

## Result C: the warmth–white photo association has both within- and between-photographic-label components

For the **exact same 1,376 original named taxa** containing ≥5 authentic photos in two original coarse PHOTO-colour categories, the shared warmth white-photo source partial slope decomposes:

- Total +0.012428
- **Between-original-photo-label** +0.008431
- **Within-original-photo-label** +0.003997
- Between contribution to the **signed mean partial slope**, 67.8% (whole-genus bootstrap marginal 95% [44.2%,85.4%]).

After one source photo per observer, exactly 1,299 taxa remain:

- Total +0.012187, between original photographic categories +0.008050, within categories +0.004137;
- Between algebraic fraction **66.1%** (whole-genus bootstrap marginal 95% [39.2%,84.3%]).

**This 66–68% is NOT the percentage of genetic morph frequency shifts, heritable color polymorphism, or selection/cost explained.** The categories are camera/photo ROI-derived and may mix petal colour, nectar-guide-like patterns, flower stage, fruit/sepals, cultivation status and measurement conditions. Both components are descriptive photo-state covariance terms.

## Limitations and value for the manuscript

- The strongest supported result is a **repeated, small PHOTO-colour geometric pattern** shared by many taxa, with the shared warmth domain pointing disproportionately along the photographed white-versus-other axis. This does not establish whether individual species' *real petals* are genetically polymorphic.
- The photo provider/measurement assay is common to all five source cohorts. Approximately half the source photographs are not colour-classifiable, which may introduce biased taxa/site representation and source camera/observer seasonality.
- The earlier independent **flower-minus-background BIO5** source validation failed to support a flower-specific climate signal (rho −0.00328, p=0.82). This new macroclimate model has **not controlled a matched background for all source photos**, so it cannot overturn that caution.
- Original WorldClim long-term BIO and solar normals are geographic habitat climate context, not temperature or moisture observed during each floral development. We remove only broad XYZ linear geography and calendar month and the other four fixed macro axes; unmeasured nonlinearity, biological history and cultivation can remain.
- Critically, grouped warmth is **not uniquely BIO5**, and grouped precipitation is **not uniquely BIO12**; the axes intentionally pool correlated original climate variables. Grouping is a separate scientifically meaningful question but not a way to claim otherwise nonidentifiable individual-climate effects.
- No measured pigment production cost, actual per-morph fitness, population gene flow, molecular pathway or pollinator preference can be derived from this photographed macroecological resource.

## Reproducibility

- Main source script: scripts/analysis/fcp4475_shared_thermal_moisture_photo_hues_20261024.py.
- Source dataset, exact climate archive SHA, analytic original-photo-month null mathematical fixture, 5 macro-domain controls, 999 whole-genus bootstrap and per-species output reproducibility: https://github.com/zuizui0223/fcp/actions/runs/40066224301 (**PASS**).
- Source artifact: fcp4475-shared-thermal-moisture-source-photo-hues-20261024 (artifact ID **12549840336**, ZIP sha256 **3f290b8028bbfcdb6fa2d1a4ae85f3155b21786d7c8d48175af4805c78c099a0**). Includes complete original and observer-thinned 4,475-taxon support CSVs, numerical original 4475×5×3×9 photo-hue coefficients/partial covariances and original ecological result JSON.
- Permanent machine-readable source result: results/fcp4475_shared_thermal_moisture_photo_hues_verified_20261024.json.
- No source photo IDs, original camera/observer IDs or individual raw coordinates are redistributed in these derived result tables.

## Scientifically bounded manuscript text

> Within the sampled photographic nine-colour distributions of 4,475 flowering-plant taxa, the shared temperature level represented by annual mean and warmest-month maximum temperature was associated with a small cross-species shift toward white-rich photo palettes, even after adjusting for precipitation amount, temperature and precipitation seasonality, mean solar radiation, broad geography and observation month. The common mean thermal colour vector predominantly followed the white–chromatic axis, whereas the precipitation-domain photo-colour direction was weaker and comparatively less aligned with that axis. The positive thermal association recurred across five species-disjoint photography cohorts but was not present in every taxon. These results support repeated photographic colour–climate geometry, not causal climate adaptation, heritable white flower morph frequencies or evolutionary pigment-maintenance cost.
