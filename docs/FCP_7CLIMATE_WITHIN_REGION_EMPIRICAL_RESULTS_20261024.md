# FCP: within-region falsification of seven-climate continuous nine-hue association

**Original-data analysis and independent paired-species audit PASSED, 24 October 2026.** This is a *source-photography observational test*, not confirmation of stable heritable flower morphs, flower-specific climatic response, adaptive benefits or genetic population-level environmental sorting.

## Original ecological question and why this matters

The upstream [seven-climate full-photo analysis, PR #265](https://github.com/zuizui0223/fcp/pull/265) found a small positive signed covariance between the photographic white-palette fraction of flowers and long-term maximum temperature of the warmest month (BIO5) after a **global 3D linear spatial coordinate and month adjustment**. Among 1,578 taxa with at least two well-supported original 4-state photographic categories, +0.005278 signed mean photo-BIO5 covariance exceeded species-month-matched null; 75.8% of that mean signed covariance came from the original photo-category-between term, not a genetic morph frequency.

The scientifically consequential remaining threat was **regional geographic structure that is not adequately absorbed by global XYZ**. This new analysis tests climate vs the continuous nine-palette at **fixed 5° and 1° geographical blocks** and calendar month in the same actual source. It does not overwrite any upstream analysis.

## Authenticated source and model

Source: 447,500 genuine measured iNaturalist photographs / 4,475 species, 227,107 original four-state-classifiable photos with 9-component ROI palette, original photo CSV SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6. WorldClim 2.1 BIO1/BIO4/BIO5/BIO12/BIO14/BIO15 and twelve-month averaged SRAD archives were byte-verified.

Each species requires ≥40 classified, valid-coordinate photos, ≥20 source observers, ≥5 occupied WorldClim 10-arcminute cells, and photo date. Seven climate axes are standardized by the **within-species observed climate SD** then residualized on indicator terms for the photo's fixed geographical block (5° or 1°) and calendar month.

For each predictor in each taxon, ≥12 actual photo rows from at least two geographical-block × calendar-month strata must contain **more than one original climate raster cell and genuinely different values**. A predictor must retain ≥5% of the original standardized climatic variance after geographical+month fixed effects, or it is explicitly **not identifiable**, not an ecological negative result. The 1° block is not a 1-km locality and the blocks are not equal-area; original location error can reach 5 km.

Then 199 spatial-block × calendar-month permutations exchange original climate residual vectors only among photo rows in the SAME block and SAME calendar month, holding original photo category and nine-palette measurement intact. Species means carry equal weight. Refit colour-blind one-original-photo-per-observer. Test total signed white fraction and 9D hue norm separately; Holm-correct within each seven-climate exploratory family. Original four-state photo-label between/within components are exact photographic covariance decomposition only.

### Geographic opportunity from original source

| Original-photo test | Taxa eligible for any regional climate effect | Taxa with BIO5 contrasts | Taxa with identifiable locally mixed photo labels for BIO5 |
|---|---:|---:|---:|
| 5° block, all available photos | 3,057 | 2,925 | 1,238 |
| 5° block, one photo per observer | 2,984 | 2,846 | 1,152 |
| 1° block, all available photos | 1,540 | 1,377 | 497 |
| 1° block, one photo per observer | 1,449 | 1,303 | 454 |

The reduction from 5° to 1° partly reflects sampled region/month/climate-cell opportunity; the less informative source taxa remain in all four 4,475-row ledgers, not converted to zero temperature responses.

## Result A: signed BIO5/white-photograph association weakens under finer region matching

Original source regional permutation results:

| BIO5 photo-white signed covariance excess above own block/month null | Taxa | Effect | Holm across 7 climate axes |
|---|---:|---:|---:|
| 5° all photos | 2,925 | +0.00215250 | 0.035 |
| 5° one per observer | 2,846 | +0.00206132 | 0.035 |
| 1° all photos | 1,377 | +0.00096935 | 0.035 |
| **1° one per observer** | **1,303** | **+0.00069396** | **0.120** |

The latter is **not supported after the seven-variable correction**. This test cannot justify a universal within-locality tendency for warmer climates to produce more photographic white hue, much less white genetic flower morphs. The 1° sample and its climatic variance differ from the 5° sample.

Within the 1° observer-thinned sample, the **9D photo hue-vector norm** nevertheless shows a tiny excess +0.00105560 over matched null, Holm7 p=0.035. This is not a signed shift toward white; it may represent other hue directions, crop/camera variation or non-homologous floral organs. Its original-photo-label **within** component in 1°-thin has 9D norm excess +0.00050774, Holm7 p=0.035, but signed-white within effect +0.00010790, Holm7 p=1.000. The original-photo-label **between** component at 1°-thin was not supported in signed white (n=454, Holm7 p=0.680) or 9D norm (Holm7 p=0.825).

Thus an apparent positive unsigned hue-distance association does **not** rescue a robust common direction of white enrichment within small geographical regions.

## Result B: control for the changing species and residual climate variance

The 5° and 1° effect magnitudes are not directly comparable because narrowing the regional fixed effects mechanically reduces climatic residual variance and may drop taxa. Independently take the **exact species intersection** with eligible BIO5 at both scales: 1,377 original and 1,303 observer-thinned species.

For each species, calculate a null-corrected slope-like observational association:

    beta_proxy = (observed signed white photo/temperature covariance minus
                  region×month permutation mean)
                 / (region+month climate residual variance),

where each species' original climate was scaled by its observed SD before regional residualization. This is not a causal slope or temperature response per °C. Bootstrap **whole genera** 1,999 times with 7-climate-axis simultaneous critical-value intervals within each statistic family.

| Same-taxon BIO5 beta proxy | Original 1,377 taxa | Observer-thinned 1,303 taxa |
|---|---:|---:|
| 5° regional estimate | +0.006485 | +0.005570 |
| 5° seven-climate simultaneous 95% interval | [+0.003375,+0.009595] | [+0.002504,+0.008637] |
| 1° regional estimate | +0.003828 | +0.002367 |
| 1° seven-climate simultaneous 95% interval | [−0.000966,+0.008623] | [−0.002379,+0.007113] |
| **5° minus 1°** | +0.002656 | +0.003203 |
| 5°–1° seven-climate simultaneous 95% interval | [−0.001330,+0.006643] | [−0.000966,+0.007373] |

The *median retained climate variance fractions*, in exactly those two matched BIO5 species sets, fell from **0.565 to 0.178** in the original photos and **0.562 to 0.174** after observer thinning.

**Conclusion of the paired audit:** estimates are numerically larger at 5° than 1°, but the paired *difference* does **not** exclude zero under seven-climate simultaneous uncertainty. A drop in raw covariance does not itself prove the biological temperature slope weakened, and a positive fine-regional estimate with an interval spanning zero does not prove absence of any true biological effect.

For seven-climate simultaneous 1° slopes, **BIO5 does not exclude zero in either photo panel**. Other axes must not be cherry-picked from these exploratory families as a universal climatic adaptation law.

## Interpretation ceiling and previous contrary evidence

1. **Primary observational positive:** species' photographic colour contains a small 5° regional BIO5 signal above region/month-matched null; the fine 1° signed-white effect fails observer-thinned seven-predictor correction.
2. **Repeated weak photographic colour geometry:** fine 1° residual 9D hue norms remain small positive but are not same as more white pigment absence or genetic morph frequency. 
3. **No robust photo colour-morph inference:** the source photo palette can include senescent flowers, buds, fruit/sepals, petal guide blotches, camera/background variation and cultivated plants; 220,393 source photos were unclassifiable. Even a conditional positive is not causation by temperature or a per-anthesis exposure experiment.
4. **Earlier independent technical negative remains:** flower-minus-background validation on a separate species-disjoint cohort was unsupported (rho −0.00328, p=0.82); strict ≤5km same-year/month distinct-observer photographic null was also unsupported. This new test does not replace those limits.
5. **No genotype, fitness, cost-benefit or true natural local maintenance** have been measured. WorldClim averages are long-term climate normals, not temperature during an individual flower's anthesis.
6. Sampled 5° and 1° blocks are unequally sized in kilometres by latitude, and source opportunity (sampling month, number of locations, photographer geography, climate cells) is nonrandom. The 199 permutation worlds cannot support more precise numerical p-values than 0.005 and are exploratory across multiple model families.

## Reproducibility and receipts

- New original photo/raster computation: scripts/analysis/fcp4475_seven_climate_within_region_falsification_20261024.py
- Authenticated seven-climate full-source run: https://github.com/zuizui0223/fcp/actions/runs/40056835373 — PASS; original 4,475-row photo ledger × 5°/1° × original/observer-thinned.
- Independent original same-species + residual climate variance audit: scripts/analysis/fcp4475_climate_region5_vs_region1_paired_audit_20261024.py
- Independent audit CI: https://github.com/zuizui0223/fcp/actions/runs/40057276809 — PASS.
- Original photo/raster source zip SHA and artifacts are immutable, separately recorded in the full machine receipts under results/fcp4475_climate_within_region_verified_20261024.json and results/fcp4475_climate_paired_region_scale_verified_20261024.json.
- Source run artifact ID 12545591935; independent audit artifact 12545748040. Original photo IDs/coordinates/observer identities are not uploaded in these derived per-species ledgers.

**Manuscript-level wording:** Within-species nine-hue photographic compositions exhibit a small temperature-related signal across coarse local regions, but the specific directional signal toward the photographic white palette is not robust to simultaneous 1° block, month, independent-observer and multiple-predictor scrutiny. This supports caution about treating the global spatial temperature gradient as a general climate-maintained genetic flower-colour polymorphism. Independent floral organ and genotype/fitness validation remains necessary.
