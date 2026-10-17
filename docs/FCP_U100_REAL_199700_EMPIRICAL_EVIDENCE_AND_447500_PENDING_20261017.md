# FCP U100 flower-colour analysis — observed 199,700-photo evidence and pending 447,500-photo replay

**Evidence freeze:** 2026-10-17 (Japan Standard Time). **Project:** FCP within-species flower-colour polymorphism; **not** MOIS seasonal colour–climate research. This document distinguishes observed biological results from source availability and unfinished image inference.

## Status: one study, unchanged original estimators

- **Original measured:** 1,499 distinct species, 149,900 authentic photographs; discovery 500, validation 500, prospective 499. These rows are byte-SHA pinned and untouched.
- **Wave4 genuinely measured:** 498 new distinct species, 49,800 authentic original H9 photo assignments with full location-blind original ROI-v4/EfficientSAM terminal inference. The first 199,700 real photo rows and their original H1/D, species-conditional 50-km and fixed white-axis H2/W estimates are released permanently.
- **Final H9 opportunity:** 4,730 taxon IDs. The original independent full 100-photo source rule initially selected 4,478 species (447,800 source *row occurrences*). A cross-five-cohort source audit found **three re-used photo IDs and one re-used observation ID** in three new-U100 taxa, overlapping the frozen discovery cohort. Exclude those three *entire* 100-photo new taxa from final scientific pooling, without altering the original measured discovery rows, replacing photos, modifying the classifier or weakening the exact 100/source-taxum rule. The resulting **final analysis-eligible source** contains 4,475 species × 100 = **447,500 photo rows**. The other 255 original frame taxa consist of 252 H9 source under-100 and three cross-cohort duplicate-source exclusions.
- **New U100 real ROI pending:** the 2,481 newly H9-qualified species are being measured from their full frozen 248,100 original source photo assignments in [the production ROI #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400). The above three exclusion decisions are applied *only during the final source-verified five-cohort join*, not to the running original image measurement. New-U100 analysis contribution after these exclusions will be 2,478 species / 247,800 real photo terminals **if all frozen image jobs succeed**. Do not claim this contribution has been measured until the terminal job completes.

The earlier [raw 447,800-source release](https://github.com/zuizui0223/fcp/releases/tag/fcp-original-U100-4730-reconciled-source-20261017) is a faithful historical source snapshot, **not** proof that all five cohorts contain 447,800 globally distinct source IDs. The independently audited [447,500-analysis-eligibility release](https://github.com/zuizui0223/fcp/releases/tag/fcp-original-4730-independent-447500-analysis-eligibility-20261017), the [actual five-cohort audit #39015322314](https://github.com/zuizui0223/fcp/actions/runs/39015322314), and [source-collision documentation](fcp_u100_4730_original_source_identity_correction_20261017.md) supersede that broader uniqueness claim.

## Genuine four-cohort empirical results — old 149,900 + wave4 49,800

**The 199,700-photo original-method analysis has actually run and been permanently released**: [audited original measurement, H1, 50-km and W archive](https://github.com/zuizui0223/fcp/releases/tag/fcp-real-original-199700-wave4-50km-H1-W-20261017). Its [original ROI #38992543115](https://github.com/zuizui0223/fcp/actions/runs/38992543115) completed all 256/256 true source-image measurement workers and the historical inference, without worker failures. The append/50-km replay [#39023781282](https://github.com/zuizui0223/fcp/actions/runs/39023781282) succeeded and reproduced all three old frozen 50-km results exactly.

### 1. Local clustering of observed within-species flower-colour labels

The exact historical statistic, held fixed, is the species-equal-weighted mean of **global pairwise colour dissimilarity minus within-radius local pairwise colour dissimilarity**, with its within-species complete-vertex label permutation (199 null assignments), the frozen 40 classifiable-photo and 30 local-pair eligibility conditions, and radius = **50 km, post hoc**.

| Cohort | Species satisfying true 50-km pair gate | Mean dissimilarity depletion | Upper Monte Carlo p |
|---|---:|---:|---:|
| Discovery (unchanged original) | 166 | 0.020529254583812922 | 0.005 |
| Validation (unchanged original) | 181 | 0.018672971642749295 | 0.005 |
| Prospective (unchanged original) | 204 | 0.01468491968437185 | 0.005 |
| **Wave4 (new authentic 49,800 pixels)** | **190** | **0.017293596294403593** | **0.005** |
| **Four-cohort pooled (199,700 pixels)** | **741** | **0.01763720788097046** | **0.005** |

For the pooled 741 evaluable species, the observed mean global pairwise D was **0.2342543723894402** and local pairwise D was **0.21661716450846974**. The original null distribution of pooled mean depletion had mean **0.00006116024056826038** and central 95% quantiles **[-0.0033654692419269484, 0.0037227124164020096]**. The across-species association of pairwise D with local depletion was **rho = 0.09640942518477592**, upper permutation p = **0.005**: statistically detectable in this particular reference null but **small in magnitude**.

**Supported interpretation:** Spatially nearby flower-colour states are more similar than these same measured states randomly reassigned *within species*, across the original four sampled cohorts. This does **not** demonstrate adaptive habitat matching, direct selection by climate or pollinators, genetic differentiation, true population colour frequencies or independent replication across photo providers. Fifty kilometres was chosen post hoc, not a newly preregistered predictive distance; the minimum attainable upper p with 199 permutations is 0.005.

### 2. Observation-independent reliability of the species-level D estimate

The unchanged historical observer-disjoint H1/D reliability analysis was applied to the wave4 real-photo cohort. Across **200 observer-disjoint splits**, primary gate (20 classifiable photographs per half):

- Median Spearman rho across taxon-level paired D values: **0.8055963704422957**; 5th and 95th percentiles: **0.7754570576905138, 0.8327328313702286**.
- Median paired species per split: **316** (minimum 307, maximum 325); median concordance correlation coefficient: **0.8441405262228767**.
- Median Spearman–Brown reliability: **0.8923327313096993**.

**Supported interpretation:** Relative photographed species-level colour diversity is substantially reproducible across disjoint observers. The observer effect is therefore not sufficient by itself to explain all between-species variation, **but** these photograph-sample statistics do not estimate wild population morph frequencies without further detection-bias calibration.

### 3. Repeated alignment with the fixed white-versus-other colour axis

The unchanged H2/W statistic is **squared axial alignment** of continuous nine-colour flower-composition differences with a frozen white-vs-eight-other-colours zero-sum axis, assessed against the **original 999 structured, coarse-morph-preserving null permutations**. The square makes W sign-invariant. Interpret this as **white-axis alignment, never as net movement toward white flowering**.

| Preset gate | Species with usable continuous vectors | Observed W | Structured-null median | Upper p |
|---|---:|---:|---:|---:|
| Primary second-mode threshold 0.10 | 134 | **0.5368155693013797** | 0.43625410562222344 | **0.001** |
| Strict second-mode threshold 0.20 | 67 | **0.5472800875548313** | 0.4580858348292451 | **0.001** |

**Supported interpretation:** In these selected within-species colour contrasts, the fixed white-versus-other compositional axis explains more squared directional alignment than expected under the specified structured null. W **does not prove more white flowers**, their evolutionary direction, ancestor/descendant state, a re-expression mechanism, or a conservation of pigment synthesis cost. The original 999-null minimum upper p is 0.001.

## The still-open main question: final five-cohort 447,500 real-photo test

These results refer to **199,700 already measured photographs and 1,997 sampled species**; they are **not** the final five-cohort 447,500-photo empirical result. The original-U100 new ROI run needs to finish its remaining 256-worker partitions and original terminal inferential step. The [five-cohort final workflow](https://github.com/zuizui0223/fcp/actions/workflows/fcp-all4730-original-pixel-join-statistics-20261017.yml) has a gate requiring both authentic 256-worker real ROI runs before:

1. Authenticate all three unchanged source SHAs; append completed wave4 and new-U100 pixel terminals.
2. Remove whole three flagged duplicate-source new taxa (taxon IDs **62367, 75633, 1045161**) *after* source-photo terminal SHA verification, resulting in **4,475 × 100** final rows.
3. Reproduce exactly the three historical 50-km effects and redo the unchanged species-equal-weighted original vertex-permutation null over all five cohorts. Report the **true 50-km local-pair eligible n**, not the broader n-classifiable>=40 photo support.
4. Run unchanged observer-disjoint D reliability and fixed-white-axis W, then SHA-archive verified real combined colour measurements and biological results.

Do not transfer the wave4 W p=0.001, wave4 50-km p=0.005, or four-cohort pooled effect to the unfinished full U100 new cohort. The final effect can change; that is what adding 2,478 independent source-eligible species is intended to measure.

## Copy-ready Results language — empirically supported interim, not final U100 result

> The original within-species floral-colour analysis was extended from 1,499 taxa (149,900 measured photographs) to a fourth, taxon-disjoint cohort of 498 taxa (49,800 additional authentic photographs) using unchanged outcome-blind H9 sourcing, ROI-v4/EfficientSAM flower isolation and original inferential definitions. Local within-species colour dissimilarity at 50 km was lower than the corresponding global pairwise dissimilarity in the new cohort (mean depletion 0.01729; 190 evaluable taxa; species-conditional vertex-permutation upper p=0.005) and in the pooled four-cohort dataset (mean depletion 0.01764; 741 evaluable taxa; upper p=0.005). Species-level colour diversity showed observer-disjoint reproducibility (median Spearman rho 0.806 across 200 partitions), and continuous within-species colour differences aligned disproportionately with the prespecified white-versus-other colour axis (W=0.5368 versus structured-null median 0.4363; n=134; upper p=0.001). These outcomes characterize photographed phenotypic patterning rather than causal selection, evolutionary direction or population morph frequencies. A larger five-cohort U100 analysis is underway and is not represented by these estimates.

**Related machine-verified evidence:** [Actual 199,700-photo release](https://github.com/zuizui0223/fcp/releases/tag/fcp-real-original-199700-wave4-50km-H1-W-20261017), [wave4 original ROI and W inference](https://github.com/zuizui0223/fcp/actions/runs/38992543115), [real four-cohort null reanalysis](https://github.com/zuizui0223/fcp/actions/runs/39023781282), [five-source ID audit](https://github.com/zuizui0223/fcp/actions/runs/39015322314).
