# FCP global flower-colour pattern — actual measured-photo results, 20 October 2026

**Status:** actual 447,500-image source verified and successfully analysed; exploratory real-image estimands. The original three-cohort prospective W, original frozen 50-km result and main branch were not modified.

**Primary run:** https://github.com/zuizui0223/fcp/actions/runs/39412659475

**Complete archived 65-cell numerical record:** ../results/fcp_world_pattern_real_447500_20261020.json

**First world figure (colour D + sampling opportunities):** figures/fcp_global_photographic_colour_opportunity_20261020.svg

## 1. What repeats among observed photographed intraspecific colour patterns?

4,475 taxa × 100 measured original photos = 447,500; 227,107 have an admissible four-state photo-colour label. Of 3,297 species with at least 40 classifiable images and three observers, 3,029 contain at least two distinct photo colour labels and 2,963 contain both white and nonwhite labels. This is photographed colour-label heterogeneity and NOT a natural-population genetic polymorphism frequency.

Photo-label counts: white 107,956; yellow/orange 70,134; red/pink 23,908; blue/purple 25,109. Species-equal photographed pair-diversity D=0.237514. White-containing pairs constitute an average 85.55% of discordant coarse photo-label pairs. This is a colour-composition-dependent descriptive statistic, NOT a fraction of species, W alignment, or evidence of molecular anthocyanin loss.

The earlier untouched prospective 158-species H2 result separately supports a repeated white-versus-chromatic direction: W=0.517246, fixed structured-null median=0.457143, p=.001. Do not conflate it with these new descriptive measures.

## 2. How is the photographed world spatially sampled?

In the 18 × 9 equal-area geographic grid, at least one classifiable photo occupies 96/162 cells. An eligible taxon-cell requires 10 classifiable photos and 3 different observers; 65 cells contain such a taxon and only 28 cells have at least 20 eligible taxa. Total: 6,704 eligible taxon-cell units across 4,112 taxa; 2,118 taxa span two or more such cells. Blank mapped cells are sampling gaps, not biological evidence of absent polymorphism.

Absolute latitude (species-equal D, with DIFFERENT taxon compositions in each band):
- 0–15°: 146 taxa, D=.33149
- 15–30°: 806 taxa, D=.29011
- 30–45°: 3,588 taxa, D=.28246
- 45–60°: 1,608 taxa, D=.26012
- 60–90°: 98 taxa, D=.24378

Within 246 identical taxa appearing in both tropical and extratropical zones, tropical-minus-extratropical D=.01213, taxon bootstrap 95% CI [−.00583, +.02981]. A species-invariant latitudinal diversity rule is unsupported. The ecological or historical causes of taxon replacement and sampling differences are untested.

## 3. How much do photograph colour-purity gates alter the apparent pattern?

A post hoc diagnostic requires ≥65% assigned coarse-colour group within segmented floral ROI, with no competing group >20%. This retains 162,003 classifiable photos; retention differs by colour class: white 80.37%, blue/purple 57.94%, red/pink 60.31%, yellow/orange 65.98%.

Among exactly the SAME 1,813 species eligible in both original and strict-purity collections, species-equal D drops from .11544 to .06468, paired mean −.05076 (95% bootstrap interval [−.05339, −.04821]). Species with two or more photo labels drop from 1,557 to 1,054. The photo-quality gate itself is class-biased and the captured flower morphology remains unverified. It would be incorrect to compare 3,297 broad species directly to 1,813 strict-purity species and call the total difference a quality effect.

## 4. Does geographic sorting exceed the sampling-and-colour-composition null?

Fixed species, same eligible regional cells, exact state composition and cell photo counts; 199 species-conditional photo-label permutations; species-equal average total variation across pairs of geographic cells. This new comparison is post hoc and its minimum attainable Monte Carlo p is .005.

| Same 1,282 taxa in identical eligible cells | Broad photos | Strict-purity photos |
|---|---:|---:|
| Observed between-cell total variation | .074785 | .055682 |
| Expected under matched vertex null | .066397 | .048409 |
| **Excess** | **+.008387** | **+.007273** |
| One-sided upper Monte Carlo p | .005 | .005 |
| Excess absolute white-frequency difference | +.006308 | +.006507 |

For all 2,118 photo-supported taxa in 2+ eligible cells, observed TV=.110117, null=.096441, excess +.013676, exploratory p=.005.

**Defensible scientific statement:** source iNaturalist photographic flower-colour composition contains a small geographically nonrandom component beyond sampled colour composition and geographical photo-count allocation. This persists under a fixed strict photo-colour-purity screen in an exactly matched 1,282-taxon, same-grid comparison. It does not prove within-population polymorphism, drift, genetic IBD/IBE, demographic colonisation history, divergent selection or speciation.

## 5. New explicit geographical colour-axis decomposition (actual source, further audits)

[The completed CI run #39422571304](https://github.com/zuizui0223/fcp/actions/runs/39422571304) used all authentic source photographs and retained the five authentic species-disjoint sourcing cohorts. All 199-source-vertex nulls and robustness statistics were recalculated; these are **post hoc exploratory additions**, not new prospective confirmation.

For four-class cell frequency vectors, the between-cell total variation has an exact nonnegative decomposition:

\[
TV(p_a,p_b)=\left|p_{a,\mathrm{white}}-p_{b,\mathrm{white}}\right|
+\frac{1}{2}\left(\sum_{k\ne\mathrm{white}}|p_{a,k}-p_{b,k}|-\left|p_{a,\mathrm{white}}-p_{b,\mathrm{white}}\right|\right).
\]

The first term captures the observed geographic white-versus-rest photo-label contrast. The second captures nonwhite colour redistribution after accounting for the required compositional balance. **These are not biochemical anthocyanin or phenotypic-vector-W measurements.**

| Same 1,282 taxa and same sampled regional grid cells | Broad photographs | Strict ROI-purity photographs |
|---|---:|---:|
| Total TV excess after source-composition-preserving vertex null | +0.008387 | +0.007273 |
| White/rest component excess | +0.006308 | +0.006507 |
| Residual nonwhite composition excess | +0.002079 | +0.000767 |
| Ratio of white-component excess to total TV excess | 75.2% | 89.5% |
| Matched upper p white | .005 | .005 |
| Matched upper p nonwhite residual | .005 | .085 |

The strict-purity 89.5% ratio measures the **relative allocation of a small observed-minus-null excess**; the bootstrap 95% interval 72.1–104.5% is unstable near a small denominator and is *not* a percentage of flowers, species or phenotype contrasts.

### Strong heterogeneity by independent taxon sampling cohort

| Cohort | Strict-purity eligible taxa | Mean excess white/rest | Original 199-world upper p |
|---|---:|---:|---:|
| Discovery | 127 | +.005351 | .035 |
| Validation | 135 | −.000934 | .645 |
| Prospective | 155 | +.002071 | .305 |
| Wave4 | 147 | +.001911 | .310 |
| New U100 | 718 | +.010009 | .005 |

Leaving out New U100 leaves 564 taxa, mean excess +.002048, upper p=.14. Hence the pooled 89.5% value is **not a consistently cross-cohort-replicated general law**. It is largely supported by the new, larger source cohort. The source cohorts are distinct *species*, but all share the same iNaturalist image acquisition and ROI model.

### Matched-null robustness and genuine photo-colour opportunity

In the whole 1,282-species strict subset, matched-null mean upper p=.005, matched-null 10%-trimmed mean upper p=.01, but the median taxon white excess equals 0 and does **not** exceed its null (p=1.0). White/excess dominance in the species-equal mean is therefore not typical of all photographed taxa.

Eligibility requiring at least five high-purity white-labelled and five nonwhite-labelled photographs **within eligible grid cells** leaves 127 taxa. In that taxon set, white-frequency component excess averaged +.05566 and had median +.01928; both mean and median matched permutation p=.005. This is a **selection on documented colour-sampling opportunity**, preserved under the within-taxon vertex null, not proof of genetic morph coexistence. The cohort counts in those 127 taxa are 9/11/24/11/72 for Discovery/Validation/Prospective/Wave4/New U100.

Excluding New U100 from that colour-opportunity-selected set leaves 55 taxa: mean white component excess +.03027, permutation p=.010, median +.00973, permutation p=.095. The remaining 55 are not sufficiently broad for a universal taxon law, and cohort-specific n=9 and n=11 are below the declared descriptive gate of 20.

**Interpretation ceiling:** white/nonwhite geographic photographic sorting is a repeatable *population-sample phenotype candidate* in several source cohorts and in a strongly colour-supported subset, but the broad claim that all plants repeatedly differentiate genetically on the white axis is not established. The white frequency indicator and prospective continuous-colour W result are distinct measures.

### Required further correction: compare each component relative to its own null baseline

The 89.5% strict-photo relative share above describes the white/rest allocation of an **absolute small geographic excess**. It is **not a claim that the white component has preferential spatial sorting relative to the nonwhite component**. A parallel FCP development [PR #248](https://github.com/zuizui0223/fcp/pull/248) already warned that 50-km white–chromatic dominance can arise because white-containing contrasts dominate the available photo-label composition; the earlier component-baseline-corrected 50-km result did **not** support a white-specific spatial clustering advantage.

In this branch, we therefore made the direct matched-taxon relative comparison at a second *equal-area interregional* scale. For every eligible species and colour component we divide the observed geographic TV component by that same component's expected TV under the original species-vertex-label null, then subtract 1. The two relative excesses are compared within identical taxa. We require a nontrivial **null expected white and null expected nonwhite residual ≥0.01** to avoid division by near-zero opportunity; 0.005 and 0.02 are explicitly post hoc sensitivity gates. The exact null worlds, species weighting, observed species/photo composition and source cohort are unchanged.

| Regional scale component-normalized comparison | Broad photos | Strict-purity photos |
|---|---:|---:|
| Eligible taxa, BOTH component null expectations ≥0.01 | 395 | 174 |
| Mean relative excess of white component | +0.1425 | +0.1747 |
| Mean relative excess of nonwhite residual component | +0.1150 | +0.0813 |
| White minus nonwhite relative excess | +0.0276 | +0.0934 |
| Matched vertex-null upper p for white component stronger | **0.370** | **0.220** |
| Species-bootstrap 95% CI for difference | [−0.1294, +0.1876] | [−0.1160, +0.3186] |

Thresholds 0.005 and 0.02 do not lead to a supported white-specific relative excess in either source-quality collection (all matched p ≥0.125). The direction even becomes negative for broad photos at the 0.02 gate. Among strict-quality taxa with baseline ≥0.01, omitting the new U100 species gives n=77, relative difference −0.1045 and p=.755.

**Verdict: white-specific proportional geographic enrichment is NOT SUPPORTED.** This result is consistent with the original 50-km component-baseline-normalized non-support on PR #248. Our pooled 89.5% white fraction must not be used in the abstract as evidence of preferential spatial selection for white flowers, repeated anthocyanin loss, molecular accessibility, cost–benefit advantage, or parallel adaptation. It can appear only as a clearly labelled *photographic compositional decomposition* with component-specific baseline and quality sensitivities.

[Source-locked numerical result](../results/fcp_world_white_axis_spatial_decomposition_real_20261020.json) — [successful exact-original-photo Actions #39423430999](https://github.com/zuizui0223/fcp/actions/runs/39423430999).

**Important remaining distinction:** the independent, species-disjoint prospective white/chromatic **phenotype-space W** confirmation is unchanged and supported. The independent empirical predictions "the white phenotype vector is common" and "white morphs are preferentially geographically structured" are NOT equivalent. The latter has failed this direct component-normalised test at present.


### Immediate evolutionary follow-up, not yet done

The next FCP empirical model should condition colour reassignment on calendar month (and separately year-month) and independent observer/locality availability, then compare **white-axis sorting effect heterogeneity across native taxa**. WorldClim BIO5/BIO14 are suitable named exploratory environment covariates, but no new environment association or historical demographic reconstruction was estimated by this run. For true history-versus-environment attribution, neutral SNPs, reconstructed colonization/connectivity, accessible habitat M, and documented heritable colour morphs remain required.


## 6. Manuscript priorities and process hypotheses

**Figure 1:** continuous white–chromatic phenotype axis from frozen prospective H2 plus new global source-photograph coverage. Ensure colour counts and W are separate quantities.

**Figure 2:** global photos show small geographically distributed intra-specific phenotype contrasts; always show the cell eligibility map and pure-photo sensitivity. The 5-km same-month strict null was previously not supported, so no claim of widespread locally maintained inherited FCP.

**Figure 3:** current photo-only IBD-like > BIO5 phenotypic IBE-like in the earlier three cohorts (mean partial rho approx .022–.028 versus .008–.010). However, the BIO5 signal fails strict background/observer controls, and the larger geographic component is not independent evidence that historical range shifts dominate. Test environment-matched, geography-different pairs and geography-matched, environment-different pairs with new biome/soil/dynamic climate annotations, spatial holdouts, photography and native-range controls.

**Figure 4:** for subsets with independently available neutral SNP/phylogeography, add ancestry, landscape migration, and post-LGM habitat connectivity to quantify history, rather than naming present-day distance alone as history.

**Figure 5:** treat floral anthocyanin expression, petal-restricted loss (PAL), and whole-plant anthocyanin loss (WAL) as an independent molecular hypothesis. Silene littorea F3h1 / SlMyb1a transcript, HPLC and PAL/WAL population-frequency results are published sources (Casimiro-Soriguer et al., 2016, DOI 10.3389/fpls.2016.00204; Del Valle et al., 2019, DOI 10.1186/s12870-019-2082-6), NOT molecular tests of these 4,475 photographic species. Direct species-level genotype and pigment profiles, fitness and local morph inheritance remain missing.

**Claim ceiling and next empirical decision:** the photographed global landscape and robust local-to-regional photo differences can be mapped now. Only wild/native, true floral-part/flower-age-validated data should be analysed as natural geographic colour sorting. IBD and phenotypic IBE are correlational; neutral genetic isolation, range history, and environmental causal fitness must be independently measured before assigning their relative importance.
