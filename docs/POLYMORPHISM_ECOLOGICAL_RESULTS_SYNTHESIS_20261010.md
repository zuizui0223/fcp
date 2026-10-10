# What the FCP atlas actually discovers ecologically: locally retained colour variation with non-universal geographic sorting

**Ecological-results synthesis for the New Phytologist paper — 2026-10-10.** This document does not add a new sensitivity analysis, benchmark or photo-colour measurement. It revises the ecological interpretation of **existing, source-verified empirical observations**, including completed draft PRs and current main. Claims in a draft PR are explicitly labelled as such and must not silently acquire confirmatory status.

## Result 1. Geographic sorting is replicated but modest compared with locally observed colour discordance

The three original cohorts are species-disjoint and used the same iNaturalist photo source and classifier. Each species' exact photographed four-state composition is fixed by its matched geographical permutation null. The empirical *species-equal pairwise* values at <=50 km are:

| Source cohort | Number of locally evaluable photographed species | Source species-wide colour-pair discordance | Observed within-50-km conspecific photo-pair discordance | Difference (species-wide - local) | Relative deficit vs mean species-wide discordance |
|---|---:|---:|---:|---:|---:|
| Discovery | 166 | 0.2213947103 | 0.2008654557 | **0.0205292546** | **9.27%** |
| Validation | 181 | 0.2286621608 | 0.2099891891 | **0.0186729716** | **8.17%** |
| Third | 204 | 0.2482035980 | 0.2335186783 | **0.0146849197** | **5.92%** |

The previously reported permutation p=0.005 in all three establishes a nonrandom **modest deficit**, not that the majority of within-species flower-colour ITV is assigned among monomorphic local populations. **Most of the baseline pair-discordance magnitude is still visible within the 50-km photo neighbourhoods**: local / range-wide paired difference is approximately 90.73% / 91.83% / 94.08%, respectively. These fractions are ratios of equal-species averaged PHOTO-PAIR discordances, not a formal genetic variance partition and not true frequencies of mating-population mixed morphs.

This distinction is the central ecological improvement. The source-supported pattern is **local colour heterogeneity with repeatable nonrandom geographic sorting layered on top**. It does not justify replacing all local diversity with an among-population partition or calling 50km neighbourhoods mostly monomorphic. The field-level result remains an evolutionary question: are nearby differently coloured flowers distinct genets, a plastic flowering response, misidentified floral organs or observational artefacts? That is not resolved.

**Source:** main source JSON results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json; current main manuscript. Direct month-conditioned geographically distributed photo-label structure remains positive in each source cohort in draft PR #129 (0.017455 / 0.015287 / 0.013324), and the nonwhite-only four-state reduction remains positive (0.01947 / 0.00729 / 0.00891) in main. These are support for an observational geographical pattern, not a field genetic mixed-population frequency.

## Result 2. Repeated geographic structure does not have one repeated geographic direction

Species-fixed absolute-latitude/elevation photographic comparisons across five absolute-latitude and four coarse elevation bands, with maxT corrections, yielded **no 3-cohort replicated same-sign white-frequency or four-state D cline** [draft PR #127](https://github.com/zuizui0223/fcp/pull/127). Both very low absolute latitude 0–15° and high 60–90° were insufficiently sampled for that stronger within-species inference. Therefore the supported finding is not "higher latitudes have more white flowers", "mountains select coloured morphs", or their reverse.

A separate source-matched within-region geographic analysis [draft PR #130](https://github.com/zuizui0223/fcp/pull/130) established a **coverage-qualified zone of positive photographic sorting**:
- In 23.5–60° absolute latitude, 152 / 169 / 197 within-region species from discovery / validation / third had the opportunity for <=50km and >50km within-region photo comparisons, while 126 / 142 / 167 had within-month label exchangeability. The month-conditional mean local-depletion residual is **+0.0191175 / +0.0157416 / +0.0114338** and passed the PR's fully corrected 3-cohort regional support gate.
- The narrower 30–60° band has positive residual point estimates +0.0217047 / +0.0162487 / +0.0093478, but **did not** meet the joint all-cohort Holm/observer support gate. It must not be described as fully replicated, even though the raw estimates have the same sign.
- The 0–23.5° tropical band has only 12 / 7 / 5 within-region geographically evaluable species, and 60–90° only 1 / 1 / 2. These are coverage HOLDs. We cannot claim the pattern is absent in the tropics/polar zones, nor assign worldwide prevalence from the 42,111 outcome-blind photo opportunity frame.

This combination is ecologically informative: photo-geographic sorting **can be recurrent in the sampled well-supported extratropics without a universal latitudinal colour direction**. It is not proven uniformly distributed across the planet; the original photos are Northern-Hemisphere-heavy. Distinguish **the existence of within-species spatial organisation** from **shared maps or shared signed environmental effects**.

## Result 3. The spatial footprint is stronger than the abiotic candidate, and the candidate does not prove selection

The current main manuscript's continuous nine-colour within-species decomposition yields, for discovery / validation / third:
- Geographic-distance-associated partial rank correlation conditioning BIO5 difference: **0.02760 / 0.02414 / 0.02225**, each p=0.005.
- BIO5-difference-associated partial rank correlation conditioning geographic distance: **0.00925 / 0.00826 / 0.01019** (nominal p 0.005 / 0.010 / 0.005).

The distance term is descriptively **~3.0 / 2.9 / 2.2 times larger**, but rank correlations are NOT variance percentages or causal proportions. In particular the BIO5 flower-minus-background control and strict observer exclusions fail to establish robust, flower-specific, observer-independent thermal sorting.

More recent separate 42,111-source ecology work must not be misrepresented as adaptation: [draft PR #148](https://github.com/zuizui0223/fcp/pull/148) found small positive heldout photo-colour classification Brier predictive increments when adding environmental blocks to original-species training mean plus an explicit 250km or 1000km lower-rank photo-geodesic spatial field (climate 20,546 photos +0.001254 / +0.000753; soil-complete 11,136 +0.003076 / +0.002333). But the neighbourhood residual spatial Moran-like statistic was not uniformly lowered, and the other exact-source 649 actual phylogenetic-tip check did not show a clear full-environment gain. [Draft PR #147](https://github.com/zuizui0223/fcp/pull/147) checked local matched original same-species photo pairs: **none** of the individual temperature, rainfall, solar, vapor/wind or ten soil features survived 25-/35-predictor maxT FWER at 100/250/500km. Fifty km was HOLD for only 22 original colour-discordant pairs in the climate-complete nearest-pair set.

**Ecological conclusion:** there is a repeatable geographical colour phenotype footprint; a **common temperature-direction rule** and a **robust fine-scale climatic selection gradient** are not supported. Remaining hypotheses include drift/founder history and restricted dispersal, lineage/locale-specific colour-related selection, ecological plasticity and technical colour confounding. The data do not rank the true adaptive-maintenance contributions of these processes.

## Result 4. Common phenotype direction is not common ecological map

The frozen achromatic-versus-chromatic displacement target was confirmed in the third, prospectively sampled species/photo-disjoint source cohort: n=158, W=0.517 versus structured null median 0.457, p=0.001. This is recurrence of a **measured phenotype-space direction beyond the four-state construction null** in a shared photo measurement system. It does **not** imply a white-pigment allele, a universal warm-to-white cline, or local adaptation.

The comparative synthesis is therefore two-axis:
- **What variation is available?** Within-species photographed colours recurrently vary along an achromatic–chromatic dimension.
- **Where does that variation occur?** The range-wide photograph-derived phenotypes have nonrandom local geographical allocation, but without one signed global latitudinal or thermal rule.

This is a stronger ecological question than seeking one universal pigment payoff: **How can a common flower-colour phenotype axis generate locality-specific frequency mosaics without globally transportable signed climate optima?** Restricted dispersal/history, species-specific biotic regimes, physiological trade-offs and nonadaptive plasticity are competing explanations. It remains a *mechanism hypothesis*, not the observed comparative mechanism.

## Interpreting evolutionary prior literature responsibly

The geography of FCP is not new in itself. Narbona et al. (2018, Plant Biology, DOI 10.1111/plb.12575) reviewed Mediterranean species composed primarily of monomorphic populations with fewer polymorphic ones, generating clinal/mosaic distributions. Source: https://pubmed.ncbi.nlm.nih.gov/28430395/.

Iris lutescens provides a source-based example where geographic floral morph structure involved neutral drift/restricted gene flow and where neighbouring genetic/pigment evidence was necessary (DOI lookup https://pubmed.ncbi.nlm.nih.gov/27084922/). More recent field/genetic work in Castilleja coccinea (2025, https://pubmed.ncbi.nlm.nih.gov/40931921/) shows an example of herbivory/habitat-related trade-offs in genetically confirmed flower/bract colour polymorphism. These individual-species findings support biological plausibility, not a global explanatory model tested here.

**The present study's distinct contribution** is source-specified, quantitative **cross-species excess geographic sorting CONDITIONAL on each species' own colour frequencies**, plus prospective confirmation of a recurrent *phenotype-space* axis. It must not claim to have discovered FCP clines, mosaics, or basic habitat-dependent morph variation for the first time.

## High-value manuscript decision

Retain the title "Species-wide flower-colour variation is geographically partitioned across plant species" if desired, but qualify that "partitioned" means **a repeatable 5.9–9.3% relative *photo-pair discordance deficit***, not most variation is among monomorphic populations.

**Main ecological result paragraph:** visible geographic sorting is stronger than expected for each source taxon's source composition, but local photograph-level colour diversity remains substantial. The result is about **local mixing and geographic sorting coexisting**, not one replacing the other.

**Ecological Discussion paragraph:** the photographed phenotype-space axis recurs across taxa whereas no universal direction of latitude/elevation/temperature response is supported; within-species local allocation is empirically best constrained for sampled 23.5–60° latitude. Absent genetic and actual colour-specific selection data, drift/history, local biotic regimes and phenotypic plasticity are alternative ecological processes with no reliable worldwide ranking. Environmental heldout predictive value must not be relabelled as ecological adaptation.

**Paper scope:** The frozen original main H1/H2 numerical receipts, original figures, independent photo-source measurement states and unopened future 2,000+730 taxa are never modified by this synthesis. No more QC pipelines or methodological appendices are being requested by this scientific redirection.
