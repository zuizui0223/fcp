# FCP: interannual local photo-colour contrast audit and verified candidate funnel

**18 October 2026 — exploratory observational evidence and next field decision.** This does **not** supersede the [original 447,500-photo FCP source](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017), the original post hoc 50-km spatial result, or the separately completed **negative** 5-km same-month spatial null (p=0.300). Their inferential roles remain distinct.

## Why this is more biologically informative than another global spatial effect

Species-range photographs can reflect patchy geography because observers sample different sites, because dispersal and drift produce autocorrelation, or because local phenotypes experience selection. **That does not establish genetically inherited morph coexistence or protected polymorphism in an interbreeding population.** A fair next question is whether two visually distinct coarse flower-colour states are observed *again* in a genuinely near-local observational setting across years, and whether the evidence survives basic checks against multi-coloured flowers and image/flower-age artefacts.

The [genuine-data workflow #39097882762](https://github.com/zuizui0223/fcp/actions/runs/39097882762) completed successfully with five synthetic tests. It acquired the immutable **447,500 actual original-photo terminals** (SHA256 of original measured CSV \`8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6\`) and used the **227,107 classifiable** four-state photograph labels. Small replay outputs are stored as the GitHub Actions artifact \`fcp-real-local-contrast-recurrence-and-photo-purity-audit\`, with output JSON SHA256 \`f5ac0c5e1b02bdec88275fc9eaecaef6f703610855470b46ab26f4610a273bce\` and anonymized candidate CSV SHA256 \`2899742b2bb737d9f36dd1f7907bb0257de1ba3b722530a949ebb598d1e47fae\`. Full reproducibility script: \`scripts/analysis/fcp_recurrent_local_photo_morph_candidates_20261018.py\`; source-locked compact Git result: \`results/fcp_fivecohort_recurrent_local_morph_candidates_20261018.json\`.

**No exact geographical coordinates, source observation IDs, photo IDs or observer IDs appear in the published species summary.** If the photographs include sensitive taxa, field sites should be resolved only in the original authorised source interface with adequate protections.

## Observational logic — actual independently observed photo pairs, not pair pseudoreplicates

For each species, a pair requires two original, independently identified observations, recorded by **different observers in the same year and calendar month**, with distinct photo IDs. The original recorded geographical accuracy of **each** observation must be at most 20% of the radius. Compare 0.25, 0.5, 1, 2 and 5 km actual geodesic distances. Each species is counted only once per radius, even if multiple photograph pairs exist.

A putatively recurrent site-colour contrast requires two *different years* showing the **same two coarse colour labels** in the **same calendar month** and the **same narrow area**. The four contributing original observation IDs and observers must all be distinct. The **maximum of all six geographical distances** among the four photo locations must remain within the stated radius. That complete-link requirement avoids claiming a local site from chained points that are farther apart. The site is only a photo-coordinate neighbourhood, **not a demarcated reproductive population**.

| Max pair distance | Taxa with any matching-month independent-observer local pair | Taxa with a different-colour photo pair | Taxa with same-month, same-colour-pair recurrence in two different years and four observers |
|---|---:|---:|---:|
| 250 m | 107 | 28 | **2** |
| 500 m | 176 | 49 | **3** |
| **1 km** | **287** | **91** | **8** |
| 2 km | 453 | 145 | **8** |
| 5 km | 801 | 300 | **17** |

These are *species with candidate photographic evidence*, not evidence that exactly 8 species have maintained genetic polymorphisms or that local FCP is biologically rare. At 1 km, only **47** species had at least five valid local pairs and **five** had at least 30; the photo atlas has sparse fine-locality temporal sampling by design. Absence of an eligible pair often means **no photographic opportunity**, not a truly single-coloured natural population.

## Source-linked first photo-review list (1 km complete-link; each with four independent observers)

These taxa require **manual flower-organ identity/colour-development review before even treating an image difference as two individual plant morphs**:

| Candidate species | Coarse photo contrast repeated in two years | Same calendar month | Example years | Best observed minimum ROI primary-state share (all four photos) |
|---|---|---:|---|---:|
| *Navarretia heterodoxa* | blue/purple – white | June | 2019 / 2026 | **0.732** |
| *Allium cratericola* | red/pink – white | March | 2021 / 2023 | 0.553 |
| *Boechera breweri* | blue/purple – red/pink | March | 2019 / 2024 | 0.552 |
| *Byrsonima lucida* | white – yellow/orange | April | 2022 / 2023 | 0.305 |
| *Castilleja densiflora* | red/pink – white | April | 2023 / 2026 | 0.540 |
| *Castilleja rubicundula* | white – yellow/orange | May | 2024 / 2025 | 0.608 |
| *Collinsia multicolor* | blue/purple – white | April | 2022 / 2023 | 0.598 |
| *Lomatium dasycarpum* | white – yellow/orange | March | 2019 / 2025 | 0.588 |

*Navarretia heterodoxa* also has a separately observed May example with the same photo colour contrast (2021 / 2024), but repeated photographic event sets and months **are not independent genetic populations**. The choice of illustrated year pair is a deterministic photo-quality/practicality screen, not a hypothesis test or selection of a favourable p-value.

## Critical ambiguity: bicolour corollas, flower age and ROI colour mixtures

Original flower-ROI output includes nine-component colour fractions. For each photo we summed palette fractions belonging to the **coarse assigned colour group** (white; yellow/orange; red/pink/magenta; blue/purple). We asked whether **all four** cross-year evidence photos exceed an *arbitrary exploratory* grouped-primary-colour fraction threshold; it is not a probability of correct morphology.

| Minimum assigned group share across each of four evidence photos | 1 km recurrent candidate species | 5 km recurrent candidate species |
|---|---:|---:|
| 50% | 7 | 14 |
| 55% | 6 | 10 |
| 60% | 2 | 7 |
| **65%** | **1** | **3** |
| 70% | 1 | 3 |
| 75% | 0 | 1 |

**Very consequential limitations:** Not one of the eight 1-km candidates had an available full four-observer, cross-year photo set in which every photograph had **<=20%** of its flower ROI in every other coarse colour group. In other words, the photos often include substantial pixels assigned to a second colour. This is compatible with a single multi-coloured corolla, local flowers at different ages, overlapping organs, or segmentation/lighting artefacts. It is *not* evidence against genetic polymorphism, but severely limits inference from the coarse classifier.

Independent natural-history information makes this a biologically relevant distinction:

- *Collinsia multicolor* is naturally bicoloured **within one corolla**: typically white upper petal lobes and purple/lavender lower lobes ([taxon description and flower photographs](https://www.calflora.org/app/taxon?crn=2283)). A switch between algorithmic "white" and "blue/purple" labels could follow which floral part dominates a photograph without *any* genetically different morphs.
- *Castilleja rubicundula* has white, pink, yellow and sometimes bicoloured flowers ([descriptive taxon account](https://en.wikipedia.org/wiki/Castilleja_rubicundula)); field photographic notes describe flowers changing from yellow to white to pale pink-purple during ageing in a local form ([original photo series](https://www.flickr.com/photos/mark_egger_castilleja/sets/72157674696410930)). This needs explicit **same-flower flower-age documentation**, not just species/hue record matching.
- *Allium cratericola* tepals are described as white, pink or purple often with contrasting dark midveins. Corolla patterning is itself another potential photo-state confounder ([flora account](https://ucjeps.berkeley.edu/eflora/eflora_display.php?tid=12550)).
- *Navarretia* has multiple distinct flower-organ/colour patterns, and cryptic taxa in related *N. divaricata* complexes have required taxonomic revision; colour labels alone are insufficient for identity/morph verification ([Johnson & Gowen 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5784239/)). The present *N. heterodoxa* candidate is **not automatically assigned to that other species complex**.

**Operational conclusion:** Most recurring "photo colour contrasts" should be routed first to botanically informed, blinded per-image scoring of **the same floral organ at the same developmental stage**, with a discrete distinction between one bicoloured flower, a transient anthesis/age change, different individuals, and observer/image artefacts.

## What would identify *evolutionary maintenance* instead of static photo-colour structure?

Even a perfectly field-measured **single cross-section of genetic morph frequencies** cannot generally distinguish local stabilising selection from migration/neutral structure. Consider a two-morph multi-population toy model for true **genetically inherited** frequency \`p_j\`:

\[
\dot p_j = p_j(1-p_j)\bigl[s_j + k(1-2p_j)\bigr] + m(\bar p-p_j).
\]

At an interior equilibrium (\`0<p_j<1\`), every observed vector of frequencies can be made an equilibrium for a **continuum** of local selection \`s_j\`, frequency-dependent stabilisation \`k\` and migration \`m\` by choosing

\[
s_j = -k(1-2p_j) - \frac{m(\bar p-p_j)}{p_j(1-p_j)}.
\]

This elementary *parameter nonidentifiability* is more fundamental than the photo classifier. Without separately measuring selection/frequency dependence, drift and gene flow, **static observed colours do not identify a mechanism of polymorphism maintenance**. Photographed non-genetic colour fractions provide even less identification.

The discriminating, feasible next work is therefore **not** merely to harvest additional photographs per species globally, but:

1. **Floral-identity gate:** Blind botanical review of original images at the same flower organ, flowering stage and light quality; record genuine candidate inter-plant colour morphs vs multi-colour organs, developmental changes, cultivar/record errors. Do not label candidates "confirmed FCP" before this.
2. **Longitudinal local population panel:** In at least two seasons, mark *individually distinct plants* within precisely delimited field populations, measure true colour morph frequencies and phenological stage, and document seed-pedigree/common-garden inheritance. Select sample size from pilot morph rarity and measurement power, rather than claiming photo-pair counts are independent replicates.
3. **Mechanism separation:** Pair frequency-specific fecundity/seed set and pollination with neutral SNP/MIG-seq population history plus source ecology; contrast frequency-dependent maintenance vs spatial divergent selection with gene flow vs drift/plasticity.
4. **Speciation gate:** Only when direct assortative pollen flow, cross-compatibility/hybrid fitness, and genetic isolation relative to geographic controls imply actual reproductive isolation should "incipient speciation" be considered. Different colour forms inside one species are not by themselves ecotypes or separate species.

The immediate field candidate count is **8 unverified taxa at 1 km**, sharply reduced under simple pixel-composition QC. The bottleneck is now **within-flower photo interpretation, individual identity and longitudinal demographic/fitness data**, not another global map of spatially autocorrelated photos.

## Separation from the active New Phytologist paper

The earlier [PR #243](https://github.com/zuizui0223/fcp/pull/243) contains a reproducible paper draft about species-wide photographed skew/segregation. It remains **held for biological scope review**: 50-km source-conditional sorting is well documented but **cannot yet be elevated to demonstrated within-population genetically maintained FCP**. This research-quality screen neither cherry-picks a new positive spatial p-value nor changes the separate original prospective white-axis W inference. The independent [PR #244](https://github.com/zuizui0223/fcp/pull/244) demonstrated that a stricter same-year-month 5-km null was **not supported** (n=194, p=0.300); preserve it unchanged.
