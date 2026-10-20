# FCP: within-versus-between regional BIO5–flower-photo association

**20 October 2026.** Exploratory analysis of authentic original 447,500 photo terminals, all from the same iNaturalist source. This is not a causal evolution result.

## Current inference

The same-species, same-calendar-month positive BIO5/white association depends strongly on the spatial scale at which comparisons are made. The new analysis exactly decomposes each qualified 5-degree species/month group's temperature contrast into a component within nested 1-degree subcells and one between subcells. A geographically structured observational association is not evidence that past range shifts, migration, neutral drift or local selection caused it.

Source photo CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6.
WorldClim BIO ZIP SHA256: 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf.
Verified initial Actions: https://github.com/zuizui0223/fcp/actions/runs/39430364175.

## Actual source results

| Same species, same month and one 5-degree cell | Broad classified photos | High-purity photos |
|---|---:|---:|
| Qualified species | 1,074 | 327 |
| Qualified taxon × cell × month strata | 1,673 | 396 |
| Total white minus nonwhite BIO5 difference | +0.19420 C | +0.32985 C |
| Among nested 1-degree subcells | +0.14477 C | +0.27158 C |
| Within nested 1-degree subcells | +0.04943 C | +0.05827 C |
| Between-1-degree original vertex-null p | 0.005 | 0.005 |
| Within-1-degree conditionally permuted p | 0.015 | 0.075 |
| Between component species-bootstrap 95% CI | +0.07141 to +0.22106 C | +0.12305 to +0.42023 C |
| Within component species-bootstrap 95% CI | +0.00774 to +0.09195 C | -0.02294 to +0.14621 C |

The covariance identity is exact: within the same five-degree month stratum, white-minus-nonwhite temperature difference equals average within-one-degree colour/temperature covariance divided by the original white-label binomial variance, plus covariance of one-degree subcell means divided by the same variance. Each decomposition is computed from the same photographed individuals, with 199 fixed-composition permutation worlds and species-equal aggregation.

### Interpretation guardrails

- Long-term WorldClim 1970–2000 temperature normals are not temperatures on flower observation dates.
- The original broad and strict subsets differ in taxonomic and photo support. Comparing them requires identical species × five-degree cell × month sampling opportunity; the separate matched-source workflow handles that.
- The photographic labels are not genetically inherited morph assignments. White-flower origin, pigment cost, local fitness, neutral genomic ancestry and demographic history remain unmeasured in these individual photographs.
- All source cohorts use the same imagery and colour measurement pipeline; source-disjoint species do not imply source-independent replication.
- A positive between-one-degree component may arise from historical population structure, contemporary habitat gradients, phenology/observer sampling or spatially structured plasticity. The current data do not identify which.

## Exact matched-quality reanalysis (now completed)

Without changing the biological stratum definitions, the source-authenticated rerun retained the same **327 species and 396 taxon x 5-degree cell x calendar-month strata** in both original and high-purity photo datasets. The retained high-purity photographs are a nested subset of the original observations, not another independently photographed study.

| Identical taxon / 5-degree cell / month opportunity | Original classified photos | High-purity photos |
|---|---:|---:|
| Entire BIO5 white-minus-nonwhite difference | +0.26495 C | +0.32985 C |
| **Between nested 1-degree cells** | **+0.21184 C** | **+0.27158 C** |
| Within nested 1-degree cells | +0.05311 C | +0.05827 C |
| Between-component 5-degree vertex-null upper p | 0.010 | 0.005 |
| Within-component 1-degree-conditioned upper p | 0.090 | 0.075 |

The difference between strict and original BIO5 effects among identical taxa is +0.06490 C for the total, +0.05974 C between nested 1-degree cells and +0.00516 C within nested 1-degree cells. However species-bootstrap 95% intervals all cross zero: total [-0.02246,+0.16193], between [-0.02341,+0.14747], within [-0.05005,+0.05978] C.

**Interpretation:** Once species and site-month observational opportunity are fixed, the between-one-degree component still dominates in **both** photo-quality representations, while the one-degree conditional effect remains unsupported. The changes caused by the strict photo-purity selection are smaller and uncertain, not a confirmed improvement in any morphological correctness. The quality filter may also change which 1-degree subcells contribute inside the matched 5-degree strata.

**Permanent real-data receipt:** [nested spatial-climate full pooled/matched source summary](../results/fcp_nested_1deg_5deg_climate_real_20261020.json). Verified actual source Action: https://github.com/zuizui0223/fcp/actions/runs/39431312956.

These outputs test geographic distribution of photographic labels and **cannot** establish that range history (past refugia, expansion or restricted genetic dispersal) dominates climate, or that photos in high-temperature cells capture genetically controlled pigment loss.

## Broader FCP synthesis

The prospective white-versus-chromatic continuous phenotype-space W is supported as a repeated direction, but the independent hypothesis of preferential proportional geographic sorting of white photographs is not supported after component-baseline normalization (PR #249). An earlier three-cohort phenotypic IBD-like distance association is larger than the BIO5-associated IBE-like residual, without establishing that past distributions have larger causal effects. Separately, PR #248 finds a weak positive average BIO5 relation on a 1,574-species panel and poor spatial transfer by the tested species-specific climate models.

**Central manuscript question:** Why is a similar white–chromatic direction of floral variation repeatedly expressed across plant lineages while its geographic organization, climatic correlates and molecular mechanisms remain heterogeneous?

**Next mechanistic discrimination:** compare photographed morphs verified as inherited in native populations, neutral genotypes and population history, current accessible climate/soil/pollinator environments, and reproductive fitness in replicated environments. Do not claim that these historic and selective mechanisms have been estimated by the current photo-only analysis.

Scripts: scripts/analysis/fcp_nested_climate_one_vs_five_degree_20261020.py and scripts/analysis/fcp_withincell_month_climate_20261020.py.
