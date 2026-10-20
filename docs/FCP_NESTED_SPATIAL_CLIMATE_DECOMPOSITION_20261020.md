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

## Broader FCP synthesis

The prospective white-versus-chromatic continuous phenotype-space W is supported as a repeated direction, but the independent hypothesis of preferential proportional geographic sorting of white photographs is not supported after component-baseline normalization (PR #249). An earlier three-cohort phenotypic IBD-like distance association is larger than the BIO5-associated IBE-like residual, without establishing that past distributions have larger causal effects. Separately, PR #248 finds a weak positive average BIO5 relation on a 1,574-species panel and poor spatial transfer by the tested species-specific climate models.

**Central manuscript question:** Why is a similar white–chromatic direction of floral variation repeatedly expressed across plant lineages while its geographic organization, climatic correlates and molecular mechanisms remain heterogeneous?

**Next mechanistic discrimination:** compare photographed morphs verified as inherited in native populations, neutral genotypes and population history, current accessible climate/soil/pollinator environments, and reproductive fitness in replicated environments. Do not claim that these historic and selective mechanisms have been estimated by the current photo-only analysis.

Scripts: scripts/analysis/fcp_nested_climate_one_vs_five_degree_20261020.py and scripts/analysis/fcp_withincell_month_climate_20261020.py.
