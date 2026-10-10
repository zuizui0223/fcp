# FCP original 85,337: scientific verdict after nearest-source species-pair full abiotic maxT null (2026-10-10)

**Decision: no robust individual thermal, precipitation, solar, wind, vapor-pressure or modelled soil predictor of original four-colour orientation among near paired localities within nominal species after photo-site geography/elevation residualization and multivariable source label-swap familywise control.**

## Source provenance and scope

- Original source 42,111 named plant taxa, separately 85,337 species×equal-area-cell photographed flower observations; **39,075** old four-state classifiable source photo cells and **46,262** unclassified kept outside the estimator. All original taxon ID + observation ID + photo ID exact matches to the real-site WorldClim BIO1–19, 12-month source solar/wind/vapor, actual site elevation and SoilGrids ten modelled edaphic properties (100,543 distinct original photo sites), no new pixels or flower-colour reclassification.
- Historical climate-complete classified 38,968 photographed species×cell rows from which exactly **5,791 nominal species** had >=2 source photographs across cells. Each of these 5,791 species was represented by ONE pair, the ORIGINAL two photographed sites with **smallest great-circle distance**, selected by geocoordinates and photo-ID tie-breaking only, before testing original hue categories. Soil-complete analysis keeps these same original pairs but excludes those missing soil; never selects a more favourable second-closest pair.
- Original four photographed hues are categorical white/yellow-orange/red-pink/blue-purple. The photo-level within-species vector difference is four one-hot components, **not** a numeric colour intensity or observed genetic anthocyanin phenotype.

## All geographic distance cutoffs and true source support

| Max distance between two photographed sites of the SAME nominal species | Climate-complete original one-pair-per-species count | Different original four-colour category | Original climate+soil complete pair count | Different colour with full soil |
|---|---:|---:|---:|---:|
| 50km | 96 | 22 | 51 | 13 |
| 100km | 293 | 71 | 163 | 37 |
| 250km | 805 | 184 | 438 | 102 |
| 500km | 1,725 | 415 | 945 | 220 |

50km analysis is **HOLD_INSUFFICIENT_INTRASPECIFIC_NEARBY_FOUR_COLOUR_PAIR_SUPPORT** under frozen >=30 genuinely colour-mismatched pairs and >=60 original source pairs; no significance p-value was computed there.

At 100/250/500km, for all fixed pairs, signed original physical environment differences were residualized against signed latitude/absolute-latitude/longitude sine/cosine/actual elevation differences using original source locations **without consulting photographed flower colour labels**. For all 25 climate-weather variables, or 35 climate+soil variables in the soil-complete subset, 999 within-pair sign flips exchange original photographed four-state colour **orientation**, retaining each species' exact two hues and physical sites. A shared flip acts on all four colour components and all environmental variables. Source conditional maxT corrects across every named variable within the radius and an additional Bonferroni ×8 conservatively spans all 4 overlapping radius ×2 complete-case selections. Full null results were saved; no winner-only source result selection.

## What actually passed the fixed null

**NONE of the 100, 250, or 500km climate-only / complete-soil support sets yielded any individual environmental variable with FWER-adjusted p below 0.05**, nor any environmental block exceeding the all-feature null; all original candidate feature values are present in the machine source receipt.

- 100km (293 original pairs /71 different hues) temperature, precipitation, solar, wind and vapor-pressure feature-block all-feature maxT p ranged ~0.998–1.0, with NO nominally positive single variable.
- 250km (805 original pairs /184 discordant): all climate feature blocks all-feature maxT p >=0.991. Soil-complete 438 pairs /102 colour mismatches showed individual SoilGrids available-water proxy nominal photo-colour orientation p=.028 but **within-radius all-variable FWER p=.505**, conservative across nested sets p=1.
- 500km climate (1725 original pairs /415 mismatches): annual water-vapor-pressure monthly climatology variable nominal p=.012, but **within-radius all-variable FWER p=.178**, conservative across nested sets p=1. Soil-complete 945/220 also had no FWER-supported block or single variable.

[Verified GitHub Actions source audit 38015951665](https://github.com/zuizui0223/fcp/actions/runs/38015951665), machine [result.json](result.json).

## How this changes the original manuscript-level claim

**Observational predictive fact retained:** original 85,337 species×cell photo data can weakly improve heldout four-state flower-colour prediction by adding environmental covariates to training-only species means and a low-rank original-site geodesic spatial covariance approximation: climate 20,546 heldout photo cells, +0.001861 Brier gain at 500km spatial kernel; full soil+climate 11,136 test photos, +0.003146. These conditional source scores are not proof of adaptation, source photo residual Moran remained nonzero, and the source pairwise null is a DIFFERENT estimand.

**New stronger same-species local test:** when exactly one outcome-blind nearest photographed site pair per species is used, and the photographed four-hue orientation is compared with signed environmental gradients remaining after geographic/elevation differences, **none of the full prespecified climate/solar/wind/vapor/soil variables have a permutation-calibrated familywise robust local directional association** at 100/250/500km. This is not a universal ecological zero: the choice of two old photographed sites per species, original ROI classifiability, incomplete true exposure, coarse hue categories, and local 50km power limits all restrict sensitivity.

**Cross-lineage scale:** directly observed LCVP dated-tree species tips cover only 342/872 and 649/1761 previously selected local original congeneric taxa, selected differently in climate/elevation; the remaining global 42,111 source taxa lack a reliable full-species dated phylogeny. Species fixed effects absorb invariant lineage components but do not simultaneously identify independent phylogenetic covariance. The 649 true-tree-tip phylogenetic+spatial-kernel test found no stable incremental improvement of full abiotic blocks across heldout genera/regions. Neither a photo-colour prediction gain nor a null can by itself demonstrate a heritable selective mechanism.

### Scientific gate

Do not claim a geography- and lineage-independent global thermal, precipitation, sunlight or soil driver of genetic flower-colour evolution on these source photos. The current warranted finding is **heterogeneous predictability of photographed flower-colour composition across observation/lineage/spatial scales and a quantitatively measured failure of fine geographic/phylogenetic identification of adaptive origins**.

The resulting p-values are conditional on within-species photographed hue exchangeability (not an exposure randomization), after linear signed geographic/elevation adjustment, using multiple outcome-exposed analyses and overlapping nested radii. Therefore this is a rigorous *sensitivity* for photographed four-colour orientation, but is NOT a generalized spatial phylogenetic GLMM or definitive falsification of any abiotic selective pressure.

No changes to old 1,499 photo-intensive H1/H2, main, or independent unexamined 2,000+730 source taxa.
