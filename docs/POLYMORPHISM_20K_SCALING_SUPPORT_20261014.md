# FCP single-paper 20,000-species support audit — 14 October 2026

**Purpose.** Clarify precisely why the global moisture–colour result uses approximately 20,000 photographed taxa while the replicate-validated within-species phenotypic-distribution and geographic-structure endpoints use fewer. This is a support/grain audit, not a new phenotype or successful biological replication. It does not change H1/H2, rerun the bootstrap or estimate genetic colour-morph frequencies.

## Outcome-blind original iNaturalist photographic capacity

Frozen canonical source: `main/docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md`. Observer-capped retained-photo capacity, measured **before inspecting colour outcomes**:

| Original source photo capacity per named taxon | Named taxa with at least this many photos |
|---|---:|
| >=1 | 41,322 |
| >=10 | 18,301 |
| >=20 | 12,985 |
| >=40 | 8,753 |
| >=60 | 6,770 |
| >=80 | 5,558 |
| >=100 | 4,730 |

All 42,111 photo-opportunity taxa, including ones without a qualifying current image, belong to an availability frame; **capacity is not successful colour classification**. The preceding numbers cannot be used as denominators for true global flower-colour polymorphism prevalence.

Consequences: for a fixed minimum **10 independent observer-capped photos per species**, 18,301 is the source-supported *upper* limit, even before checking botanical visibility, ambiguous flower organs, colour classification, overlap in individual localities and continuous-colour minority cluster. At **100 retained photos per species**, 4,730 is the source-supported opportunity limit. The existing high-depth series deliberately selected 500 discovery + 500 validation + 499 prospective species with 100 images each (149,900 images total), not all 4,730 eligible candidate species.

## Two dissimilar measurement grains, both belong in ONE New Phytologist Full Paper

1. **Between-species worldwide flower colour–environment directions (the paper's major macroecological result).** The sparse global photographic sample contains 100,543 *unique source photo IDs* from the 42,111 nominal taxa. Four-state photographic colour is classifiable for 45,953, including 45,836 with the 26 nonsoil climatic predictors from 25,028 taxa. The signed **genus-centred/geo-adjusted** worldwide BIO12-to-photo-white and mean vapour-pressure-to-photo-white analyses actually use **20,694 named species and 1,934 qualifying genera**, not all 25,028 climate-complete taxa; nominal genera with too few distinct sampled species have insufficient comparative support. Main quantitative effects: BIO12 +2.255 pp / original environment SD [scope-family maxT 95% +0.319,+4.190]; mean atmospheric vapour pressure +3.952 [+1.199,+6.704]. Climate26 45,836/25,028 and all36+soil 31,925/20,786 are *complete-case denominators*, not the genus-eligibility counts or the fit/test taxa of every predictive model. Held-out genus-equal Brier gains have 17,183 climate-only and 13,961 all36 test taxa.
2. **Within-species colour richness, two-mode achromatic–chromatic axis and local colour mixing.** These are *distributions*, not the colour of a representative taxon. The frozen high-depth resources obtained 149,900 photos from 1,499 intentionally disjoint sampled species (100 photo IDs per species). In the 499-species prospective cohort, 377 species satisfied the **>=40 four-state classifiable photos** measurement gate, and 158 entered the primary frozen continuous-phenotype W axis test; W=0.517 versus null median 0.457, P=.001. For the geographically distributed phenotype test at 50 km the source-qualifying sets were 166 discovery, 181 validation and 204 prospective taxa. Geographic inference further requires adequate pairs inside the radius and enough different sites, not just a high photo count.

**These populations are not equivalent.** The 100,543-photo macro sample has on average only about 2.4 recorded source photos per named taxon before classifiability; even on the climate-complete classifiable support there are 45,836/25,028 ~= 1.83 qualifying images per species on average. A photo-grain likelihood can model the dominant colour of such a species, but these averages do not permit a reliable **individual species' two-mode contrast direction** or **local colour mixture vs among-site differentiation**. The sparse and high-depth samples also arise from overlapping iNaturalist opportunity/sampling systems and should not be added as independent species confirmations.

## Why low within-species n creates genuine identifiability problems

If a putative minority hue has real local frequency `p`, the chance it appears in `n` idealized independent *classifiable* photos is `1-(1-p)^n`. For `p=0.05` this is only **40.1%** at n=10 and **99.4%** at n=100. At p=0.10 the chance at n=10 is only **65.1%**. Observers/flower organs/season and source-photo colour selection violate independence, so these illustrative binomial probabilities are **optimistic**, not a calibrated guarantee for current FCP photographs. A species observed only in one photographic colour is not necessarily monomorphic.

For an achromatic–chromatic W axis, the two continuous-colour clusters must be recoverable *within each species* before a unit phenotype direction exists. Shrinking a low-count species' composition toward a genus-wide mean is feasible as a Bayesian partial-pooling estimator, but **the resulting alignment can be imposed by the model**; this cannot be counted as independently confirmed within-species recurrence unless calibrated against disjoint rich-photo taxa and a matching construction-preserving null.

Likewise, two photographs in distant places do not establish local co-occurrence or geographic partitioning. The 50-km local comparison requires enough **observed local pairs** and a species-composition-preserving null. Similar photographed flowers at nearby coordinates are a photographic geographical pattern, not a census of individual genets or heritable morph ratios.

## Design for scaling the paper without diluting scientific claims

- **Full worldwide line (already near 20k):** retain ~20,694 eligible taxa for signed photographed colour–water associations, plus broad 26- and 36-variable models with their own complete-case/heldout denominators. Retain photo-classification missingness and original taxonomic cluster covariance.
- **Species-level uncertain variation:** a *separate* planned extension could analyse the <=18,301 names with >=10 observer-capped photos using count-aware multinomial/hierarchical likelihoods, report posterior intervals and probability of an undetected minority state, and label low support as `UNRESOLVED` rather than monomorphic. This is a **future analysis**, not an executed current 18,301-species result.
- **Confirmed two-mode phenotype direction:** keep the frozen source-based prospective W inference on its **158 actually evaluable species**. Extending W to more candidates needs new, observer-diverse continuous-colour photo measurements, explicit minor-cluster support and new frozen/heldout null designs. Neither eligibility (4,730 taxa at >=100 photos) nor total source image count equals successful W eligibility.
- **Local mixture versus geographic partitioning:** keep the three repeated-cohort qualified sets **166 / 181 / 204**, or require source-verified locality pair opportunities and matched nulls for any expanded analysis. Habitat and genotype/local population interpretations require additional data.

### Single-paper presentation rule

**One large-scale NP paper, two necessary inference resolutions.** High-depth 1,499 for valid distribution measurement/pattern and about 20,694 different-source-photo genus-eligible taxa for global environmental composition. Do **not** describe either as a failed attempt to cover exactly 20,000 species in *every* endpoint; the sampling effort needed for a distributional phenotype is qualitatively different from measuring a genus-level colour–environment correlation.

Source evidence is frozen in FCP original claim ledger (`main/docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`), capacity frame (`main/docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md`), manuscript-integrated source ledger (`docs/POLYMORPHISM_UNIFIED_EVIDENCE_AND_CLAIM_LEDGER_20261014.md`), and the two original photo/environment SHA-256 macro receipts under `results/fcp_macro_direction_20261013/`. No author's outcome selected threshold was changed by this note.
