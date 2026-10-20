# Supplementary source-method record: previous long-form methods

This is a provenance-preserving copy of the original technical methods passages condensed in the hypothesis-led New Phytologist main manuscript (21 October 2026). Each item is copied from the PR #248 snapshot. These are not independent analyses or new findings.

## Original method 1

### Study design, data provenance and why each stage was required

The study uses two sequential high-depth image resources from a common outcome-blind iNaturalist sampling frame: **100,000 photographs from 1,000 species for discovery and species-disjoint validation**, followed by **49,900 newly sampled photographs from 499 species for prospective confirmation**. Later technical, climatic and phylogenetic data annotate these resources rather than forming additional flower-colour cohorts.

#### Stage 0 — global outcome-blind opportunity frame

Metadata-only discovery across an 18 × 9 equal-area grid identified **42,111 unique iNaturalist plant species** before candidate image pixels or flower-colour outcomes were opened; **4,730 species** could supply at least 100 retained photographs after observer capping. This outcome-blind frame prevents species assembly around the focal phenotype. It is a sampling opportunity frame, not a flower-colour dataset or prevalence denominator. Full filters and hashes are in `docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md`.

#### Stage 1 — discovery and species-disjoint validation image cohorts

The first resource contained **500 discovery + 500 validation species**, each with 100 photographs. Acquisition used Research Grade species-rank iNaturalist observations with photographs, georeferences, flowering annotation, positional accuracy <=5 km, open coordinates and permitted CC licences. Selection was colour-blind, each observer contributed at most two retained photographs per species, and deterministic geographic maximin sampling fixed 100 photographs per species. With the frozen >=40-classifiable rule, **369 discovery** and **363 validation** species entered D-based inference.

Repeated observations estimate within-species diversity, while the species-disjoint split separates discovery from validation on different species. No native-range restriction or explicit captive/wild filter was imposed; spatial estimands therefore refer to the observed community-photograph records.

#### Stage 1b — fresh-image transport of the D phenotype

A later same-system transport check used fresh photo IDs and the frozen four-state D definition; **136 species** were evaluable in both executions. Fresh images test transport beyond the original photo set, but not independent-source replication.

#### Stage 2 — prospective species-disjoint confirmation cohort

The white-versus-nonwhite axis was localized only after the original geometry was opened. Excluding all previously used high-depth species left an outcome-blind **3,230-species candidate frame**; deterministic hash selection froze 500 species and fresh retrieval yielded **499 species × 100 photographs = 49,900 new photographs**, with no replacement.

**A new species- and photo-disjoint cohort was therefore necessary** because a data-derived target cannot be prospectively confirmed in its discovery data. Before biological opening we froze the axis, W statistic, thresholds, structured null, support gates and one-shot execution rules. Because the cohort uses the same iNaturalist source and measurement system, it is prospective species-disjoint confirmation rather than independent-source replication.

#### Stage 3 — targeted technical and ecological annotations

After H2 terminalization, prospective photographs were reacquired for response-blind highlight metrics and annotated with **WorldClim 2.1 BIO5, BIO14 and solar radiation**; validation was linked to **V.PhyloMaker2**, and sampled span came from frozen coordinates. Highlight, climate and phylogeny test image exposure, a simple abiotic explanation and broad ancestry. **These are alternative-explanation and mechanism filters**; they cannot retroactively strengthen prospective H2.

#### Stage 4 — post-freeze unchanged-method fourth-cohort extension

After terminalization of the original three cohorts, the fourth outcome-blind draw added 498 previously unused taxa with 100 original photographs each, processed by the unchanged, blind ROI-v4/EfficientSAM pipeline (49,800 measured photographs). The unchanged H1 observer-disjoint splits, 50-km local-pair statistic (>=40 classifiable photos and >=30 local pairs; 199 fixed-composition vertex permutations), and the fixed white-axis W were applied. The intermediate original-method total was **1,997 species and 199,700 measured photographs**. This same-provider extension did not alter the original three-cohort frozen prospective H2 decision.

#### Stage 5 — completed authentic 447,500-photo extension

The fifth species-disjoint cohort contributed **2,478 additional species × 100 original photographs = 247,800 measured rows**. The complete five-cohort atlas therefore contains **4,475 species and 447,500 genuinely measured photographs**. Of these, **227,107** received a classifiable original four-state colour label. The original >=40-classifiable gate retained **3,297 species**; the additional >=30 conspecific photo pairs within 50 km retained **1,761 species**, consisting of 166 discovery, 181 validation, 204 prospective, 190 fourth-cohort and 1,020 fifth-cohort species. The original complete-composition 199-permutation 50-km null was replayed and additional pooled distributional tests and a white–chromatic versus chromatic–chromatic decomposition were performed *post hoc*. The 17 October statement that the fifth-cohort source was an uncompleted target was an interim evidence freeze, not the current status.

## Original method 2

### Complementary test: species-level D and within-species geographic organization

To distinguish structured geographic variation from unstructured colour noise or broader photographic coverage, **all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity**. Species-level organization was

`rho_i = Spearman(d_geo_ij, d_colour_ij)`.

Within each species, complete colour vectors were permuted among fixed coordinates, preserving geographic geometry and colour-vector composition; each species had **999 matched within-species null values**. Across species, D was correlated with observed `rho_i`; each frozen null realization was analysed identically, with upper-tail probability

`(1 + # {rho_null >= rho_obs}) / 1000`.

For the primary robustness analysis, **Rank(D) and rank(`rho_i`) are separately residualized** on ranked sampled span and clear technical-failure rate before correlation. Validation additionally used paired flower/background measurements as `Spearman(d_geo_ij, d_flower_ij - d_background_ij)`; **It is not the difference between separate flower and background Spearman coefficients**. Exact diversity-minimizing/maximizing ambiguity completions provide endpoint stress tests. These analyses test a structural correlate, not a causal maintenance mechanism.

## Original method 3

### Additional worldwide photographic-opportunity and geographic-sensitivity analyses

For the unchanged **447,500-photo original-method source** (identical archived CSV SHA256), we classified 227,107 photographs into the original four states and summarized eligible species using equal-taxon weights, without imputing unclassifiable photographs. We examined 18 longitude by nine equal-area latitude bands, retaining exact taxon-specific photographic and observer opportunities, and contrasted geographically sampled colour frequencies with taxon-conditional vertex label permutations that preserved each taxon's colour multiset and grid-cell sample sizes. A strict, *post hoc* flower-region purity screen (assigned coarse colour occupying at least 65% of the measured floral region, each competing colour no more than 20%) was applied to the same taxa and matched cells; it is a diagnostic for photographic colour composition, **not** an independent botanical identification of white morphs. Later sampling-month and year-month-stratified permutations fixed dates as well as the original geographic sampling arrangement. All geographic visualizations display the number of eligible taxa; empty cells are not biological absences.

A supplementary decomposition represented between-cell four-state total variation as the absolute difference in white frequency plus a residual component attributable to redistribution among nonwhite states. These components were compared both as absolute observed-minus-null differences and *relative to their own permutation-null expected opportunity*. The latter comparison is essential for testing any assertion that white contrasts are preferentially structured, because white contributes disproportionately to the raw number of mismatching photo pairs. The geographic-scale and technical quality screens were formulated after earlier patterns had been seen, and their p-values are exploratory. Scripts and frozen source receipts are documented in the 20 October world-pattern branch (Supporting Information).

### Prediction transport and nested thermal association diagnostics

In a separate exploratory assessment of the same selected 1,574-species climate panel, ten genus-held-out folds assessed whether shared BIO5 or BIO5–BIO14–radiation slopes discriminated within-species white/chromatic pairs in taxa from untrained genera. Five spatial folds held out complete 1° cells, with additional comparisons excluding conspecific training photos within 50 km, and evaluated species-specific penalized climatic models against intercept and geographic controls on identical successful held-out rows. Log-loss differences are defined as **reference loss minus augmented-model loss**; negative values indicate worse generalization. The genus-held-out shared-slope and spatially held-out species-specific comparisons are different estimands and must not be subtracted.

In a separately source-authenticated conditional-photo analysis, we contrasted the long-term normal maximum temperature of the warmest month (BIO5) at white-labelled versus nonwhite-labelled photographs within identical species, 5° geographic cells and calendar month, requiring multiple distinct observers and adequate photos of each label. With the same field and photographic colour-purity eligibility, we algebraically separated the average 5° contrast into variation **between** nested 1° cells and residual variation **within** those cells, with composition-conditioned label permutations. This decomposition does not reconstruct annual flowering-weather exposures or the history of genetic populations; smaller within-cell contrasts can also reflect reduced environmental or photographic variation.

## Original method 4

### H3a: broad phylogenetic signal

Shared ancestry could generate cross-species similarity without repeated ecological organization. We therefore tested validation-cohort D using Blomberg's K with 9,999 tip-label permutations under three frozen V.PhyloMaker2 placement scenarios (341 species each). Support required p<0.05 in all three raw-D scenarios; Pagel's lambda and opportunity-adjusted traits were sensitivities and could not rescue the primary test.

### H3b: sampled photographic span

Species photographed across wider areas may appear more diverse simply because more environments were sampled. We therefore tested whether the discovery association between D and log sampled span reproduced in the species-disjoint validation cohort using a 20,000-permutation Spearman test, with observer/classifiability-adjusted partial-rank and rank-PGLS sensitivities. Sampled photographic span is not interpreted as true biological range size.

