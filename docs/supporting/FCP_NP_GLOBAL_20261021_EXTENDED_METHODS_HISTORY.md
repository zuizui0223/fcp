# New Phytologist global synthesis — expanded source-method contracts (historical stage copy)

This manuscript-preserving supplement records technical method paragraphs removed from the shorter main paper on 21 October 2026. These passages originated from the versioned PR #248 manuscript and do not constitute independent results. The authoritative measured source and frozen tests are unchanged.

## Archival original subsection 1

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

## Archival original subsection 2

### Complementary test: species-level D and within-species geographic organization

To distinguish structured geographic variation from unstructured colour noise or broader photographic coverage, **all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity**. Species-level organization was

`rho_i = Spearman(d_geo_ij, d_colour_ij)`.

Within each species, complete colour vectors were permuted among fixed coordinates, preserving geographic geometry and colour-vector composition; each species had **999 matched within-species null values**. Across species, D was correlated with observed `rho_i`; each frozen null realization was analysed identically, with upper-tail probability

`(1 + # {rho_null >= rho_obs}) / 1000`.

For the primary robustness analysis, **Rank(D) and rank(`rho_i`) are separately residualized** on ranked sampled span and clear technical-failure rate before correlation. Validation additionally used paired flower/background measurements as `Spearman(d_geo_ij, d_flower_ij - d_background_ij)`; **It is not the difference between separate flower and background Spearman coefficients**. Exact diversity-minimizing/maximizing ambiguity completions provide endpoint stress tests. These analyses test a structural correlate, not a causal maintenance mechanism.

## Archival original subsection 3

### H3a: broad phylogenetic signal

Shared ancestry could generate cross-species similarity without repeated ecological organization. We therefore tested validation-cohort D using Blomberg's K with 9,999 tip-label permutations under three frozen V.PhyloMaker2 placement scenarios (341 species each). Support required p<0.05 in all three raw-D scenarios; Pagel's lambda and opportunity-adjusted traits were sensitivities and could not rescue the primary test.

### H3b: sampled photographic span

Species photographed across wider areas may appear more diverse simply because more environments were sampled. We therefore tested whether the discovery association between D and log sampled span reproduced in the species-disjoint validation cohort using a 20,000-permutation Spearman test, with observer/classifiability-adjusted partial-rank and rank-PGLS sensitivities. Sampled photographic span is not interpreted as true biological range size.

## Archival original subsection 4

### Additional worldwide photographic-opportunity and geographic-sensitivity analyses

For the unchanged **447,500-photo original-method source** (identical archived CSV SHA256), we classified 227,107 photographs into the original four states and summarized eligible species using equal-taxon weights, without imputing unclassifiable photographs. We examined 18 longitude by nine equal-area latitude bands, retaining exact taxon-specific photographic and observer opportunities, and contrasted geographically sampled colour frequencies with taxon-conditional vertex label permutations that preserved each taxon's colour multiset and grid-cell sample sizes. A strict, *post hoc* flower-region purity screen (assigned coarse colour occupying at least 65% of the measured floral region, each competing colour no more than 20%) was applied to the same taxa and matched cells; it is a diagnostic for photographic colour composition, **not** an independent botanical identification of white morphs. Later sampling-month and year-month-stratified permutations fixed dates as well as the original geographic sampling arrangement. All geographic visualizations display the number of eligible taxa; empty cells are not biological absences.

A supplementary decomposition represented between-cell four-state total variation as the absolute difference in white frequency plus a residual component attributable to redistribution among nonwhite states. These components were compared both as absolute observed-minus-null differences and *relative to their own permutation-null expected opportunity*. The latter comparison is essential for testing any assertion that white contrasts are preferentially structured, because white contributes disproportionately to the raw number of mismatching photo pairs. The geographic-scale and technical quality screens were formulated after earlier patterns had been seen, and their p-values are exploratory. Scripts and frozen source receipts are documented in the 20 October world-pattern branch (Supporting Information).

### Prediction transport and nested thermal association diagnostics

In a separate exploratory assessment of the same selected 1,574-species climate panel, ten genus-held-out folds assessed whether shared BIO5 or BIO5–BIO14–radiation slopes discriminated within-species white/chromatic pairs in taxa from untrained genera. Five spatial folds held out complete 1° cells, with additional comparisons excluding conspecific training photos within 50 km, and evaluated species-specific penalized climatic models against intercept and geographic controls on identical successful held-out rows. Log-loss differences are defined as **reference loss minus augmented-model loss**; negative values indicate worse generalization. The genus-held-out shared-slope and spatially held-out species-specific comparisons are different estimands and must not be subtracted.

In a separately source-authenticated conditional-photo analysis, we contrasted the long-term normal maximum temperature of the warmest month (BIO5) at white-labelled versus nonwhite-labelled photographs within identical species, 5° geographic cells and calendar month, requiring multiple distinct observers and adequate photos of each label. With the same field and photographic colour-purity eligibility, we algebraically separated the average 5° contrast into variation **between** nested 1° cells and residual variation **within** those cells, with composition-conditioned label permutations. This decomposition does not reconstruct annual flowering-weather exposures or the history of genetic populations; smaller within-cell contrasts can also reflect reduced environmental or photographic variation.

### Additional checks of broad ancestry and photographic span

The validation sample's original phylogenetic signal diagnostic for species-level D used the three pre-frozen tree placements (Blomberg's K), and the previously frozen sampled-span regression was compared across discovery and validation rather than optimized after seeing validation. These **secondary original-stage diagnostics** do not estimate neutral within-species ancestry, migration or range expansion. Full placement variants and numerical checks are retained in Supporting Information.

### Reproducibility and frozen decisions

Frozen protocols, machine-readable results, input hashes and no-rescue rules govern interpretation. Prospective H2 used one authorized execution with no replacement, threshold changes, axis refitting or rerun-based selection after biological opening. Exact recovery routes are given in the data-lineage map and Supporting Information.

---

## Table 1. Data sources, lineage and inferential necessity

| Stage / data resource | Source and selection | Scale used in the paper | Why it was required | Inferential role |
|---|---|---:|---|---|
| Global opportunity frame | Metadata-only iNaturalist discovery across a fixed 18 × 9 equal-area grid; no flower-colour outcomes used | 42,111 species; 4,730 with capacity for >=100 retained photos | Define the candidate universe independently of flower colour and quantify high-depth sampling opportunity | Sampling frame only; not a polymorphism-prevalence denominator |
| Discovery image cohort | Colour-blind high-depth iNaturalist sampling under fixed filters, observer cap and geographic maximin design | 500 species × 100 photos; 369 D-eligible species | Estimate within-species distributions at high depth and discover candidate recurrent colour-space structure | Discovery/calibration |
| Species-disjoint validation image cohort | Same acquisition contract, different 500 species | 500 species × 100 photos; 363 D-eligible species; 341 in phylogenetic analysis | Test reliability and spatial associations on species not used for discovery | Validation/replication within the same source system |
| Fresh-image D transport | Newly sampled iNaturalist photo IDs measured under the same frozen four-state system | 136 species evaluable in both executions | Test whether D transports beyond the original observer partitions to a fresh photograph set and separate execution | Same-source measurement transport |
| Prospective confirmation cohort | All previously used high-depth species excluded; 500 species frozen outcome-blind from 3,230 candidates; fresh photo IDs, no replacement | 499 species × 100 photos = 49,900; 377 measurement-evaluable; 158 in primary H2 | The white axis was identified after Stage 1 was opened, so a new untouched species/photo cohort was required for prospective confirmation | Primary confirmatory H2 test |
| Highlight-validity reacquisition | Reacquisition of the same 49,900 prospective rows after H2 terminalization; response-blind exposure metrics frozen before joining white outcomes | 49,900 reacquisitions; highlight metrics for 44,098 rows | Quantify whether white classification is coupled to image exposure, which the H2 biological null cannot establish | Post-confirmatory measurement-validity bound |
| Climate annotation | WorldClim 2.1 BIO5, BIO14 and solar radiation assigned to prospective-cohort coordinates; post hoc observer-balanced/paired sensitivities; fixed BIO5 transport to discovery/validation cohorts | Primary panel: 281 species / 12,583 rows; observer-balanced: 352 species; observer-paired: 106 species / 144 paired species-observer strata | Test a simple abiotic explanation, its observer sensitivity and whether it transports across species-disjoint cohorts | Secondary mechanism/alternative-explanation test |
| Phylogenetic annotation | V.PhyloMaker2 S1–S3 placements joined to validation-cohort D | 341 validation species | Test whether broad shared ancestry accounts for between-species D differences | Alternative-explanation filter |

| Post-freeze fifth-cohort and complete source photographs | Authentic ROI-v4/EfficientSAM measurements from all five species-disjoint iNaturalist cohorts, with original fixed 50-km matched vertex null | 4,475 species / 447,500 genuine photographs; 227,107 classifiable; 1,761 spatially eligible species | Extend geographic mismatch organization beyond the four-cohort 199,700-photo interim freeze | Later post hoc five-cohort spatial and composition inference; not a new frozen H2 confirmation |
| Expanded same-species environmental panel | WorldClim 2.1 location normals, joint within-species Firth regression and two-stage normal-REML synthesis | 1,574 species / 86,114 photographs / 727 genera | Identify average photographed white–climate association in the complete atlas, keeping species identity and separating colour from climatic composition | Exploratory environmental association; no causal selection estimate or date-specific weather |

The central sampling logic is therefore **frame → discovery/validation → frozen target → prospective confirmation**. The later highlight, climate and phylogenetic analyses answer different validity or explanatory questions and do not create additional independent flower-colour cohorts. The prospective confirmation cohort remains within the same iNaturalist opportunity universe and measurement system, so it is species- and photo-disjoint confirmation rather than independent-source replication.


---

## Results

### Worldwide photographic variation: a common repertoire, not a global genetic-polymorphism census

The outcome-blind five-cohort source produced **447,500 actual photograph measurements of 4,475 species**; **227,107 photographs** were assignable to the fixed four-state classifier. Among **3,297 species** with at least 40 classifiable photographs and three observers, **3,029** exhibited two or more photographic colour labels and **2,963** included both white and nonwhite labels. These are evidence of labelled photograph heterogeneity, not estimates of the prevalence of genetically inherited floral polymorphisms in world angiosperms. A high-purity flower-colour filter retained 162,003 classified photographs. In **the same 1,813 eligible species**, mean finite-sample-corrected colour-pair diversity decreased from **0.11544** to **0.06468**, exposing sensitivity to photograph inclusion even without a change in taxon composition.

Coverage limited geographical inference. Of 162 equal-area worldwide grid cells, 96 contained at least one classifiable source photograph, 65 supported at least one eligible species-cell comparison, and only **28** contained at least 20 eligible taxa. The apparent higher photo-colour diversity at lower latitudes was not robust as a same-taxon difference: among **246 species** observed in both tropical and extratropical bands, mean paired tropical-minus-extratropical diversity was **+0.0121**, species-bootstrap 95% interval **[−0.0058, +0.0298]**. The descriptive atlas therefore establishes an observational phenotype distribution and its sampling gaps, **not** a biological global latitudinal polymorphism gradient.

### H1: species-level polymorphism is reproducible across observer-disjoint photo sets

The first-frozen repeated-partition H1 test supported observer-disjoint reproducibility in the validation cohort. Across 200 partitions, the median number of paired species was 329. Median split Spearman rho was **0.7891**, with a 5th percentile of **0.7652** and a 95th percentile of 0.8109. Median Lin CCC was **0.8548**, and median Spearman-Brown projected reliability was **0.8821**.

The later deterministic stress test retained 363 validation species with zero observer leakage and yielded rho = **0.7927** (bootstrap 95% interval **0.7419–0.8324**) and CCC = **0.8474**. This split missed its deliberately stricter prespecified rho = 0.80 floor by 0.0073. We therefore retain the first-frozen H1 support while explicitly rejecting a claim of near-perfect or split-invariant reliability.

The later fresh-image transport check provided a stronger same-system replication of D. Among **136** species evaluable in both the earlier frozen analysis and the fresh execution, D showed Spearman rho = **0.9681**, Lin CCC = **0.9720**, and a calibration slope of **0.9691**. The median absolute change in D was **0.0148**, while mean signed fresh-minus-prior change was **+0.0032**. Thus the between-species ordering and scale of D transported strongly across a fresh photo set and a separate measurement execution within the same iNaturalist/FCP measurement system.

These results admit D as a reproducible high-depth species phenotype for the subsequent geometry analyses, but they do not estimate global polymorphism prevalence.

### Finite-sample correction does not alter D-based conclusions

Finite-sample correction had negligible influence on the species-level phenotype. Raw and corrected D rankings were almost identical in discovery (Spearman rho = **0.999977**) and validation (**0.999978**); mean absolute changes were **0.00486** and **0.00469**, respectively. The primary D–spatial result was unchanged after correction: discovery partial rho changed from 0.126637 to **0.126256** with p = **0.007**, and validation from 0.099288 to **0.099561** with p = **0.025**. The validation flower-minus-background sensitivity likewise remained supported (corrected partial rho = **0.116155**, p = **0.010**).

The sampled-span conclusion was also unchanged (validation raw rho = -0.0026, p = 0.9586; corrected rho = **-0.0024**, p = **0.9629**), and corrected validation Blomberg-K effect sizes remained nearly identical to raw values across S1-S3. Thus the main D-based conclusions are not explained by finite-sample bias in the plug-in diversity estimator.

### H2 discovery and audit: recurrent geometry localizes to white versus nonwhite

The original H2 analysis first established recurrent label-free geometry before naming its biological direction. At the primary 0.10 tier, discovery had 152 vector species with leading-axis concentration lambda1 = **0.541412**; validation had 129 vector species with lambda1 = **0.535045**. Validation mean squared projection onto the frozen discovery axis was **0.524157**, and the independently fitted discovery and validation leading axes had absolute alignment of approximately **0.986**. Discovery concentration and validation frozen-axis transport each exceeded both the isotropic reference and the subsequent coarse-state-preserving structured null (structured-null p = **0.001** in both cases). The strict 0.20 tier showed the same direction of support.

Audit then localized that recurrent component to white versus nonwhite. In the original discovery and validation cohorts, the fixed white-axis statistic exceeded the construction-preserving structured null at both admissibility tiers.

At the primary 0.10 tier, discovery contained 152 eligible species with W = **0.514625** (null median 0.430808, p = **0.001**) and validation contained 129 species with W = **0.514586** (null median 0.466546, p = **0.001**).

At the strict 0.20 tier, discovery contained 75 species with W = **0.542355** (null median 0.443626, p = **0.001**) and validation contained 65 species with W = **0.510517** (null median 0.469943, p = **0.008**).

When the white contrast was projected out, residual directional concentration no longer exceeded the same construction-preserving null. A non-white-only diagnostic likewise did not support a recurrent chromatic hue axis. Thus the original-cohort H2 signal was narrow: the repeatable component was concentrated along an achromatic–chromatic white-versus-nonwhite direction rather than a general shared hue direction.

Because q_white was identified after the original broad H2 geometry had been opened, these results establish target localization but not untouched prospective confirmation.

### Prospective confirmation measurement gate

The prospective-confirmation run completed all **256 / 256** terminal partitions and produced **49,900** terminal rows from **499** species. All measurement IDs were unique; no species or rows were replaced; image pixels were not persisted.

Of the 49,900 rows, **25,788** were classifiable into the four biological states and 24,112 were nonclassifiable. Using the prespecified minimum of 40 classifiable rows per species, **377** species were measurement-evaluable, exceeding the required minimum of 250. The measurement-support gate therefore passed before H2 was opened.

### Prospective H2: the frozen white axis is confirmed

At the primary 0.10 tier, **158** species produced admissible nonzero displacement vectors. Observed W was **0.517**. Across 999 structured-null worlds, the null median was **0.457** and the 95% interval was **0.436–0.475**. Observed W was 1.13 times the null median, with upper-tail p = **0.001**.

At the strict 0.20 tier, **86** species were eligible. Observed W was **0.533**, compared with a null median of **0.459** and a 95% interval of **0.433–0.487**. Observed W was 1.16 times the null median, again with p = **0.001**.

The frozen terminal verdict was therefore

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

This is an untouched prospective test of a previously frozen axis in a species-disjoint cohort. Because the prospective confirmation cohort was drawn from the same iNaturalist source/opportunity universe and processed with the same measurement system, it is not described as an independent-source replication.

### What the prospective H2 contrast confirms, and post hoc validity diagnostics

The structured-null baseline was itself strongly white-axis aligned: the isotropic eight-dimensional expectation is 0.125, whereas the frozen primary structured-null median was **0.457** and observed W was **0.517**. The confirmatory result therefore concerns the **increment above a coarse-state-preserving construction baseline**, not the entire difference from isotropy. Of the 158 primary vector species, **137 (86.7%)** had white as one of their two leading coarse morphs; their mean per-species white-axis contribution was 0.562, compared with 0.224 among the remaining 21 species.

A post-confirmatory null diagnostic re-applied the continuous minor-cluster gate within each of 299 structured-null worlds, starting from the 185 species that passed the coarse-state gate. The null median was **0.447** (95% interval **0.429–0.470**; plus-one p = **0.0033**), lower than the frozen conditional-null median of 0.457. Thus the preregistered null was slightly more conservative with respect to this gate asymmetry.

The direct one-shot digital-highlight control completed all 49,900 frozen reacquisitions with zero acquisition failures and zero source-byte drift. Highlight metrics were available for 44,098 rows. The response-blind q95 rule set the high-clip threshold at near_clip_fraction > **0.4933**, identifying **2,205** rows. In the within-species conditional logistic model, 23,320 originally classifiable rows from 461 species contributed to inference. A one-SD increase in within-species near_clip_fraction increased the odds of frozen white classification by **1.444** (95% CI **1.389–1.502**). The entire interval lies above the predeclared negligible-coupling upper bound of 1.25, providing direct evidence that the measured white state is exposure-coupled.

Removing the response-blind high-clip set did not remove H2 support. The sensitivity retained **142** primary-tier vectors, with W = **0.503**, structured-null median = **0.454**, and upper-tail p = **0.001**. This corresponds to **89.9%** retention of the original 158 vectors, one vector below the prespecified >=90% requirement. The frozen executable checked this retention criterion before its coupling-CI branch and therefore returned the terminal machine state **INDETERMINATE**. The separately frozen prose protocol would classify the coupling interval itself as a FLAGGED condition because its complete 95% CI exceeds 1.25. Because both the prose rule and executable precedence were frozen before reacquisition, we do not recode the terminal machine state after observing the result. Instead, we report the two facts directly: near-clipping is substantially coupled to white classification, while excess H2 alignment remains supported after removal of the response-blind high-clip set.

### The prospective BIO5 association is observer-sensitive and does not transport as a common rule

The post-confirmatory environmental-filter test retained **281** prospective-cohort species and 12,583 rows after the response-blind high-clip exclusion. BIO5 was the only one of the three frozen environmental predictions to pass its full gate. The median within-species white-minus-nonwhite BIO5 contrast was **+0.0690 SD**; 57.3% of species had a positive contrast. The species-level Wilcoxon p-value was 0.0118 and remained significant after Holm correction across BIO5, BIO14 and solar radiation (**Holm-adjusted p = 0.0354**). In the species-stratified conditional logistic model, a one-SD within-species increase in BIO5 was associated with OR = **1.073** for white classification (95% CI **1.029–1.119**, p = **0.000919**) after continuous near-clip adjustment. BIO14 and mean solar radiation did not pass their frozen gates (Holm-adjusted p = 0.692 for each).

Post hoc observer controls weakened that BIO5 signal. In the strict observer-paired design, **106 species** contributed 144 species-observer strata containing both colour states; the median BIO5 contrast was **0.000 SD** (Wilcoxon p = **0.485**) and the conditional estimate reversed direction without statistical support (OR = **0.787**, 95% CI **0.544–1.138**, p = **0.202**). In the broader observer-balanced design, **352 species** retained a positive median contrast of **+0.0541 SD**, but the Wilcoxon test was not significant (p = **0.0750**; sign test p = **0.0484**). These sensitivities were post hoc and do not overwrite the frozen within-cohort gate, but they show that the magnitude and support of the BIO5 association depend on observer conditioning.

The fixed species-disjoint BIO5 transport test did not reproduce that association across both original high-depth cohorts. In discovery, **271** eligible species had a median contrast of **-0.0068 SD** (Wilcoxon p = **0.743**), and the species-stratified estimate was OR = **1.033** (95% CI 0.990–1.078, p = **0.131**). In validation, **260** species had a directionally concordant median contrast of **+0.0662 SD**, but the species-level test missed the frozen criterion (p = **0.0541**); the row-level model was positive (OR = **1.048**, 95% CI 1.004–1.094, p = **0.0319**). Because both cohorts were required to pass at both inferential levels, the prespecified transport criterion was not met.

Thus the frozen primary prospective analysis provides a within-cohort association between white states and warmer maximum-temperature environments, but post hoc observer controls weaken that association and the fixed transport test does not support a common cross-cohort BIO5 rule (Fig. S9a).

### Greater D is associated with stronger within-species geographic colour organization

The positive D–spatial association reproduced across the two species-disjoint high-depth cohorts. The raw association was rho = **0.0892133** (p = **0.034**) in discovery and rho = **0.1016008** (p = **0.025**) in validation.

After controlling for sampled geographic span and clear ROI/flip technical-failure rate, the geometry-preserving analysis remained positive in discovery (partial rho = **0.1266367**, p = **0.007**) and validation (partial rho = **0.0992877**, p = **0.025**). In validation, the matched flower-minus-background response was also positive (partial rho = **0.1162411**, p = **0.010**).

The validation result remained supported under exact uniform ambiguity-endpoint completions: primary D_min4 rho = **0.0970781** (p = **0.029**) and D_max4 rho = **0.1252858** (p = **0.008**); flower-minus-background p-values were **0.009** and **0.006**. Species with greater species-wide sampled colour-state diversity therefore tend to show stronger internal geographic organization, although this association neither partitions D into within- versus among-location components nor identifies the causal process that creates or maintains that organization.

### Species-wide flower-colour diversity is geographically partitioned

The distance-decay result does not by itself say whether range-wide diversity reflects local coexistence or differentiation among localities. A post-outcome composition-preserving test addressed this directly. At the fixed 50-km primary scale, nearby observations contained less four-state colour diversity than expected from each species' exact overall colour composition in all three species-disjoint cohorts. Mean species-wide minus local pairwise diversity was **0.02053** in discovery (n = 166), **0.01867** in validation (n = 181) and **0.01468** in the third cohort (n = 204); the matched vertex-null upper-tail p-value was **0.005** in each cohort. Local depletion was also supported at 25, 100 and 250 km in all three cohorts.

The result survived five direct falsification tests at 50 km. After excluding all same-observer photo pairs, mean depletion remained **0.01998**, **0.01690** and **0.01457** in discovery, validation and the third cohort, respectively (p = **0.005** in each). Removing every white observation and analysing only yellow/orange, red/pink and blue/purple retained positive depletion in all three cohorts (0.01947, p = 0.005; 0.00729, p = 0.025; 0.00891, p = 0.040). Replacing coarse states with continuous nine-colour Jensen–Shannon dissimilarity also retained local homogeneity (mean species-wide minus local JSD = **0.01518**, **0.01366** and **0.01301**; p = **0.005** in each cohort). Calendar-quarter-stratified label permutations likewise retained local depletion in discovery, validation and the third cohort (**0.02287**, **0.02061**, **0.01560**; p = **0.005** in each), as did restricting local comparisons to photographs from different years (**0.02095**, **0.01437**, **0.01377**; p = **0.005** in each). Thus the geographic partitioning is not restricted to the exposure-coupled white classifier, coarse morph assignment, same-observer photo pairs, calendar-quarter composition or same-year sampling.

A stronger claim that higher-D species allocate a larger fraction of their diversity among localities was not scale invariant. At 50 km the D–depletion correlation was unsupported in discovery (rho = 0.0390, p = 0.135) but positive in validation (rho = 0.1171, p = 0.025) and the third cohort (rho = 0.1118, p = 0.005). We therefore treat **local depletion itself**, rather than a universal D–partitioning slope, as the replicated result.

### Post-freeze fourth-cohort replication under the unchanged estimators

After the three-cohort results were frozen, **498 additional species and 49,800 genuinely measured photographs** were processed with the original image and statistical methods. At **50 km (post hoc)**, **190** new species met both unchanged classifiability and local-pair gates, with species-equal mean depletion **0.0172936** (199-permutation upper-tail **p = 0.005**). Across all four measured cohorts, **741** eligible species gave pooled depletion **0.0176372** (p = **0.005**) from **199,700** photographic measurements. The original discovery, validation and prospective-cohort effects reproduced without changes to their numerical values.

In the new fourth cohort, D reproducibility across 200 observer-disjoint splits had median Spearman rho **0.8056** (median 316 paired species; median concordance correlation **0.8441**). The fixed white-axis test gave W **0.5368** versus structured-null median **0.4363** (primary 0.10 gate, n = **134**, p = **0.001**); the strict 0.20 gate gave W **0.5473** versus **0.4581** (n = **67**, p = **0.001**). This is a *post-freeze, taxon-disjoint, same-iNaturalist-source extension*, not independent-source replication or a new untouched prospective H2 confirmation. W measures squared axis alignment, not a directional rise in white morphs. The completed fifth-cohort and full 447,500-photo results are reported below. The four-cohort results are preserved as their chronologically earlier analysis snapshot and are not re-described as the frozen prospective H2 test.

### Post hoc robust spatial effect distribution under matched nulls

A further species-level audit asked whether the positive *cohort mean* was driven solely by unusually large effects. Of the **741** eligible species across four cohorts, **417** (56.3%) had positive individual depletion, **61** zero, and **263** negative. The species median was **0.005974**, below the species-equal mean **0.017637**. After excluding the highest-effect **149 species** (20%), the raw mean was **−0.006553**. This sign change by itself does not refute spatial sorting: when the **same top-20% exclusion was performed anew in each of the 199 original species-conditional permutation worlds**, the null mean was **−0.018278**, giving exploratory upper p = **0.005** (Holm p = **0.020** within four newly examined robust statistics). The median, fraction positive (observed **0.56275** vs null mean **0.45007**), and top-10% trimmed mean (**0.001911** vs null **−0.011139**) also exceeded their matched nulls (all pooled within-four Holm p = **0.020**). Across individual cohorts, these Holm p-values ranged **0.020–0.030** (Supporting Information S5e). This analysis remains post hoc and does not control all prior exploratory tests.

### Geographic sorting across all 447,500 authentically measured photographs

The final fifth cohort added **2,478 new species and 247,800 original measured photographs**, bringing the five-cohort atlas to **4,475 species, 447,500 photographs**. Of 227,107 classifiable photographs, the unchanged >=40-usable rule retained 3,297 species; **1,761 species** additionally met the original >=30 within-50-km pair criterion. In this species-equal panel, mean global four-state photo-pair mismatch was **0.224305** compared with **0.206690** within 50 km: a global-minus-local depletion of **0.0176149** (matched label-permutation null mean 0.0000446; 199-permutation upper-tail p = **0.005**). The fifth cohort separately showed mean depletion **0.017599** in 1,020 eligible species, consistent in sign and scale with the first four cohorts. These are photographic differences rather than estimates of genetically segregating morph frequencies.

The full pooled species median was **0.006871**, and **984 of 1,761 species (55.88%)** had positive depletion, versus a matched-null mean of **44.58%**. Both statistics exceeded the matched permutation null (Holm p = **0.020** within the four *post hoc* pooled distributional checks). Removing the highest-effect 20% in both the observations and every resampled null left an **observed mean of −0.006247**, still greater than the equally trimmed null mean **−0.018278** (within-family Holm p = **0.020**). Accordingly, the aggregate signal is not confined to the highest-effect species, but neither do these results show every species to be spatially differentiated. The separate stricter 5-km, same-site/month and independent-observer analysis did *not* support local within-population polymorphism; the broader 50-km result does not reverse that inference (Fig. S12).

In a further, exploratory equal-area-grid audit of **1,281 taxa matched across ordinary and strict-purity photographs**, conditioning the label permutations on the **same year and month** reduced the strict-purity regional total-variation excess from **0.007263** (unrestricted shuffle) to **0.003515** (year–month conditioning), which remained above the source-conditional permutation null (p = **0.005**). This small effect was not uniform across source cohorts: only two of five cohorts individually supported a positive mean under the strictest condition, and omission of the largest new-U100 cohort gave p = **0.075**. Temporal/photo-sampling conditioning therefore attenuates observed organization; it does not supply evidence of genetically persistent local colour morphs. In the separately frozen stricter 5-km analysis the matched same-month statistic was **not supported** (n = 194, p = **0.300**).

### The white–chromatic axis contributes most absolute sorting but is not disproportionately strong

Exact partitioning of the full 1,761-species effect yields **0.0147632** from unequal white–chromatic photograph pairs and **0.00285166** from chromatic–chromatic pairs. Their sum equals **0.0176149**; white–chromatic contrasts therefore constitute **83.81%** of the aggregate geographic mismatch depletion. However, they already account for **79.57%** of global unequal-pair mismatch, and *both* components individually exceed their matched within-species label permutations (exploratory two-component Holm p = **0.010** for each).

The abundance of white-labelled differences alone does not establish a special ecological advantage of a white phenotype. Among **636 species** with each type of global unequal-pair mismatch above 0.01, the respective normalized null-adjusted depletion averages were **0.05696** (white–chromatic) and **0.06945** (chromatic–chromatic); the paired difference was **−0.01249**, with genus-cluster 95% interval **[−0.06813, +0.04591]**. This difference was not resolved, and the direction was sensitive to the exploratory baseline threshold (Fig. S13). The recurrent continuous white–chromatic displacement direction established by the original prospectively tested W statistic is a different, independently defined phenotype-space property; neither W nor this absolute decomposition estimates pigment biosynthetic costs or natural selection.

At a second, exploratory **interregional scale**, the matched 1,282-taxon strict-purity subset likewise yielded a large white-containing proportion of the *absolute* excess colour turnover (**89.5%**). However, when both white and nonwhite components were normalized by their **own** permutation-null expected regional turnover, the white-minus-nonwhite difference was **+0.0934** for 174 sufficiently informative taxa (95% species-bootstrap interval **[−0.1160, +0.3186]**, matched p = **0.220**). The corresponding broad-quality comparison also failed (n = 395, p = 0.370). Hence recurrent phenotype-space alignment with the achromatic–chromatic direction **does not imply preferential geographic sorting of white-labelled photo variants**.

### A small positive BIO5–white association in the expanded within-species dataset

The distinct **post hoc** five-cohort analysis of **1,574 species and 86,114 photographs** found a small positive average association between white classification and BIO5 (maximum temperature of the warmest month), after simultaneous adjustment for BIO14 and mean solar radiation within each species. The equal-species mean log-odds slope in common units was **+0.03233 per 1°C**, with 5,000-draw genus-bootstrap 95% interval **[+0.01786, +0.04728]**. The separate two-stage REML synthesis gave **+0.01708 per 1°C** (modified-HKSJ 95% interval **[+0.00897, +0.02518]**); on each species' sampled-climate-SD scale its mean was **+0.05467** [0.03149, 0.07785]. Each number refers to a distinct summary estimand, not a percentage-point difference in white-flower prevalence.

BIO14 and solar radiation had mean intervals spanning zero. The small BIO5 direction persisted in specified observer-thinned, spatially thinned and geographically adjusted sensitivity fits, but these do not fully eliminate spatial, seasonal, cultivation, lighting or observer biases. In particular, this larger-sample *exploratory* mean association does **not** retrospectively pass the earlier *frozen* two-cohort BIO5 transport gate, and it does not demonstrate white-morph fitness advantage in warm environments (Fig. S14).

### Weak global climatic direction does not predict novel regions reliably

Shared BIO5 coefficients transported with a **very small** positive discrimination gain in ten genus-held-out folds (mean conditional pair log-loss gain **+0.000699 nats**, AUC **0.5161**), whereas a three-variable common-climate model offered no clear additional benefit. By contrast, when **complete 1° cells** were held out, species-specific climate models performed worse than within-species intercept baselines (mean log-loss gain **−0.0599**, 1,570 species) and models augmenting geographic predictors with climate performed worse than the geographic baselines (**−0.1498**, 1,237 species). Both losses worsened with an additional 50-km conspecific exclusion. Genus-held-out common slopes and geographically held-out species-specific fits address different prediction tasks, but jointly show why a small average climate association should not be converted into a universal spatially transferable rule.

Matching species, 5° cells and photographic calendar month gave a positive BIO5 difference between white- and nonwhite-labelled photo localities in **327 strict-purity species** (mean **+0.32985°C**, 199-permutation p = **0.005**). An exact, nested spatial decomposition allocated **+0.27158°C** of the signed contrast to *between* 1° geographic cells and **+0.05827°C** to the *within*-1° component; the latter did not pass the stricter within-subcell permutation test (p = **0.075**). This is a decomposition of photographed locality associations, **not** a partition of causal selection, dispersal history or pigment physiology.



### Geographic colour turnover contains both IBD-like and bounded IBE-like components

Using continuous nine-colour dissimilarity, BIO5 difference explained a small positive component of colour turnover after geographic distance was rank-partialled out in all three cohorts: mean phenotypic IBE-like rho was **0.00925** in discovery (p = **0.005**), **0.00826** in validation (p = **0.010**) and **0.01019** in the third cohort (p = **0.005**). Geographic distance also retained a positive association after BIO5 difference was removed, and this IBD-like component was larger in every cohort: **0.02760**, **0.02414** and **0.02225**, respectively (p = **0.005** in each). The mean IBE-minus-IBD contrast was negative in all three cohorts.

The BIO5 component is technically bounded. Although continuous flower-colour measurements reproduced it across all three cohorts, the validation flower-minus-background differential was not positive (mean rho = -0.00328, p = 0.82), and strict same-observer-pair sensitivities were unsupported in all three cohorts. The correct interpretation is therefore that the raw continuous measurements contain a reproducible BIO5-associated residual component, **not** that temperature-driven local adaptation has been demonstrated.

A broader hypothesis that environmentally heterogeneous species simply maintain more total flower-colour diversity was not supported. The primary multivariate heterogeneity test did not replicate discovery to validation, and a solar-radiation heterogeneity signal selected from the 500+500 screen failed fixed transport to the third cohort (partial rho = 0.0171, p = 0.382).
