# FCP frozen source-stage Methods preserved for New Phytologist Supporting Information

This is the complete parent PR #248 original-stage Methods section, preserved verbatim as historic reproducibility provenance. The current manuscript consolidates only the subset required for the full 4,475-species primary estimand; these underlying source algorithms and original 158-species prospective W and 1,761-species 50-km inference are unchanged.

**Source Git content SHA:** 5a6a440b349b459af16ffc6c5212da3c12fb530b

## Materials and Methods

The workflow is sequential because each stage addresses a different inferential failure mode: outcome-dependent species selection, observer-dependent measurement, target-selection overfitting, and post hoc causal interpretation. We first ask whether repeated photographs define a stable within-species colour phenotype, then whether its geometry recurs across species. A location-blind pipeline yields four-state frequencies for D and continuous nine-colour displacement; geography enters only after phenotypes are defined.

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



### Photographic measurement and outcome firewall

Geography is a biological outcome, so colour measurement was completed without access to geography or species-level spatial hypotheses, and failed rows could not be replaced after colour outcomes were known. All rows used a frozen location-blind pipeline (source commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`) with ROI-v4 flower detection, EfficientSAM segmentation, a normalized nine-colour palette and four coarse biological states. Failures became terminal nonclassifiable states.

For prospective confirmation, photo IDs were disjoint from previous rows, manifests were checksum-verified, and biological identity and geography remained sealed until all **256** terminal partitions completed. There was no early stopping or replacement, and image pixels/masks were not retained after sealing. Detailed firewall fields and hashes are in the prospective protocol and data-lineage map.

### Image-level biological states and continuous polymorphism score

The frozen biological coarse states are:

- white;
- yellow/orange;
- red/pink;
- blue/purple.

Rows classified as `mixed_uncertain` are not treated as a fifth biological colour state. Species enter D-based inference when at least 40 rows pass the global classifiability rule.

For a species with class proportions (p_k) over the four biological states, we calculate

[
D = 1 - sum_k p_k^2.
]

D is the Gini–Simpson diversity of the four frozen biological colour states (Simpson 1949). It combines richness and evenness and has a direct probability interpretation: under independent draws from the observed species-wide sample, D is the probability that two retained observations belong to different colour states. D therefore summarizes **overall sampled colour-state diversity across the observed species range**; it is not a within-population polymorphism measure and can reflect local coexistence, differentiation among sampled locations, or both. It is not converted to a global binary polymorphic/monomorphic outcome for the primary analyses.

Raw D remains the primary estimand. To assess finite-sample plug-in bias, we additionally use the finite-sample-corrected sensitivity

[
D_{\mathrm{corr}} = 1-\sum_k \frac{n_k(n_k-1)}{n(n-1)}
= D\frac{n}{n-1}.
]

Under independent multinomial sampling this correction targets the finite-sample bias of the plug-in Gini–Simpson estimator. It does not correct community-science sampling representativeness, observer effects, image formation or classifiability. With the frozen eligibility rule n >= 40, the largest possible multiplier is 40/39 = 1.0256.

### H1: observer-disjoint reproducibility

A comparative species phenotype must be repeatable across observation subsets rather than reflect particular observers. The primary H1 protocol therefore split observers, not photographs, into 200 deterministic balanced partitions without using colour outcomes or geography. D was estimated independently in each half; the validation rule required adequate paired-species support, median split Spearman rho >=2/3 and 5th-percentile rho >=0.5. CCC and absolute differences were agreement diagnostics.

A later deterministic split imposed a stricter rho >=0.80 stress test. Separately, fresh photo IDs were remeasured under the same frozen D definition to test same-system transport. The latter does not alter the first-frozen H1 decision and is not independent-source replication.

### H2: discovery of recurrent continuous colour geometry

Testing only pre-labelled categories would build the expected contrast into the analysis, so H2 first asked whether continuous within-species colour displacement shared any direction across species. Nine-colour compositions were normalized, Hellinger-transformed and partitioned within species by deterministic unlabeled two-means; the sign-invariant unit displacement was (u_i).

Two admissibility tiers required minor-mode fractions >=0.10 (primary) or >=0.20 (strict). The decisive construction-preserving null permuted normalized nine-colour rows across selected species within each frozen coarse morph while preserving every species × coarse-morph row count, then refitted the same construction.

The original discovery/validation geometry localized the recurrent component to the fixed white-versus-equal-nonwhite contrast:

[
q_{\mathrm{white}} = operatorname{normalize}(1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8).
]

For eligible species,

[
W = mean_i (u_i^T q_white)^2.
]

Because (q_{white}) was named only after the broad geometry was opened, these cohorts provide target localization rather than prospective confirmation. Projection-removal tests ask whether a residual recurrent hue direction remains.

### Prospective H2 confirmation in newly sampled species

Once the white-versus-nonwhite target had been identified, a new species/photo cohort was required to separate confirmation from target selection. Before biological opening we fixed species selection, no replacement, the location-blind pipeline, >=40-row species and >=250-species support gates, (q_{white}), W, the 0.10/0.20 tiers, a 999-replicate structured null, no axis refitting and one-shot durable terminalization.

Fresh retrieval yielded 499 species and 49,900 rows. All 256 partitions had to terminate before metadata–colour joining; insufficient support was predefined as underidentification. The result also had to serialize, validate as H2_COMPLETE and be committed immutably.

For each tier,

[
p = (1 + #(W_null >= W_obs)) / 1000.
]

With 999 null worlds, p=0.001 is the minimum plus-one Monte Carlo value (Phipson & Smyth 2010). Primary support required an evaluable 0.10 tier with p<0.05; the 0.20 tier was prespecified sensitivity.

### Post-confirmatory H2 validity diagnostics

A construction-preserving null does not show that the measured white state is free from image-exposure effects, so technical validity was evaluated only after H2 terminalization. We first compared observed W with isotropic and coarse-state-preserving baselines and re-applied the continuous gate within 299 additional null worlds.

For a direct exposure control, all 49,900 prospective rows were reacquired without replacement through the same pipeline. Before white/nonwhite outcomes were joined, flower-mask pixels were reduced to clip_fraction, near_clip_fraction and luminance_q99 and a response-blind high-clip set was frozen. Species-stratified conditional logistic regression tested within-species coupling between near_clip_fraction and white classification. The response-blind high-clip rule was near_clip_fraction > max(0.01, q95); the negligible-coupling interval was OR 0.80–1.25. A second sensitivity removed the high-clip set and reran the unchanged H2 construction, (q_{white}), W, 999 null worlds and seed; the executable gate required >=90% retention of the original 158 vectors.

### Post-confirmatory environmental filter and BIO5 transport test

A recurrent phenotype-space axis does not identify why it recurs. After H2 terminalization and high-clip freezing, prospective-cohort coordinates were annotated with WorldClim 2.1 maximum temperature of the warmest month (BIO5), precipitation of the driest month (BIO14) and mean solar radiation. The primary panel required >=5 white and >=5 non-white rows per species; 281 species qualified. Within-species white-minus-nonwhite contrasts were tested with Holm correction across the three variables and corroborated by species-stratified conditional logistic models including near-clip. The gate required >=100 species, the prespecified direction, Holm-adjusted p<0.05 and same-direction logistic p<0.05.

After the prospective BIO5 result was opened, the same fixed BIO5 contrast was transported to the original discovery and validation cohorts. Support required both cohorts to show the prespecified positive species-level and conditional-logistic results. This tests a simple abiotic explanation and whether it generalizes; it is not an independent-source causal test.

Two post hoc observer sensitivities bounded the prospective BIO5 association. Observer-paired analysis retained species-observer strata containing both colour states; observer-balanced analysis first averaged environment within species × observer × colour, then compared colour means across observers. These post hoc analyses cannot upgrade or redefine the frozen primary gate.

### Complementary test: species-level D and within-species geographic organization

To distinguish structured geographic variation from unstructured colour noise or broader photographic coverage, **all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity**. Species-level organization was

`rho_i = Spearman(d_geo_ij, d_colour_ij)`.

Within each species, complete colour vectors were permuted among fixed coordinates, preserving geographic geometry and colour-vector composition; each species had **999 matched within-species null values**. Across species, D was correlated with observed `rho_i`; each frozen null realization was analysed identically, with upper-tail probability

`(1 + # {rho_null >= rho_obs}) / 1000`.

For the primary robustness analysis, **Rank(D) and rank(`rho_i`) are separately residualized** on ranked sampled span and clear technical-failure rate before correlation. Validation additionally used paired flower/background measurements as `Spearman(d_geo_ij, d_flower_ij - d_background_ij)`; **It is not the difference between separate flower and background Spearman coefficients**. Exact diversity-minimizing/maximizing ambiguity completions provide endpoint stress tests. These analyses test a structural correlate, not a causal maintenance mechanism.

### Post hoc localization of species-wide diversity among geographic localities

The frozen D–spatial analysis established distance-dependent organization but did not directly distinguish local coexistence from among-locality differentiation. We therefore defined a post-outcome distributed-polymorphism diagnostic and applied it separately to discovery, validation and the third cohort. For each species, finite-sample species-wide pair diversity was the fraction of all unordered observation pairs assigned to different frozen coarse states. At radii 25, 50, 100 and 250 km, local diversity was the corresponding fraction among pairs separated by at most the radius; **50 km was fixed as the primary post hoc scale before these outcomes were opened**, and species required at least 30 local pairs.

Local depletion was

`D_pair - D_local(r)`.

Within each species, complete coarse-state labels were permuted among fixed coordinates 199 times. This matched vertex null preserves the exact state counts, sample size, geography and species-wide pair diversity, so positive depletion cannot arise merely because high-diversity species have more possible mismatched pairs. Cohort inference used the equal-species mean depletion against matched null means.

Three falsification analyses challenged the primary 50-km result: (i) all same-observer local pairs were excluded; (ii) every white record was removed and the test repeated on yellow/orange, red/pink and blue/purple only; and (iii) coarse states were replaced by continuous nine-colour Jensen–Shannon dissimilarity. Calendar-quarter-stratified and cross-year sensitivities additionally test whether seasonal or short-lived temporal sampling can account for local homogeneity. These analyses are post hoc and cannot alter earlier frozen decisions.

### Post hoc phenotypic IBD versus IBE-like decomposition

To distinguish geographic separation from environmental association, each species' continuous nine-colour Jensen–Shannon dissimilarities were paired with great-circle geographic distance and absolute WorldClim BIO5 difference. We calculated two symmetric partial-rank statistics:

`rho_IBE = Spearman(colour distance, BIO5 difference | geographic distance)`

and

`rho_IBD = Spearman(colour distance, geographic distance | BIO5 difference)`.

Complete colour vectors were permuted among fixed coordinates 199 times, preserving geography, BIO5 values, geography–BIO5 covariance and the exact colour-vector multiset. Species contributed equally to cohort means. We also report a descriptive rank-based commonality partition into unique geographic, unique BIO5 and shared components; pairwise rows are not treated as independent inferential replicates.

Because a positive phenotype–environment association can reflect image/background structure or observer geography rather than biological sorting, BIO5 interpretation was additionally challenged with continuous-colour, same-observer-pair and matched flower-minus-background sensitivities. We use **phenotypic IBE-like** only as a descriptive analogy to landscape-genetic IBE; these tests do not measure gene flow, fitness or local adaptation.

### Exploratory five-cohort environmental and colour-component analyses

Following the original prospective 281-species BIO5 test, we fitted a distinct **post hoc** within-species model to **1,574 species, 86,114 classifiable photographs and 727 genera** with at least five white and five chromatic records per species. Firth-penalized logistic models simultaneously included BIO5 (WorldClim 2.1 maximum temperature of the warmest month), BIO14 (driest-month precipitation) and mean solar radiation, standardized by each species' sampled climate SD. Two-stage normal-REML synthesis accounted for first-stage slope covariance; this is not a one-stage binomial mixed model. We independently report the equal-species arithmetic mean after transforming the slope and full covariance to units per 1°C, with 5,000 genus-cluster bootstrap draws, and a separate REML mean. Observer thinning, 1°-cell thinning and limited linear geographic adjustment were sensitivity analyses. The 1970–2000 WorldClim climatologies describe locations, not weather on the photograph date.

An additional **post hoc** composition check separated white–chromatic from chromatic–chromatic contributions to global-minus-50-km unequal-colour photo-pair frequency and divided each null-adjusted component by its own global mismatch within species where *both* denominators exceeded exploratory thresholds. All findings concern photographs and their selected climate normals rather than genetic morph frequencies or reproductive fitness; the earlier frozen BIO5 observer and cohort-transport decisions remain unchanged.

### H3a: broad phylogenetic signal

Shared ancestry could generate cross-species similarity without repeated ecological organization. We therefore tested validation-cohort D using Blomberg's K with 9,999 tip-label permutations under three frozen V.PhyloMaker2 placement scenarios (341 species each). Support required p<0.05 in all three raw-D scenarios; Pagel's lambda and opportunity-adjusted traits were sensitivities and could not rescue the primary test.

### H3b: sampled photographic span

Species photographed across wider areas may appear more diverse simply because more environments were sampled. We therefore tested whether the discovery association between D and log sampled span reproduced in the species-disjoint validation cohort using a 20,000-permutation Spearman test, with observer/classifiability-adjusted partial-rank and rank-PGLS sensitivities. Sampled photographic span is not interpreted as true biological range size.

### Reproducibility and frozen decisions

Frozen protocols, machine-readable results, input hashes and no-rescue rules govern interpretation. Prospective H2 used one authorized execution with no replacement, threshold changes, axis refitting or rerun-based selection after biological opening. Exact recovery routes are given in the data-lineage map and Supporting Information.

---

