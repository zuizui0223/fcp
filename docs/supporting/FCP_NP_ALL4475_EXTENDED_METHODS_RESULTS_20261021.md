# New Phytologist global 4,475-taxon methods and results: extended audit copy

**Frozen before Full Paper length editing (21 October 2026).** This keeps all analytical detail from the first 4,475-taxon post hoc manuscript reframe. The main paper can be shortened for the journal without discarding eligible-denominator definitions, post-freeze exploratory controls, earlier scientific failures or provenance. This is not an independent result source; source-locked JSON and the full 4,475-taxon CSV remain controlling.

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

### Frozen photographic observer-disjoint reproducibility test

A comparative species phenotype must be repeatable across observation subsets rather than reflect particular observers. The primary H1 protocol therefore split observers, not photographs, into 200 deterministic balanced partitions without using colour outcomes or geography. D was estimated independently in each half; the validation rule required adequate paired-species support, median split Spearman rho >=2/3 and 5th-percentile rho >=0.5. CCC and absolute differences were agreement diagnostics.

A later deterministic split imposed a stricter rho >=0.80 stress test. Separately, fresh photo IDs were remeasured under the same frozen D definition to test same-system transport. The latter does not alter the first-frozen H1 decision and is not independent-source replication.

### Frozen discovery of recurrent continuous colour geometry

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

### Prospective confirmation of the frozen white–chromatic geometry

Once the white-versus-nonwhite target had been identified, a new species/photo cohort was required to separate confirmation from target selection. Before biological opening we fixed species selection, no replacement, the location-blind pipeline, >=40-row species and >=250-species support gates, (q_{white}), W, the 0.10/0.20 tiers, a 999-replicate structured null, no axis refitting and one-shot durable terminalization.

Fresh retrieval yielded 499 species and 49,900 rows. All 256 partitions had to terminate before metadata–colour joining; insufficient support was predefined as underidentification. The result also had to serialize, validate as H2_COMPLETE and be committed immutably.

For each tier,

[
p = (1 + #(W_null >= W_obs)) / 1000.
]

With 999 null worlds, p=0.001 is the minimum plus-one Monte Carlo value (Phipson & Smyth 2010). Primary support required an evaluable 0.10 tier with p<0.05; the 0.20 tier was prespecified sensitivity.

### Post-confirmatory photographic-axis validity diagnostics

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

### Frozen exploratory test of broad phylogenetic signal

Shared ancestry could generate cross-species similarity without repeated ecological organization. We therefore tested validation-cohort D using Blomberg's K with 9,999 tip-label permutations under three frozen V.PhyloMaker2 placement scenarios (341 species each). Support required p<0.05 in all three raw-D scenarios; Pagel's lambda and opportunity-adjusted traits were sensitivities and could not rescue the primary test.

### Frozen exploratory sampled-range-span test

Species photographed across wider areas may appear more diverse simply because more environments were sampled. We therefore tested whether the discovery association between D and log sampled span reproduced in the species-disjoint validation cohort using a 20,000-permutation Spearman test, with observer/classifiability-adjusted partial-rank and rank-PGLS sensitivities. Sampled photographic span is not interpreted as true biological range size.

### Reproducibility and frozen decisions

Frozen protocols, machine-readable results, input hashes and no-rescue rules govern interpretation. Prospective H2 used one authorized execution with no replacement, threshold changes, axis refitting or rerun-based selection after biological opening. Exact recovery routes are given in the data-lineage map and Supporting Information.

### Complete 4,475-species post hoc observational estimands and missingness

To avoid making originally selected technical gates into unrelated ecological study populations, we constructed a master roster containing **each of the 4,475 source taxa** and all 100 original source photo records per taxon. Four-state diversity was finite-sample corrected whenever ≥2 classifiable photos were available (the observed minimum was eight), and fixed nine-colour-palette Hellinger two-mode white-axis alignment W was recomputed for every taxon. A newly exploratory, source-cohort × coarse-four-state preserving null shuffled full nine-colour vectors *across species only within those strata* and refitted the same two-mode axis in every one of 199 worlds, keeping the full 4,475 species unchanged in every replicate. We report that full-source W as a distinct post hoc estimator from the prospectively frozen 158-species W experiment. All unmatched genomic or direct biochemical annotations remained explicitly unmeasured.

We computed geographic pairwise colour-label mismatch using all classifiable geocoded photos of each species. For every species with one or more 50-km photo pairs, the estimator was the species' whole sampled-range pairwise mismatch minus the 50-km mismatch. For each of 199 Monte Carlo null worlds, colour labels were shuffled within the original species, leaving geography and the four-colour composition fixed. We report the 4,475-species denominators, species with no local pairs as missing, and opportunity sensitivity at ≥1/5/30/100 local pairs, as well as the exact originally frozen ≥30 pairs and ≥40 classifiable images gate. Every sample carries the same real-photo SHA-256.

Six WorldClim 2.1 bioclimatic normals (BIO1, BIO4, BIO5, BIO12, BIO14 and BIO15) were sampled at all source-valid photo sites from an authenticated frozen raster archive. We separately evaluated the full 4,475-taxon roster with **seven source-locked abiotic axes**: BIO1, BIO4, BIO5, BIO12, BIO14, BIO15 and the mean of 12 monthly WorldClim 2.1 solar-radiation normals. The within-taxon white-minus-nonwhite source-photo climate-location difference was scaled by that species' observed climate SD; undefined contrasts for taxa without both photo labels remained missing. The seven-axis family used a 1,999-resample genus-cluster, simultaneous maximum-absolute-studentised-statistic 95% bootstrap, with all 4,475 taxa in the registry (4,112 evaluable for contrasts). This exploratory univariate descriptive comparison is distinct from conditional multivariable Firth coefficients, flowering-date temperatures, and adaptive fitness (run #39456497702). An additional **adversarial photo-classifiability sensitivity** kept all 100 original photo opportunities per species but hypothetically assigned any of four coarse labels to previously unclassifiable images; sharp lower/upper pair-diversity bounds were computed by maximising/minimising the completed squared class-count sum. Exact extrema were tested exhaustively on small integer fixtures, with 4,475 genuine taxa analysed on the original source SHA (run #39518450603).

We also fit a post hoc two-parameter Firth-penalised logistic model \(\mathrm{logit}\Pr(\mathrm{photo\ white})=\alpha_i+\beta_i\mathrm{BIO5}\) independently within every species with both photographed colour labels, at least three climate-linked classifiable photographs, a BIO5 sampled standard deviation ≥0.05°C and range ≥0.1°C. Predictors were centred and scaled within each species for numerical fitting and coefficients were exactly transformed back to log odds per °C; fit failures remained explicit and no colour label was imputed. Species estimates were compared by species-equal mean, median, 10%-trimmed mean and source cohort, with 999 genus-cluster bootstrap draws conditional on completed fits. This full-opportunity **univariate** estimand is different from the previously reported selected-panel BIO5+BIO14+radiation model and cannot be equated to climate-based fitness or genetic selection. The descriptive source-wide observer-disjoint D check deterministically split independent observer IDs of all 4,475 species into disjoint halves without using flower-colour outcomes; its species-rank Spearman statistic is a *single split*, unlike the previously frozen 200-split validation gate. The complete-population primary descriptive signed environmental difference was the within-species mean BIO5/BIO12 of white-labelled minus nonwhite-labelled photo locations, requiring at least one photo in each class but never imputing a missing class. Species-level genus-cluster bootstrap intervals were calculated conditional on that selected support. Pairwise signed differences in four-state colour mismatch were also related descriptively to geography and BIO5 distance by partial Pearson correlations when at least 20 independent-source photo pairs supplied nonconstant response and predictors; these are not population-genetic IBD/IBE tests. Thus species remain in the same parent data frame while the estimable denominator can differ between *measurements*.

Local-morph photographic opportunities were audited in the full register as source-photo pairs ≤5 km apart, recorded in the same calendar year and month by distinct observers. Mere contrasting photo labels did not constitute evidence for inherited morphs. Molecular literature coverage was joined from a separately curated nonexhaustive taxon-ID table and explicitly distinguished from genotype, transcript and petal/leaf chemistry of the photographed individuals. No genetic, demographic or fitness outcomes were imputed to fill 4,475 rows.

### Post hoc joint test of repeated colour geometry against photographed geographic and climatic fates

To connect H1 with the geographic and abiotic estimands without changing the 4,475-species parent roster, we merged the four successful original-photo source-analysis artifacts by exact accepted iNaturalist taxon ID. The resulting permanent matrix has **4,475 rows × 114 source-verified columns**: continuous W geometry, categorical D, locality colour mismatch, source observer-disjoint D and univariate Firth BIO5, seven WorldClim photo-location climate contrasts, and sharp hypothetical bounds on unclassified photo labels. Counts, taxon identity, source cohorts and original classifiable/white photo numbers were verified for every joined table. The complete species CSV and its SHA-256 receipt are stored at `results/fcp_np_all4475_unified_species_ledger_20261021.csv` and `results/fcp_np_all4475_unified_species_ledger_receipt_20261021.json`.

Each pair of species-level traits was compared by ranked partial correlation after adjustment for D, white-label frequency, classifiable image count, independent observer count, original cohort and additional opportunity covariates relevant to that pair. Genus clusters were resampled for 1,999 bootstrap intervals on the **fixed nuisance residuals**. These conditional intervals ignore uncertainty in residual-model coefficients and do not address biological kinship beyond genus clustering or source-photo measurement reuse. Species with missing local or climatic estimates were retained in the 4,475-species master as NA. The four joint comparisons test empirical predictive coherence, not genetic isolation or a causal cost–benefit mechanism (source: `results/fcp_np_all4475_W_geography_climate_fate_20261021.json`).

### Supplementary post-freeze regional, seasonal and predictive falsification audits

After the original prospectively frozen white-axis test, we examined whether photographed geographic sorting could be a by-product of colour-label composition, photo purity, or heterogeneous observation times. In a separate, source-locked regional analysis we allocated classifiable conspecific records to equal-area cells; 1,282 taxa and the same eligible cells could be compared in full and strict flower-ROI-colour-purity photographs. Four-state between-cell total variation was partitioned exactly into absolute white-frequency difference and the residual nonwhite-composition component. For each component, a source-species label-permutation null preserved the original photo-state multiset and cell counts. The relative excess of the white component versus the nonwhite component was calculated **only** when both components had non-negligible expected null differences; thresholds of 0.005, 0.01 and 0.02 were evaluated post hoc. Calendar-month and exact-year-month constrained label permutations were performed separately. These are technical sensitivity and observation-structure controls, not tests of pigment biosynthesis or morph fitness (source: PR #249; run 39425739017).

We also distinguished the *mean conditional association* of flower photographs with climate from prediction at independently withheld locations. The expanded five-cohort BIO5/BIO14/radiation Firth and normal-REML species-slope summaries, ten-fold leave-genus-out common-effect score test, and one-degree spatial-cell and 50-km-buffer withheld-site predictions constitute different prediction targets. These were analysed post hoc on a selected 1,574-species, 86,114-photograph panel and must not be treated as a new preregistered test of a temperature-adaptation mechanism (source: PR #248 and the stored transport model outputs).

A further nested geographic climate comparison used frozen WorldClim 2.1 long-term BIO5 normals on the true 447,500-photograph input. Among same-species, same-calendar-month comparisons in five-degree geographic cells, we decomposed the white-versus-other-temperature contrast algebraically into differences among nested one-degree cells and differences remaining within those nested cells. It was recalculated for original and high-purity photo labels; the conditions excluded cases with insufficient photo-label or independent-observer opportunities. Long-term climate normals were not mistaken for weather at flowering (source: run 39430364175). These exploratory extensions are scientifically separate from the preregistered continuous-colour W test.

---

## Table 1. Source design and distinct inferential targets

| Component | Authentic measurement / comparison | Principal inference | Boundary |
|---|---|---|---|
| Photographic source | Five disjoint species cohorts; 4,475 taxa × 100 photos = 447,500; 227,107 classifiable | Range-wide photographed phenotype distributions | Not a global census of heritable morphs |
| Frozen prospective axis | 158 evaluable species; fixed W and 999 structured-null worlds | Repeated white–chromatic phenotype-space direction | Same photo provider; no genotype |
| Geographic sorting | 1,761 taxa with fixed 50-km local pair opportunity; 199 original species-conditional permutations | Nonrandom local vs whole-range photo-colour composition | Does not identify selection or history |
| Climate and pair-distance | 1,574 taxa / 86,114 photos post hoc; three original cohorts for continuous IBD-/IBE-like | Conditional climate association and geographic distance structure | Not climate-driven fitness or genetic IBE |
| Technical/narrative falsifications | Month/observer/photo-ROI purity, component-baseline nulls, 1° spatial holdout, nested geography | Which empirical regularities survive specified filters | Exploratory after original W freeze; no independent source |

---

## Table 2. Five evolutionary hypotheses and the highest available level of evidence

| Hypothesis and distinguishing prediction | What the current multi-species photographs establish | Decisive independent test still required |
|---|---|---|
| **H1. Repeated phenotypic accessibility:** floral displacement directions recur beyond label-composition nulls | **Supported for source-measured phenotype space:** frozen prospective W across 158 species, with a technically indeterminate stringent exposure gate | Independent calibrated flower-organ reflectance and genotype/morph identity; not another reshuffling of the same photo set |
| **H2. Historical allocation:** neutral ancestry and past connectivity order colour differences even among climatically matched sites | **History unidentified:** phenotypic IBD-like signal exceeds BIO5 IBE-like, but geographic distance is not ancestry | Neutral SNP distances, independent range-history constraints and environment-matched population-pair predictions |
| **H3. Conditional ecological sorting:** present habitat predicts colour after history, and fitness depends on morph × environment | **Photographic correlation only:** small positive BIO5 association, weakened background/sampling controls and negative species-specific spatial transfer | Neutral-ancestry-matched environments plus common-garden reproductive or survival-fitness differences |
| **H4. Local maintenance:** frequency dependence, variable selection or migration retains inherited morphs within breeding populations | **Not identified:** strict 5-km month-matched photo test is unsupported; biological C/S crosswalk lacks S-only spatial matches | Marked plants, inherited morph verification, pollen/seed fitness and repeated years in true interbreeding populations |
| **H5. Tissue-specific pigment constraint:** petal-restricted biochemical changes persist differently from whole-plant losses | **Published case evidence only:** petal-only versus whole-plant loss in Silene; not the photographed individuals | Paired petal/leaf chemical assays, genotypes and fitness across independent origins |

These five hypotheses are an **evolutionary synthesis proposed after the principal photographic tests**. They must not be presented as five pre-registered positive tests or five experimentally verified causes.

---

## Results

### Repeated directions of colour variation across species

Five outcome-blind photo cohorts provided 447,500 original measurements from 4,475 species. Of these, 227,107 photographs could be assigned to one of the four frozen flower-colour labels. **All 4,475 species had at least eight classifiable photos**, so each could contribute a continuous nine-palette geometric estimate and four-state finite-sample-corrected diversity, with technical uncertainty retained. The previously selected ≥40-classifiable gate admitted 3,297 species; it is now a sensitivity/confirmatory eligibility stratum rather than the global observational denominator. These are photo-label data selected for sampling opportunity, not a prevalence survey of genetic polymorphism in natural populations. Within the complete 4,475-species frame, the finite-sample-corrected species-equal four-state diversity was D = **0.29187** (median 0.26869); 4,189 species had at least two coarse colour labels and 4,112 contained white and nonwhite photo labels. These counts are **not** frequencies of genetic colour polymorphism in natural populations. Under the former ≥40-classifiable eligibility, only 3,297 species remained (3,029 with multiple photographic colour labels), and their mean D was 0.23751; the 1,178 species with fewer than 40 classifiable images had mean D = **0.44401**. The classifiability–D Spearman association across all 4,475 was **−0.58989**; a genus-resampled 95% interval for the low-minus-high-classifiability D difference (+0.20650) was [+0.19068, +0.22237]. This could reflect differential ROI/classification errors or genuine differences in flowers and image acquisition; the finite-sample correction does not make classification missingness random.

An exact missing-colour sensitivity analysis retained all 4,475 taxa and treated the **220,393 unclassified original photographs** as unknown rather than zero-diversity observations. If every such photograph were hypothetically assignable to one of the four original colour classes, the **sharp attainable bounds** on the species-equal mean all-100-photo mismatch D would be **0.15282–0.67063**, compared with **0.29187** estimated *conditional on classified photos*. The extrema allocate missing colours to the most common state (minimum D), or balance the counts (maximum D). These are not confidence intervals or biological population-colour bounds: an unclassifiable source photo may fail the flower-organ assumption entirely. They demonstrate that global D estimates are not insensitive to photo-label missingness (Fig. 1; Source workflow #39518450603).

A further source-opportunity audit showed that apparent regional trends are sensitive to the species and photographic opportunities sampled. An unpaired, equal-species 0–15° absolute-latitude subset had greater photo-colour diversity (D = 0.3315; 146 taxa) than the 60–90° subset (D = 0.2438; 98 taxa), but the same 246 taxa recorded across both tropical and extratropical zones yielded only a mean D difference of +0.0121 (95% species bootstrap interval −0.0058 to +0.0298). Thus the unpaired global gradient does not establish a replicated within-species tropical polymorphism advantage. An exploratory strict photo-colour-purity sensitivity (≥65% flower ROI in the assigned state and ≤20% in every alternative state) reduced mean categorical D from 0.1154 to 0.0647 in the *same* 1,813 photo-eligible species, underlining the importance of flower-part identity and mixed-colour ROI classification.

The complete-population geometry refit applied the identical Hellinger two-mode algorithm to every taxon, including all 1,178 species previously below the 40-photo gate. Across all 4,475 species the descriptive mean white-axis squared alignment W was **0.478872** (genus-bootstrap 95% interval **0.46874–0.48978**). The complete global post hoc null shuffled entire nine-palette vectors only within each of the five original source cohorts and frozen four-state labels, preserving each taxon's photo-label counts, identity and measurement opportunity, and refitted the Hellinger two-mode axis of all **4,475 species in each of 199 worlds**. The observed W exceeded the matched-null median **0.432039** (95% null range **0.42847–0.43587**; upper P = **0.005**). Each of the five species-disjoint source cohorts separately exceeded its own matched null at Monte Carlo P = 0.005. Under the original 40-photo/second-colour/minor-cluster gates, **1,217 species** remained, W = **0.519455** versus matched null median **0.494150** (P = 0.005). Thus the new all-species construction-preserving test supports a recurring *photo-phenotype* direction without discarding the 1,178 low-classifiability taxa. It was designed after the source data and phenotype axis were opened; it is exploratory, not a replacement for the separately authorised frozen prospective result in 158 species.

The independent-observer reliability check demonstrated that the sampled colour diversity is not arbitrary noise: the **all-4,475-species** source-wide deterministic observer-ID split yielded a between-split species-level Spearman rank correlation for D of **0.82259**, with **all 4,475 species estimable**. The earlier frozen 363-species validation stress test yielded split-half Spearman 0.7927 (bootstrap 95% interval 0.7419–0.8324) and concordance 0.8474. However, this only establishes within-provider photographic repeatability; it does not distinguish genetic morphs from consistent flowering-stage, petal-structure or observer-source effects.

At the primary 0.10 tier, **158** species produced admissible nonzero displacement vectors. Observed W was **0.517**. Across 999 structured-null worlds, the null median was **0.457** and the 95% interval was **0.436–0.475**. Observed W was 1.13 times the null median, with upper-tail p = **0.001**.

At the strict 0.20 tier, **86** species were eligible. Observed W was **0.533**, compared with a null median of **0.459** and a 95% interval of **0.433–0.487**. Observed W was 1.16 times the null median, again with p = **0.001**.

The frozen terminal verdict was therefore

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

This is an untouched prospective test of a previously frozen axis in a species-disjoint cohort. Because the prospective confirmation cohort was drawn from the same iNaturalist source/opportunity universe and processed with the same measurement system, it is not described as an independent-source replication.

The result is a prospective *species-disjoint* test of a phenotype-space direction, not an independent imaging-platform replication and not a reconstruction of repeated molecular transitions. The frozen null already had a substantial baseline white-axis alignment (median W = 0.457), making a composition-preserving comparison indispensable. The independent image-highlight sensitivity retained 142/158 primary species under its stringent clipping gate; that technical decision remained indeterminate rather than an unconditional proof that every white photo is an independent white flower morph. No analogous recurrent nonwhite-hue axis was supported by the frozen analyses.

### Geographic redistribution is detectable but heterogeneous

A newly source-authenticated post hoc analysis retained **all 4,475 source species** in a single row-level table, permitting local mismatch to be marked missing rather than excluding the species. Of these, **3,600** had at least one pair of classified photos within 50 km; **875** had no such pair and therefore had no observed local mismatch estimate. The species-equal mean global-minus-local difference among the 3,600 with any local pair was **+0.016116** (199 matched photo-label permutations, null mean −0.0001093, upper P = **0.005**; observed median +0.008726). Requiring five or more local pairs retained 3,178 species (mean +0.017610), 30 or more retained 2,084 (+0.017421), and 100 or more retained 1,117 (+0.013746); all were above their opportunity-matched nulls at Monte Carlo P = 0.005. The original 40-classifiable and ≥30-local-pair eligibility produced exactly **1,761 species and +0.017614902**, thereby reproducing the previously frozen numerical result from the same complete source without dropping taxa from the parent registry.

The final fifth cohort added **2,478 new species and 247,800 original measured photographs**, bringing the five-cohort atlas to **4,475 species, 447,500 photographs**. Of 227,107 classifiable photographs, the unchanged >=40-usable rule retained 3,297 species; **1,761 species** additionally met the original >=30 within-50-km pair criterion. In this species-equal panel, mean global four-state photo-pair mismatch was **0.224305** compared with **0.206690** within 50 km: a global-minus-local depletion of **0.0176149** (matched label-permutation null mean 0.0000446; 199-permutation upper-tail p = **0.005**). The fifth cohort separately showed mean depletion **0.017599** in 1,020 eligible species, consistent in sign and scale with the first four cohorts. These are photographic differences rather than estimates of genetically segregating morph frequencies.

The full pooled species median was **0.006871**, and **984 of 1,761 species (55.88%)** had positive depletion, versus a matched-null mean of **44.58%**. Both statistics exceeded the matched permutation null (Holm p = **0.020** within the four *post hoc* pooled distributional checks). Removing the highest-effect 20% in both the observations and every resampled null left an **observed mean of −0.006247**, still greater than the equally trimmed null mean **−0.018278** (within-family Holm p = **0.020**). Accordingly, the aggregate signal is not confined to the highest-effect species, but neither do these results show every species to be spatially differentiated. The separate stricter 5-km, same-site/month and independent-observer analysis did *not* support local within-population polymorphism; the broader 50-km result does not reverse that inference (Fig. S12).

The geographic effect replicated across species-disjoint cohorts under the original matched-null procedure and was not restricted to colour-label contrasts involving white. Coarse-label within-50-km depletion persisted after excluding same-observer pairs and white photos and after replacing categorical mismatch with continuous nine-colour Jensen–Shannon dissimilarity in the three original cohorts (details and statistics in the frozen supporting analyses). These comparisons detect spatial organisation of photographed phenotype labels, **not** local persistence of two inherited morphs. The stricter within-5-km, same-year-month, distinct-observer test did not show additional local colour sorting relative to its matched null (194 eligible species, P = 0.300).

### Recurrent phenotype-space geometry does not predict a common spatial or climatic outcome

To test the study's central distinction **on the same source species**, rather than juxtaposing estimates from differently gated populations, we analysed the permanently joined **4,475 × 114** taxon-by-estimand table. We compared each species' Hellinger two-mode achromatic–chromatic W value with its photo-colour global-minus-local 50-km mismatch and its own BIO5 Firth slope or BIO4 standardised mean location contrast. In each comparison, we residualised species-rank variables against the same source-photo classifiability, D, white-label fraction, observer opportunities, source cohort and, where relevant, local-pair count or climate SD. We then calculated genus-cluster bootstrap intervals (1,999 repetitions) *conditional on these source-based nuisance residualisations*. This was planned after all original measurements and is an **exploratory, same-provider association**, not a new independent genetic or environmental test.

Across all **3,600 species** with geographic photo-pair opportunity, W was very weakly related to the amount of 50-km spatial mismatch depletion after the source-photo controls (partial rank r = **+0.03375**, genus-cluster 95% interval **[−0.00226, +0.06920]**). Across **4,112 taxa** with a species-specific Firth BIO5 slope, the corresponding conditional association between W and the signed temperature response was **−0.00927 [−0.03874, +0.02151]**. W and the standardised BIO4 white-versus-other climate-location contrast likewise had **−0.02484 [−0.05423, +0.00498]** (4,112 taxa). Omitting each source cohort in turn did not generate a large W–BIO5 association (five leave-cohort-out partial rank r between −0.0144 and +0.0047). These intervals constrain large *observed rank associations* conditional on the source controls, but do **not** establish literal evolutionary independence, equivalence, or absence of environment-specific selection.

A distinct comparison of the geographical photo-colour depletion with the species-level BIO5 photo-location logistic slope was weakly positive in the 3,292 species where both were estimable (conditional partial rank r = **+0.04888 [ +0.01580, +0.08160 ]**). This value is small, noncausal, and may partly reflect the common original photograph labels; it is not a history-versus-climate variance decomposition. Together with the negative white-component-normalisation and weak spatial-transfer results, the central outcome is an **empirical separation of a recurrent phenotype direction from strong, transferable geographic/climatic predictions**, rather than evidence that all species follow one evolutionary fate.

### The repeated white–chromatic direction does not have a demonstrable spatial advantage

Exact partitioning of the full 1,761-species effect yields **0.0147632** from unequal white–chromatic photograph pairs and **0.00285166** from chromatic–chromatic pairs. Their sum equals **0.0176149**; white–chromatic contrasts therefore constitute **83.81%** of the aggregate geographic mismatch depletion. However, they already account for **79.57%** of global unequal-pair mismatch, and *both* components individually exceed their matched within-species label permutations (exploratory two-component Holm p = **0.010** for each).

The abundance of white-labelled differences alone does not establish a special ecological advantage of a white phenotype. Among **636 species** with each type of global unequal-pair mismatch above 0.01, the respective normalized null-adjusted depletion averages were **0.05696** (white–chromatic) and **0.06945** (chromatic–chromatic); the paired difference was **−0.01249**, with genus-cluster 95% interval **[−0.06813, +0.04591]**. This difference was not resolved, and the direction was sensitive to the exploratory baseline threshold (Fig. S13). The recurrent continuous white–chromatic displacement direction established by the original prospectively tested W statistic is a different, independently defined phenotype-space property; neither W nor this absolute decomposition estimates pigment biosynthetic costs or natural selection.

At a second, equal-area between-region scale, a matched comparison of 1,282 taxa in identical geographic cells showed a small positive species-equal photographic sorting excess under both broad and strict photo-colour-purity rules (strict total-variation excess = 0.00727; matched 199-world upper P = 0.005). The white-versus-other-photo component contributed 89.5% of that small *absolute* excess, but white colour contrasts were already frequent in the source. Under a component-specific null-opportunity denominator ≥0.01 in both components, the putative *additional proportional* spatial enrichment for white was unsupported (174 strict-purity taxa, excess difference +0.0934, 95% bootstrap interval −0.1160 to +0.3186; matched P = 0.220). The broad photo subset was likewise unsupported (395 taxa, P = 0.370). The 89.5% value therefore cannot be interpreted as preferential spatial selection for white flowers.

Year-month conditioning also reduced the strict-purity regional sorting excess (from 0.007263 unrestricted to 0.003515 under exact-year-month shuffling in 1,281 matched species), although the pooled matched null remained positive (P = 0.005). This source-dependent effect did not attain significance in three of the five source cohorts, and removal of the large new-U100 cohort weakened the pooled comparison (P = 0.075). The inference is a detectable but heterogeneous same-provider photo-label pattern, not a cross-provider law of historical or adaptive white-flower differentiation.

### Geographical colour separation and climate association are distinct, small signals

In the complete 4,475-species registry, **4,112 species** had both white-labelled and nonwhite-labelled photographs at sites with valid BIO5/BIO12 normals. Without imposing the previous ≥40-classifiable / ≥5 photos-per-state gate, their species-equal mean white-minus-nonwhite BIO5 contrast was **+0.128217°C** (499-draw genus-bootstrap 95% interval **+0.07905 to +0.17844°C**) and the corresponding BIO12 annual-precipitation contrast was **−4.237 mm**. The very small observations per class in some taxa require caution, and these are within-species photographic *mean climate-location differences*, not conditional Firth logistic slopes. An additional univariate within-species Firth logistic regression of photographic white labels against BIO5 was estimable in **4,112 of all 4,475 species**, with the other 363 species remaining in the complete source table as one-photo-colour-only cases. The mean log-odds slope was **+0.021888 per 1°C** (999-draw genus-bootstrap 95% interval **+0.00846 to +0.03489**); median +0.01851, 10%-trimmed mean +0.01982, and 54.45% of fitted species had a positive slope. Source-cohort means were −0.00495 in discovery, −0.00925 in validation, +0.01838 in prospective, +0.03390 in wave4 and +0.03199 in new-U100. This is a weak, heterogeneous average directional relation, **not** a five-cohort shared-direction result. The prior selected 1,574-species temperature regression jointly conditioned on other climate variables and retains a different modelling target; it is a more restrictive multivariate comparator, not the common global denominator.

In the full-source seven-variable abiotic comparison, **4,112 of the original 4,475 taxa** supplied a white–chromatic climate-location contrast. After a simultaneous seven-axis, genus-cluster 95% adjustment, only **BIO4** (temperature seasonality: +0.03237 within-species SD, 95% interval **+0.00421 to +0.06053**) and **BIO5** (maximum temperature of the warmest month: +0.05069 SD, **+0.02596 to +0.07542**) excluded zero. BIO1 (+0.01854), BIO12 (−0.00948), BIO14 (+0.00750), BIO15 (−0.02068) and solar radiation (−0.01029) had adjusted intervals spanning zero. The BIO4 direction was not homogeneous across source cohorts: the prospective and new-U100 cohort means were positive (+0.09286 and +0.03921 SD), while discovery and validation were −0.00281 and −0.00483 SD. This climate-normal seasonality index measures average *within-year* temperature variation, not year-to-year selective fluctuations. These are photographic, univariate post hoc correlations, with no evidence here of a genotype-by-environment fitness mechanism.

Also from the complete source, four-state photo-pair mismatch yielded descriptive partial Pearson correlations with geographical distance and BIO5 separation for **4,189 taxa** where both components were identifiable (mean geography-only partial r = **+0.01358**, BIO5-only partial r = **+0.00645**). This exploratory *coarse-label, individual-photo-pair* estimator differs from the earlier three-cohort continuous nine-colour Jensen–Shannon partial-rank statistic and is not an independent source replication or a test of genetic isolation.

The distinct **post hoc** five-cohort analysis of **1,574 species and 86,114 photographs** found a small positive average association between white classification and BIO5 (maximum temperature of the warmest month), after simultaneous adjustment for BIO14 and mean solar radiation within each species. The equal-species mean log-odds slope in common units was **+0.03233 per 1°C**, with 5,000-draw genus-bootstrap 95% interval **[+0.01786, +0.04728]**. The separate two-stage REML synthesis gave **+0.01708 per 1°C** (modified-HKSJ 95% interval **[+0.00897, +0.02518]**); on each species' sampled-climate-SD scale its mean was **+0.05467** [0.03149, 0.07785]. Each number refers to a distinct summary estimand, not a percentage-point difference in white-flower prevalence.

BIO14 and solar radiation had mean intervals spanning zero. The small BIO5 direction persisted in specified observer-thinned, spatially thinned and geographically adjusted sensitivity fits, but these do not fully eliminate spatial, seasonal, cultivation, lighting or observer biases. In particular, this larger-sample *exploratory* mean association does **not** retrospectively pass the earlier *frozen* two-cohort BIO5 transport gate, and it does not demonstrate white-morph fitness advantage in warm environments (Fig. S14).


Using continuous nine-colour dissimilarity, BIO5 difference explained a small positive component of colour turnover after geographic distance was rank-partialled out in all three cohorts: mean phenotypic IBE-like rho was **0.00925** in discovery (p = **0.005**), **0.00826** in validation (p = **0.010**) and **0.01019** in the third cohort (p = **0.005**). Geographic distance also retained a positive association after BIO5 difference was removed, and this IBD-like component was larger in every cohort: **0.02760**, **0.02414** and **0.02225**, respectively (p = **0.005** in each). The mean IBE-minus-IBD contrast was negative in all three cohorts.

The BIO5 component is technically bounded. Although continuous flower-colour measurements reproduced it across all three cohorts, the validation flower-minus-background differential was not positive (mean rho = -0.00328, p = 0.82), and strict same-observer-pair sensitivities were unsupported in all three cohorts. The correct interpretation is therefore that the raw continuous measurements contain a reproducible BIO5-associated residual component, **not** that temperature-driven local adaptation has been demonstrated.

A broader hypothesis that environmentally heterogeneous species simply maintain more total flower-colour diversity was not supported. The primary multivariate heterogeneity test did not replicate discovery to validation, and a solar-radiation heterogeneity signal selected from the 500+500 screen failed fixed transport to the third cohort (partial rho = 0.0171, p = 0.382).

Despite the positive mean association, generalisation to withheld data was weak or negative. When held-out genera were used with a pooled BIO5 component, mean within-species conditional-pair log-loss improvement was only +0.000699 nats (AUC = 0.516). In contrast, species-specific climate models trained on other geographic cells performed worse than an intercept-only model on withheld one-degree cells (mean species-equal log-loss gain = −0.059924 nats); adding climate to a species-specific geography model also deteriorated withheld-cell loss (−0.149777 nats). The corresponding 50-km-buffer comparisons were negative as well (−0.073453 and −0.237316 nats). These two prediction programmes use different estimators and loss definitions and are not exchangeable effect estimates. They show that a weak average within-species climatic association need not yield robust prediction at novel localities.

In an additional exploratory same-month climate audit on the original source, strict-purity white-labelled photographs were recorded at a mean BIO5 normal 0.32985°C higher than nonwhite photos within the same species and five-degree cell (327 taxa). Exact decomposition gave +0.27158°C between nested one-degree subcells and +0.05827°C within those subcells; the one-degree conditional permutation for the within-subcell component did not reach P < 0.05 (P = 0.075). This spatial-scale decomposition is *not* a division of historical and adaptive causes. It does not establish species-independent thermal benefits of flower depigmentation.

### Observational evidence does not yet identify morph maintenance or a shared molecular origin

All **4,475 taxa** were also examined for geographically close, same-year-and-month, **distinct-observer** colour-labelled photo pairs. **840 species** had at least one qualifying pair; **325 species** showed at least one photographic colour-label difference in such pairs. These are candidates for blind petal-part/flower-age validation, **not** evidence of stable inherited sympatric morphs, because separate observers within 5 km and a common month need not sample the same interbreeding population or different biological morphs. The matched strict locality null remained unsupported in the previously frozen subset.

The complete taxon register was joined by verified taxon ID to the separate, curated published-molecular-evidence inventory. Only **11 taxa** had any curated external molecular/colour mechanism evidence and **nine** were annotated at its T2-or-higher level. The remaining taxa are **not missing a molecular mechanism biologically**; they simply have no matched source evidence in this nonexhaustive literature screen. None of the 4,475 photographic individuals had paired phenotype, genotype, petal/leaf pigment chemistry and fitness in this source.

A separate positive-evidence literature table identified 34 species in which local colour-morph coexistence, geographic segregation or both had been independently documented. Only 16 matched the 4,475-species source-photograph roster, and only five matched the 1,761-species 50-km analysis, including **zero** with literature-documented geographic segregation *only*. Thus a quantitative comparison of the biological local-maintenance versus spatial-segregation states cannot be estimated from the current photograph-plus-literature overlap without strong ascertainment bias.

Published tissue chemistry and gene-expression studies, including petal-limited and whole-plant anthocyanin-loss states in *Silene littorea*, provide independent molecular examples of non-equivalent phenotypes that can appear white to a camera. They are not individual-matched genetic data from our 4,475-species photographic dataset, and no pigment investment, seed fitness or common-garden genotype-by-environment effect was measured in the global dataset.

---


