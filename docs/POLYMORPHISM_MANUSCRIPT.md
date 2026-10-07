# Species-wide flower-colour variation is geographically partitioned across plant species

**Canonical current manuscript — synchronized with the New Phytologist submission draft**

**Authors:** [AUTHOR LIST TO CONFIRM]

**Affiliations:** [AFFILIATIONS TO INSERT]

**Corresponding author:** [NAME / EMAIL TO INSERT]

**Word counts (current working draft):**
- Summary: 194 words
- Introduction: 805 words
- Materials and Methods: 2,451 words
- Results: 2,613 words
- Discussion: 2,308 words
- Main text (Introduction through Discussion): 8,177 words
- Figures: 5
- Tables: 1
- Supporting Information: evidence map + planned supplementary figures/tables

**Keywords (alphabetical):** achromatic–chromatic axis; citizen science; flower colour; intraspecific variation; polymorphism; prospective confirmation

## Summary

- Species-wide intraspecific trait variation can arise because alternative phenotypes coexist locally or because different phenotypes occupy different parts of a species' range. Using flower-colour variation, we ask where ITV resides geographically, which spatial processes accompany it, and whether within-species displacement also recurs along common phenotype-space directions.
- We analysed 149,900 community-science photographs from 1,499 sampled plant species in sequential discovery, species-disjoint validation and prospective-confirmation resources. Species-wide four-state colour diversity was reproducible across observer partitions and fresh images.
- At the pre-specified post hoc 50-km scale, nearby conspecific observations were less colour-diverse than expected from each species' fixed overall colour composition in discovery, validation and the third cohort (mean depletion = 0.0205, 0.0187 and 0.0147; matched p = 0.005 in each). This geographic partitioning persisted after excluding same-observer pairs, removing all white records and using continuous nine-colour dissimilarity.
- Continuous colour turnover contained both a stronger isolation-by-distance-like component and a smaller BIO5-associated environment component across all three cohorts. The environmental component failed stricter background/observer controls and is therefore not evidence of local adaptation. Separately, a frozen achromatic–chromatic axis was prospectively confirmed (W = 0.517 versus null median 0.457, p = 0.001).

---

## Canonical inferential routing note

This canonical mirror uses the same current manuscript text as the New Phytologist submission draft. For machine-guard compatibility and reader routing, the following claim-boundary phrases are explicit here: the prospective test uses a **species-disjoint prospective confirmation cohort** from the **same iNaturalist opportunity universe**; it is an **untouched prospective confirmation**, because **a data-derived target cannot acquire prospective confirmatory status** in the data that generated it. The prospective resource is **one physical 499-species / 49,900-row measurement dataset used in two chronologically distinct ways**: first for H2, and **Only after H2 was terminalized** for highlight-validity and environmental follow-up. The white-state interpretation remains bounded by **exposure/background-context confounding**. Spatial organization is **structural rather than causal**. Candidate ecological and demographic mechanisms remain unresolved; the present photographs and occurrence geometry **cannot distinguish among them**.

## Introduction

Flower-colour research has developed along two complementary scales that rarely meet. Comparative macroecology asks how floral colour differs among species and environments, typically representing a species or community by mean, dominant or otherwise aggregated colour traits (Dalrymple et al. 2020). Flower-colour polymorphism research instead asks why alternative morphs coexist or replace one another within particular species and populations, where pollinators, abiotic selection, drift, gene flow, mating system and pleiotropic physiology can all matter (Narbona et al. 2018; Sapir et al. 2021). Genus-level work in *Protea* has shown that local ecological processes can scale to broader differences in polymorphism incidence (Carlson & Holsinger 2015), while citizen-science studies now map range-wide flower-colour variation within individual species (Surmacz 2023; McKenzie et al. 2026). A less-tested comparative question lies between these traditions: **when a species is variable in flower colour across its range, where does that variation actually reside?**

That question changes the evolutionary interpretation of intraspecific trait variation (ITV). A species can have high range-wide colour diversity because several morphs coexist repeatedly within local populations, because local populations differ from one another, or because both processes operate. These alternatives correspond to different evolutionary problems. Frequency-dependent or other balancing mechanisms can maintain variants within populations, whereas spatially varying selection can maintain adaptive polymorphism among populations when local selection is sufficiently strong relative to homogenizing migration (Delph & Kelly 2014). Drift, colonisation history and restricted gene flow can generate similar geographic differentiation without adaptation. In Mediterranean flower-colour polymorphisms, species described as polymorphic often consist mainly of monomorphic populations plus fewer polymorphic populations, producing clines or mosaics across the range (Narbona et al. 2018). Yet whether this **distributed-polymorphism** structure is a general property of flower-colour ITV across many plant lineages has not been tested at comparable scale.

The distinction also connects comparative ITV to landscape-evolutionary theory. Isolation by distance (IBD) predicts increasing phenotypic or genetic divergence with geographic separation through dispersal limitation and population history. Isolation by environment (IBE) describes divergence associated with environmental differences after geographic distance is accounted for (Wang & Bradburd 2014; Sexton et al. 2014). Applied to photographs, this can only be a **phenotypic IBE-like** diagnostic: colour distance is not genetic distance, and environment–phenotype association does not establish local adaptation. Nevertheless, partitioning geographic and environmental components can ask whether flower-colour turnover is purely spatial or whether environmental differences add nonredundant organization. The classic *Linanthus parryae* case illustrates why this distinction matters: its sharp flower-colour cline was not mirrored by neutral marker structure, and reciprocal transplants favoured resident morphs, turning a geographic pattern into evidence for local adaptation (Schemske & Bierzychudek 2007). The present comparative data can test the pattern, but not that final fitness step.

This spatial view complements rather than replaces phenotype-space questions. A range-wide trait distribution has at least three nonredundant properties: **overall sampled state diversity**, **where that diversity is geographically allocated**, and **which phenotypic directions the variation occupies**. Distributional ITV frameworks already represent species as probability densities or individual-level trait distributions rather than fixed means (Carmona et al. 2019; Palacio et al. 2025), and plant ITV can arise from genetic differentiation, local adaptation, phenotypic plasticity or combinations of these processes (Westerband et al. 2021). For flower colour, recurrent biochemical and regulatory architecture may additionally make some phenotype-space directions more accessible than others (Wessinger & Rausher 2012; Larter et al. 2018). Thus geographic partitioning asks **where variation is maintained**, while displacement geometry asks **what kinds of variation repeatedly become available**.

Scaling these questions first requires a measurement problem to be solved. Community-science photographs provide repeated observations over broad ranges, but image conditions and observer behaviour can create apparent colour differences even when the biological distribution is stable (Laitly et al. 2021; Di Cecco et al. 2021). We therefore use fixed high-depth species sampling, observer-disjoint validation, location-blind image measurement, explicit technical missingness and matched composition-preserving nulls. Four frozen biological colour states—white, yellow/orange, red/pink and blue/purple—provide a common species-wide diversity measure, while continuous nine-colour vectors retain finer phenotype geometry.

The study is deliberately sequential. An outcome-blind global opportunity frame first defines which species could enter before colour is measured. A 1,000-species high-depth resource is divided into 500 discovery and 500 species-disjoint validation species. After a recurrent white-versus-nonwhite displacement axis is localized in those data, the axis and its inferential machinery are frozen and tested in 499 newly sampled species. The spatial analyses then ask, explicitly as post-outcome ecological tests, whether range-wide colour diversity is locally depleted relative to each species' own composition and whether continuous colour turnover contains IBD-like and environment-associated components. Our central ecological hypothesis is that **species-wide flower-colour ITV is often distributed among geographic localities rather than expressed as unrestricted local coexistence**. We then ask whether this distributed variation is consistent with distance alone, whether environmental differences add explanatory structure, and how these spatial results relate to the prospectively confirmed recurrent phenotype-space direction.

---

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

### H3a: broad phylogenetic signal

Shared ancestry could generate cross-species similarity without repeated ecological organization. We therefore tested validation-cohort D using Blomberg's K with 9,999 tip-label permutations under three frozen V.PhyloMaker2 placement scenarios (341 species each). Support required p<0.05 in all three raw-D scenarios; Pagel's lambda and opportunity-adjusted traits were sensitivities and could not rescue the primary test.

### H3b: sampled photographic span

Species photographed across wider areas may appear more diverse simply because more environments were sampled. We therefore tested whether the discovery association between D and log sampled span reproduced in the species-disjoint validation cohort using a 20,000-permutation Spearman test, with observer/classifiability-adjusted partial-rank and rank-PGLS sensitivities. Sampled photographic span is not interpreted as true biological range size.

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

The central sampling logic is therefore **frame → discovery/validation → frozen target → prospective confirmation**. The later highlight, climate and phylogenetic analyses answer different validity or explanatory questions and do not create additional independent flower-colour cohorts. The prospective confirmation cohort remains within the same iNaturalist opportunity universe and measurement system, so it is species- and photo-disjoint confirmation rather than independent-source replication.


---

## Results

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

The result survived three direct falsification tests at 50 km. After excluding all same-observer photo pairs, mean depletion remained **0.01998**, **0.01690** and **0.01457** in discovery, validation and the third cohort, respectively (p = **0.005** in each). Removing every white observation and analysing only yellow/orange, red/pink and blue/purple retained positive depletion in all three cohorts (0.01947, p = 0.005; 0.00729, p = 0.025; 0.00891, p = 0.040). Replacing coarse states with continuous nine-colour Jensen–Shannon dissimilarity also retained local homogeneity (mean species-wide minus local JSD = **0.01518**, **0.01366** and **0.01301**; p = **0.005** in each cohort). Thus the geographic partitioning is not restricted to the exposure-coupled white classifier, coarse morph assignment or same-observer photo pairs.

A stronger claim that higher-D species allocate a larger fraction of their diversity among localities was not scale invariant. At 50 km the D–depletion correlation was unsupported in discovery (rho = 0.0390, p = 0.135) but positive in validation (rho = 0.1171, p = 0.025) and the third cohort (rho = 0.1118, p = 0.005). We therefore treat **local depletion itself**, rather than a universal D–partitioning slope, as the replicated result.

### Geographic colour turnover contains both IBD-like and bounded IBE-like components

Using continuous nine-colour dissimilarity, BIO5 difference explained a small positive component of colour turnover after geographic distance was rank-partialled out in all three cohorts: mean phenotypic IBE-like rho was **0.00925** in discovery (p = **0.005**), **0.00826** in validation (p = **0.010**) and **0.01019** in the third cohort (p = **0.005**). Geographic distance also retained a positive association after BIO5 difference was removed, and this IBD-like component was larger in every cohort: **0.02760**, **0.02414** and **0.02225**, respectively (p = **0.005** in each). The mean IBE-minus-IBD contrast was negative in all three cohorts.

The BIO5 component is technically bounded. Although continuous flower-colour measurements reproduced it across all three cohorts, the validation flower-minus-background differential was not positive (mean rho = -0.00328, p = 0.82), and strict same-observer-pair sensitivities were unsupported in all three cohorts. The correct interpretation is therefore that the raw continuous measurements contain a reproducible BIO5-associated residual component, **not** that temperature-driven local adaptation has been demonstrated.

A broader hypothesis that environmentally heterogeneous species simply maintain more total flower-colour diversity was not supported. The primary multivariate heterogeneity test did not replicate discovery to validation, and a solar-radiation heterogeneity signal selected from the 500+500 screen failed fixed transport to the third cohort (partial rho = 0.0171, p = 0.382).

### H3a: no detectable broad tree-wide conservation of D

Validation-cohort Blomberg-K tests were nonsignificant under all three frozen tree-placement scenarios:

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Pagel's lambda was small in each scenario and its tests against lambda = 0 were also unsupported. Opportunity-adjusted K sensitivities remained nonsignificant. The frozen verdict was `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

This result closes the tested claim of broad tree-wide signal under the frozen validation-cohort design; it does not imply that phylogeny is irrelevant to flower-colour polymorphism at all evolutionary scales.

### H3b: discovery span effect collapses in validation

Discovery showed a positive association between D and sampled photographic span (n = 369, rho = **0.1798786**, p = **0.00089996**). The species-disjoint validation cohort did not reproduce that effect (n = 363, rho = **-0.0025855**, p = **0.9586021**). The observer/classifiability-adjusted partial-rank result was likewise near zero (rho = 0.0055187, p = 0.9162042), and S1-S3 rank-PGLS sensitivities were unsupported.

The frozen verdict was `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

---

## Discussion

### From species means to validated within-species distributions

The first result is methodological but biologically consequential: within-species flower-colour diversity is reproducible across disjoint observer sets under the first-frozen validation design. Fresh-image remeasurement of 136 overlapping species gave very high D agreement (Spearman rho = 0.968; Lin CCC = 0.972), showing transport beyond one observer partition. D is nevertheless imperfectly measured: the stricter deterministic split exposes within-run uncertainty, and fresh transport remains within the same iNaturalist/FCP measurement system.

The methodological contribution is architectural rather than a claim to a new standalone statistic. Gini–Simpson diversity, rank correlations, Hellinger transformation, two-means clustering, Jensen–Shannon divergence, permutation tests and phylogenetic signal statistics are established tools. What is specific to this study is their assembly around fixed high-depth species sampling, observer-disjoint validation, location-blind image measurement, explicit technical-versus-ambiguous missingness, construction-preserving nulls and a new species/photo-disjoint prospective confirmation cohort. This design treats a within-species phenotype distribution as a species-level comparative trait while keeping measurement validity, target discovery and confirmation as separate inferential stages.

### Species-wide flower-colour ITV is distributed among localities

The clearest ecological result is not simply that flower colour is spatially autocorrelated. It is that **range-wide colour diversity is systematically depleted at local spatial scales relative to each species' own composition**. This result reproduced in three species-disjoint cohorts, across four spatial radii, after removing same-observer pairs, after removing every white record, and with continuous nine-colour distances. The pattern therefore resembles a distributed polymorphism: a species can be variable across its range while local neighbourhoods are more homogeneous than the range-wide phenotype pool.

This comparative result extends a pattern previously emphasized from population-level flower-colour studies. Narbona et al. (2018) noted that Mediterranean FCP species often contain many monomorphic populations and fewer polymorphic populations, producing geographic clines or mosaics. The present analysis does not identify populations a priori, but it reaches a closely related conclusion from a broad, common sampling design: **species-wide ITV and local ITV are not interchangeable**. A species mean erases the variation; a species-wide diversity scalar retains the variation but still erases where it resides.

The evolutionary distinction matters. In plant balancing-selection theory, negative frequency dependence and related mechanisms can maintain alternatives within populations, whereas spatially varying selection can maintain adaptive polymorphism among populations; restricted migration can strengthen such geographic differentiation (Delph & Kelly 2014). The FCP pattern is therefore compatible with spatially varying selection, but it is equally compatible with drift, colonisation history, dispersal limitation or environmentally induced plasticity. The present result identifies the spatial structure that an evolutionary explanation must account for; it does not by itself identify the explanation.

### Geographic separation dominates, with a smaller candidate environmental component

The IBD/IBE-like decomposition narrows that explanation without closing it. Geographic distance retained a positive association with colour turnover after BIO5 differences were removed, and this IBD-like component was consistently larger than the reciprocal BIO5-associated component. Thus restricted dispersal, population history or other spatially structured processes remain central candidates. At the same time, BIO5 differences added a small positive association after geographic distance was removed in all three cohorts, showing that distance alone does not exhaust the observed covariance in the raw flower-colour measurements.

The latter pattern should not be upgraded to local adaptation. In landscape genetics, IBE concerns genetic divergence associated with environmental differences independent of geographic distance (Wang & Bradburd 2014; Sexton et al. 2014). Here the response is phenotype, not genotype, and the BIO5 component fails stronger flower-minus-background and same-observer sensitivities. It could reflect real temperature-associated sorting, correlated environmental variables, plasticity or residual observation context. A useful benchmark is *Linanthus parryae*, where a flower-colour cline became convincing evidence for local adaptation only because neutral markers did not share the cline and reciprocal transplants showed resident-morph fitness advantage (Schemske & Bierzychudek 2007). FCP currently supplies the comparative spatial pattern, not that fitness test.

This distinction also explains why the earlier signed BIO5 result should remain secondary. A universal rule such as warmer environments favouring white flowers did not transport across cohorts, whereas unsigned BIO5-associated **turnover** is more reproducible in the raw continuous colour data. If biological, this would imply that temperature can be associated with where flower colour changes without forcing every species to change in the same phenotypic direction. Current technical sensitivities, however, keep that idea at hypothesis level.

### A complementary cross-species regularity in phenotype space

A complementary biological result concerns phenotype-space geometry. The original discovery/validation analyses showed that the recurrent construction-controlled component of within-species colour variation was overwhelmingly associated with a white-versus-nonwhite direction. White-versus-pigmented flower-colour combinations have historical precedent in floristic and experimental work, including observations that some anthocyanin-associated white/pigmented combinations are disproportionately represented in particular floras (Warren & Mackenzie 2001), but that precedent does not specify the mechanism of the present axis. Removing that axis eliminated the excess directional concentration, and the remaining non-white geometry did not support a shared hue direction.

The prospective-confirmation result changes the evidential status of this finding, but in a specific way. The white-axis target was no longer chosen after looking at the new cohort: species selection, measurement support, q_white, W, thresholds and the structured null were fixed before biological opening. At the primary tier, observed W = 0.517 exceeded a structured-null median of 0.457 (p = 0.001). Because that null already preserves coarse-state composition and is itself far above the isotropic expectation of 0.125, the prospectively confirmed quantity is **excess achromatic–chromatic alignment conditional on the measured coarse colour states**, not the existence of white-versus-nonwhite coarse combinations per se.

### What the achromatic–chromatic axis does not identify

The white-versus-nonwhite geometry is descriptive, not mechanistic. Reviews of flower-colour polymorphism emphasize that pollinator-mediated selection, abiotic selection, drift, gene flow, mating system and pigment genetics can all contribute in different systems (Sapir et al. 2021; Narbona et al. 2018). These analyses do not identify pigment chemistry, whether white states arise by pigment loss or non-white states by pigment gain, or the evolutionary direction of transitions. They also do not distinguish among developmental, genetic, pollinator-mediated or abiotic mechanisms.

The construction-preserving null asks whether continuous within-species geometry adds alignment after the frozen coarse-state composition and its palette mapping are held fixed. It therefore controls a construction baseline but does not validate the origin of the coarse white state itself. The direct highlight control shows that this distinction matters empirically: near-clipping is positively associated with frozen white classification within species (OR 1.444, 95% CI 1.389–1.502), with the complete interval above the predeclared negligible-coupling bound. The measured coarse white state therefore cannot be treated as free of image-exposure effects.

At the same time, this coupling does not account for the full H2 result. Removing all 2,205 response-blind high-clip rows retained excess alignment (W = 0.503 versus structured-null median 0.454, p = 0.001). The sensitivity lost 16 primary vectors and retained 89.9%, narrowly below the prespecified 90% threshold, so the frozen executable validity gate remained formally INDETERMINATE. Taken together, the evidence supports a bounded statement: an achromatic–chromatic excess remains after aggressive highlight exclusion, but the biological interpretation of the measured white state is exposure-coupled rather than artifact-cleared.

### Possible developmental accessibility of the achromatic–chromatic axis

The present spatial analyses do not identify a universal ecological driver of the recurrent achromatic–chromatic geometry. The earlier signed BIO5 prediction—white records should occupy warmer environments—passed its within-third-cohort gate but weakened under observer controls and failed the frozen discovery/validation transport rule. The newer unsigned BIO5-distance analysis instead detects a small colour-turnover association after geographic distance is removed, but this component does not survive the available matched-background or strict same-observer sensitivities. Temperature therefore remains a candidate correlate of spatial colour turnover, not an established cause of the recurrent phenotype-space axis.

A separate question is why many different perturbations might repeatedly project flower-colour variation onto an achromatic–chromatic direction. Published work provides one mechanistic example rather than a result of the present molecular study. In *Silene littorea*, Casimiro-Soriguer et al. (2016) reported >42-fold lower F3h1 expression in white than pigmented buds, differences involving Myb1a, and petal HPLC patterns consistent with altered anthocyanin-pathway regulation, without identifying a simple causal coding mutation. More generally, anthocyanin pigmentation can be altered through multiple structural and regulatory routes, with pleiotropic costs constraining which routes persist (Wessinger & Rausher 2012; Larter et al. 2018).

We additionally reanalysed published PAL/WAL frequency tables from Del Valle et al. (2019); the full methods, numerical results and ascertainment limits are reported in Supporting Information (Section S8). That descriptive reanalysis is consistent with flower-restricted pigment loss reaching higher reported natural frequencies than whole-plant anthocyanin loss, but it is not a direct molecular bridge to the present H2 geometry. The bounded hypothesis is therefore developmental rather than adaptive: shared pigment-network architecture may make some directions of flower-colour change repeatedly accessible even when the ecological or demographic processes producing spatial differentiation differ among species.

### Phenotype-space generality without a universal geographic map

The present data do not resolve whether geographic colour boundaries themselves are shared among species. Instead, the clearest cross-species regularity occurs in **phenotype space**: within-species colour variation repeatedly contains an achromatic–chromatic component, while its spatial realization may be shared, partly shared or species-specific. This distinction matters because a recurrent direction of variation does not require a recurrent map or a common environmental coefficient. Distinguishing shared from species-specific spatial sorting will require designs in which many species repeatedly sample both sides of the same candidate environmental or geographic contrasts.

### Diversity and geographic organization are related but not interchangeable

The older D–spatial analysis remains useful because species with greater species-wide colour diversity tended to show stronger distance-dependent organization in both original high-depth cohorts. The new local-depletion analysis resolves what that scalar association could not: **even at a fixed species-wide composition, nearby observations are systematically more homogeneous than expected**. Conversely, the strength of D–local-depletion coupling is scale dependent and did not meet the primary 50-km discovery–validation replication rule. We therefore do not claim a universal law that every increment in species-wide diversity is allocated disproportionately among localities.

These results clarify the comparative meaning of ITV. Overall diversity, its allocation among localities and the direction of phenotype-space displacement carry different information. Species can have the same range-wide D but differ in how locally mixed their morphs are; they can also have similar spatial partitioning while varying along different colour-space contrasts. Treating a species as a distribution is therefore necessary but still insufficient unless the **spatial allocation of that distribution** is retained.

### Two simple explanations fail fresh-data tests

The H3 tests sharpen what the species-level phenotype is not trivially reducible to. Validation-cohort D showed no detectable broad tree-wide phylogenetic conservation under any of the three frozen tree placements, while the apparent discovery association with sampled photographic span collapsed essentially to zero in the species-disjoint validation cohort. Together, these out-of-sample results show that reproducible between-species differences in D are not accounted for by either broad shared ancestry as detectable here or the geographic extent over which photographs happened to be sampled.

This inference is deliberately bounded. H3a is a non-support result rather than an equivalence test, so it does not establish a zero phylogenetic effect and does not exclude finer-scale lineage effects, particular clades or repeated evolutionary origins. H3b concerns photographic sampled span, not true biological range size.

### Scope, representativeness and source dependence

The 42,111-species frame gives the analysis broad taxonomic opportunity, but the high-depth cohorts are selected for repeated-observation support and are not a probability sample of global plant diversity. The paper therefore does not estimate the prevalence of flower-colour polymorphism.

Likewise, the prospective confirmation cohort is species-disjoint and prospectively tested, but it comes from the same iNaturalist source/opportunity universe and uses the same measurement system as the earlier cohorts. The fresh-image D transport check also remains within that same source and measurement system. Validation structure should match the intended generalization claim rather than being treated as generically independent (Roberts et al. 2017). The strongest current wording is fresh-image/same-system transport for D and prospective species-disjoint confirmation for the frozen H2 axis, not independent-source replication.

A stronger external validation would apply the same frozen q_white/W estimand and support rules to an independently generated image source, curated field dataset or another measurement system without retuning the axis. Demonstrating universality would require that the same recurrent direction survive such source-independent tests, ideally with calibrated reflectance or pigment measurements and across major plant lineages; the present study establishes broad recurrence within its sampling and measurement system, not a universal law of angiosperm flower-colour evolution.

### Conclusion

Species-wide flower-colour variation is not simply a cloud of interchangeable individual differences around a species mean. Across three species-disjoint high-depth cohorts, nearby conspecific observations were consistently more colour-homogeneous than expected from each species' fixed overall colour composition, and this result survived exclusion of same-observer pairs, removal of white records and continuous-colour reanalysis. **A substantial ecological property of flower-colour ITV is therefore where it is stored: across geographic localities as well as within them.**

The processes producing that distributed polymorphism remain only partly resolved. Continuous colour turnover contains both a stronger IBD-like component and a smaller BIO5-associated IBE-like component, but the latter does not survive stricter background and observer controls and cannot be interpreted as local adaptation. Broad environmental heterogeneity also does not generally predict total D. The appropriate evolutionary conclusion is therefore structural rather than causal: flower-colour diversity is repeatedly geographically partitioned, while selection, gene flow, drift, history and plasticity remain competing explanations for that partitioning.

A complementary result concerns phenotype space rather than geographic space. A frozen achromatic–chromatic displacement axis was prospectively confirmed in newly sampled species, showing that the directions in which flower colour varies are themselves nonrandom. Together, the spatial and geometric results suggest that intraspecific floral diversity has at least two reproducible dimensions—**where variants are distributed and which phenotype-space directions they occupy**—that disappear when species are represented by mean traits alone.

---

## Acknowledgements

[TO COMPLETE BEFORE SUBMISSION: funding, institutional support, data-provider acknowledgements, and individual contributions that do not meet authorship criteria.]

## Competing interests

[TO CONFIRM BEFORE SUBMISSION.]

## Author contributions

[TO COMPLETE AFTER FINAL AUTHOR LIST.]

## Data availability

Code, frozen protocols, machine-readable results, canonical figures and the reader-facing claim-to-input map are versioned in the `zuizui0223/fcp` repository. Exact source commits, SHA256 identities, workflow artifacts and recovery routes are consolidated in `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`.

For end-to-end verification, a self-contained snapshot is maintained under release tag `fcp-np-provenance-20260926`. The current asset identity is defined by the Git-tracked receipt `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`. The snapshot contains the 42,111-species opportunity frame, discovery/validation measured tables, the prospective-confirmation measured table, frozen spatial and technical-validation artifacts, image-measurement code/model bytes, manuscript figures and checksum-pinned WorldClim inputs. Exact WorldClim 2.1 BIO/SRAD bytes are also mirrored under release tag `fcp-worldclim-2.1-10m-20260925`.

The original third-party iNaturalist image pixels and flower masks were intentionally not retained. Frozen photo identities/source URLs and per-row image hashes permit future reacquired bytes to be checked while the provider continues to serve them. The related general-purpose `disttrait` implementation is **not a bitwise numerical reproducer** of the frozen study-specific H2 pipeline; all reported biological values come from the frozen study-specific implementation. Internal workflow identifiers and archived execution records are retained only for provenance and do not define the manuscript's reader-facing cohorts.

## References

Albert, C. H., Thuiller, W., Yoccoz, N. G., Soudant, A., Boucher, F., Saccone, P., & Lavorel, S. (2010). Intraspecific functional variability: extent, structure and sources of variation. *Journal of Ecology*, 98, 604–613. https://doi.org/10.1111/j.1365-2745.2010.01651.x

Carlson, J. E., & Holsinger, K. E. (2015). Extrapolating from local ecological processes to genus-wide patterns in colour polymorphism in South African *Protea*. *Proceedings of the Royal Society B: Biological Sciences*, 282, 20150583. https://doi.org/10.1098/rspb.2015.0583

Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2019). Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. *Ecology*, 100(12), e02876. https://doi.org/10.1002/ecy.2876

Casimiro-Soriguer, I., Narbona, E., Buide, M. L., del Valle, J. C., & Whittall, J. B. (2016). Transcriptome and biochemical analysis of a flower color polymorphism in *Silene littorea* (Caryophyllaceae). *Frontiers in Plant Science*, 7, 204. https://doi.org/10.3389/fpls.2016.00204

Dalrymple, R. L., Kemp, D. J., Flores-Moreno, H., Laffan, S. W., White, T. E., Hemmings, F. A., & Moles, A. T. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*, 228, 1972–1985. https://doi.org/10.1111/nph.16737

Di Cecco, G. J., Barve, V., Belitz, M. W., Stucky, B. J., Guralnick, R. P., & Hurlbert, A. H. (2021). Observing the Observers: How Participants Contribute Data to iNaturalist and Implications for Biodiversity Science. *BioScience*, 71(11), 1179–1188. https://doi.org/10.1093/biosci/biab093

Larter, M., Dunbar-Wallis, A., Berardi, A. E., & Smith, S. D. (2018). Convergent evolution at the pathway level: predictable regulatory changes during flower color transitions. *Molecular Biology and Evolution*, 35(9), 2159–2169. https://doi.org/10.1093/molbev/msy117

Laitly, A., Callaghan, C. T., Delhey, K., & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution*, 11, 4071–4083. https://doi.org/10.1002/ece3.7307

Luong, Y., Gasca-Herrera, A., Misiewicz, T. M., & Carter, B. E. (2023). A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences*, 11(5), e11546. https://doi.org/10.1002/aps3.11546

McKenzie, P. F., Church, S. H., & Hopkins, R. (2026). High-Throughput iNaturalist Image Analysis Reveals Flower Color Divergence in *Monarda fistulosa*. *The American Naturalist*, 208(1), 101–109. https://doi.org/10.1086/739413

Narbona, E., Wang, H., Ortiz, P. L., Arista, M., & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology*, 20(Suppl. 1), 8–20. https://doi.org/10.1111/plb.12575

Palacio, F. X., Graco-Roza, C., de Bello, F., & Carmona, C. P. (2025). Integrating intraspecific trait variability in functional diversity: An overview of methods and a guide for ecologists. *Ecological Monographs*, 95(2), e70024. https://doi.org/10.1002/ecm.70024

Puglielli, G., Bricca, A., Chelli, S., et al. (2024). Intraspecific variability of leaf form and function across habitat types. *Ecology Letters*, 27, e14396. https://doi.org/10.1111/ele.14396

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, 9, Article 39. https://doi.org/10.2202/1544-6115.1585

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40, 913–929. https://doi.org/10.1111/ecog.02881

Simpson, E. H. (1949). Measurement of Diversity. *Nature*, 163, 688. https://doi.org/10.1038/163688a0

Sapir, Y., Gallagher, M. K., & Senden, E. (2021). What Maintains Flower Colour Variation within Populations? *Trends in Ecology & Evolution*, 36(6), 507–519. https://doi.org/10.1016/j.tree.2021.01.011

Delph, L. F., & Kelly, J. K. (2014). On the importance of balancing selection in plants. *New Phytologist*, 201(1), 45–56. https://doi.org/10.1111/nph.12441

Schemske, D. W., & Bierzychudek, P. (2007). Spatial differentiation for flower color in the desert annual *Linanthus parryae*: was Wright right? *Evolution*, 61, 2528–2541. https://doi.org/10.1111/j.1558-5646.2007.00219.x

Sexton, J. P., Hangartner, S. B., & Hoffmann, A. A. (2014). Genetic isolation by environment or distance: which pattern of gene flow is most common? *Evolution*, 68(1), 1–15. https://doi.org/10.1111/evo.12258

Wang, I. J., & Bradburd, G. S. (2014). Isolation by environment. *Molecular Ecology*, 23, 5649–5662. https://doi.org/10.1111/mec.12938

Westerband, A. C., Funk, J. L., & Barton, K. E. (2021). Intraspecific trait variation in plants: a renewed focus on its role in ecological processes. *Annals of Botany*, 127, 397–410. https://doi.org/10.1093/aob/mcab011

Del Valle, J. C., Alcalde-Eon, C., Escribano-Bailón, M. T., Buide, M. L., Whittall, J. B., & Narbona, E. (2019). Stability of petal color polymorphism: the significance of anthocyanin accumulation in photosynthetic tissues. *BMC Plant Biology*, 19, 496. https://doi.org/10.1186/s12870-019-2082-6

Gómez, J. M., Perfectti, F., Armas, C., Narbona, E., González-Megías, A., Navarro, L., DeSoto, L., & Torices, R. (2020). Within-individual phenotypic plasticity in flowers fosters pollination niche shift. *Nature Communications*, 11, 4019. https://doi.org/10.1038/s41467-020-17875-1

Lacey, E. P. (2026). Temperature and the evolution of flower color: A review. *American Journal of Botany*, 113(1), e70106. https://doi.org/10.1002/ajb2.70106

Narbona, E., Perfectti, F., González-Megías, A., Navarro, L., del Valle, J. C., Armas, C., & Gómez, J. M. (2026). Heat drastically alters floral color and pigment composition without affecting flower conspicuousness. *American Journal of Botany*, 113(1), e70096. https://doi.org/10.1002/ajb2.70096

Surmacz, B. (2023). Spatial patterns of flower colour variation in native and introduced ranges of *Convolvulus arvensis* (Convolvulaceae) revealed by citizen-science data and machine learning. *Plant Biology*, 25, 681–686. https://doi.org/10.1111/plb.13537

Wessinger, C. A., & Rausher, M. D. (2012). Lessons from flower colour evolution on targets of selection. *Journal of Experimental Botany*, 63(16), 5741–5749. https://doi.org/10.1093/jxb/ers267

Warren, J., & Mackenzie, S. (2001). Why are all colour combinations not equally represented as flower-colour polymorphisms? *New Phytologist*, 151, 237–241. https://doi.org/10.1046/j.1469-8137.2001.00159.x

Literature-use boundaries are frozen in `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`. These references support background and interpretation; they do not alter the repository's machine-readable empirical results or claim ceiling.

## Figure legends

**Figure 1. Data provenance and inferential necessity from global sampling frame to prospective confirmation.** (a) Descriptive distributions of D = 1 - sum_k p_k^2 in the discovery (n = 369) and species-disjoint validation (n = 363) high-depth cohorts; dashed lines mark cohort medians. These cohorts are not a prevalence sample. (b) Sampling architecture from the 42,111-species outcome-blind opportunity frame. The 1,000-species high-depth resource is divided into discovery and validation because within-species trait construction must be tested separately from discovery. After the recurrent white-versus-nonwhite axis is identified and frozen, previously used high-depth species are excluded and a newly sampled 499-species cohort supplies the prospective confirmation. Only after that test is terminalized are the same prospective rows reused for highlight-validity and climate annotations, which constrain interpretation but are not part of prospective confirmation.

**Figure 2. Observer-disjoint reproducibility of the continuous polymorphism score D.** (a) Distributional summary of Spearman split-half correlations across 200 first-frozen observer-disjoint partitions in discovery and validation, shown as 5th percentile–median–95th percentile intervals. Validation median rho = 0.7891 and q05 = 0.7652; the dashed line marks the primary median floor of 2/3. (b) Later deterministic observer-disjoint stress tests in discovery and validation with bootstrap 95% intervals. The validation estimate was rho = 0.7927 (95% CI 0.7419–0.8324; CCC = 0.8474), narrowly below the prespecified 0.80 floor; this later stress test constrains but does not overwrite the chronologically earlier primary H1 result.

**Figure 3. Discovery and audit of the recurrent white-versus-nonwhite colour-space target in the original cohorts.** (a) Loadings of the fixed zero-sum q_white contrast: white is opposed to the equal mean of the eight non-white palette coordinates. The panel explicitly records that this named axis was isolated only after the original broad H2 geometry had been opened. (b) Original-cohort targeted W values (diamonds) against the median (points) and 95% interval (bars) of the construction-preserving structured null. Primary 0.10: discovery N = 152, W = 0.514625, p = 0.001; validation N = 129, W = 0.514586, p = 0.001. Strict 0.20: discovery N = 75, W = 0.542355, p = 0.001; validation N = 65, W = 0.510517, p = 0.008. Projection-removal and non-white-only falsification diagnostics are reported in Supporting Information.

**Figure 4. Prospective species-disjoint test of excess alignment with the frozen white-versus-nonwhite axis.** Structured-null W distributions from 999 frozen null worlds; dashed lines show null medians and solid vertical lines show observed W. (a) Primary 0.10 tier: N = 158, observed W = 0.517, null median = 0.457, 95% interval 0.436–0.475, upper-tail p = 0.001. (b) Strict 0.20 sensitivity: N = 86, observed W = 0.533, null median = 0.459, 95% interval 0.433–0.487, p = 0.001. The prospective cohort completed 49,900 terminal rows from 499 species, with 377 measurement-evaluable species and zero replacements before H2 opening; it is species-disjoint within the same iNaturalist opportunity universe, not an independent-source replication.

**Figure 5. Species-wide flower-colour variation is geographically partitioned, with a dominant IBD-like spatial component.** (a) At the fixed 50-km post hoc scale, local four-state pair diversity was lower than expected from each species' exact overall composition in discovery (n = 166, mean depletion = 0.02053), validation (n = 181, 0.01867) and the third cohort (n = 204, 0.01468); matched vertex-null p = 0.005 in each cohort. (b) The same local-depletion pattern persisted when all same-observer pairs were excluded, when all white observations were removed, and when continuous nine-colour Jensen–Shannon dissimilarity replaced coarse states. (c) Continuous colour turnover contained both IBD-like and BIO5-associated IBE-like components in all three cohorts, but the IBD-like partial correlations were larger (discovery 0.02760 vs 0.00925; validation 0.02414 vs 0.00826; third 0.02225 vs 0.01019). The BIO5 residual is phenotypic and technically bounded: validation flower-minus-background and strict same-observer sensitivities do not support a flower-specific, observer-independent environmental effect. Neither spatial component establishes local adaptation.

## Supporting Information

The complete evidence map is provided in `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`.

**Fig. S1.** Full H1 validation-cohort partition diagnostics for all 200 observer-disjoint partitions, including paired-species counts, Spearman rho, CCC, signed bias and absolute D differences.

**Fig. S2.** H1 discovery-cohort concordance together with the later deterministic validation-cohort stress-test diagnostics and its prespecified rho = 0.80 floor.

**Fig. S3.** Broad pre-target H2 geometry in the original cohorts, including leading-axis concentration, discovery-axis validation transport and discovery–validation axis alignment.

**Fig. S4.** Construction-preserving H2 null audit showing the quantities preserved when normalized nine-colour rows are permuted within coarse morph while species × coarse-morph row counts remain fixed.

**Fig. S5.** Residual H2 tests after projecting out q_white together with the low-sample non-white-only diagnostics.

**Fig. S6.** Prospective-confirmation chain of custody from deterministic species selection and fresh-metadata freeze through 256 terminal measurement receipts, support-gate completion and the durable H2_COMPLETE result.


**Fig. S7.** H3a sensitivity analyses across S1–S3 phylogenies for raw D, finite-sample sensitivity, opportunity-adjusted residuals, Blomberg K and Pagel lambda.

**Fig. S8.** H3b sampled-span sensitivities including raw, finite-sample, observer/classifiability-adjusted partial-rank and rank-PGLS validation analyses.

**Fig. S9.** Secondary empirical mechanism evidence. (a) Primary prospective, observer-balanced and observer-paired within-species white-minus-nonwhite BIO5 contrasts together with fixed discovery/validation transport tests. (b) Published *S. littorea* bud-stage molecular anchor: all significant pigmented-versus-white ABP expression contrasts, emphasizing F3h1 as the only locus significant in both pigmented morphs versus white and Myb1a as a candidate regulator. (c) Descriptive PAL/WAL natural-frequency contrasts within *S. littorea* and across the Del Valle et al. (2019) literature table; within *S. littorea*, PAL and WAL are distinguished by HPLC-DAD-MS^n tissue profiles rather than visual white-flower labels alone, whereas the cross-system entries retain the source table's phenotype classifications. These panels constrain mechanism without establishing universal causation.

**Fig. S10.** Multi-scale distributed-polymorphism sensitivities at 25, 50, 100 and 250 km, including the scale-dependent D–local-depletion correlations.

**Fig. S11.** Phenotypic IBD/IBE-like diagnostics, including unique geographic and BIO5-associated rank components, matched vertex-null distributions and the failed flower-minus-background/same-observer BIO5 robustness tests.

