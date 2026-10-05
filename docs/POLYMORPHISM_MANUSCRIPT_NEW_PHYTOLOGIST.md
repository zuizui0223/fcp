# Within-species flower-colour variation shows recurrent achromatic–chromatic geometry across plant species

**New Phytologist Full Paper — submission-format working draft**

**Authors:** [AUTHOR LIST TO CONFIRM]

**Affiliations:** [AFFILIATIONS TO INSERT]

**Corresponding author:** [NAME / EMAIL TO INSERT]

**Word counts (current working draft):**
- Summary: 181 words
- Introduction: 756 words
- Materials and Methods: 2,361 words
- Results: 2,538 words
- Discussion: 1,816 words
- Main text (Introduction through Discussion): 7,471 words
- Figures: 5
- Tables: 1
- Supporting Information: evidence map + planned supplementary figures/tables

**Keywords (alphabetical):** achromatic–chromatic axis; citizen science; flower colour; intraspecific variation; polymorphism; prospective confirmation

## Summary

- Species means can erase the structure of intraspecific trait variation. Using flower-colour polymorphism, we ask whether within-species phenotype distributions are reproducible comparative traits, whether they vary along recurrent phenotype-space directions, and whether greater diversity is more strongly organized geographically.
- We analysed 100,000 photographs from 1,000 discovery–validation species and then prospectively tested a pre-frozen white-versus-nonwhite axis in 49,900 newly sampled photographs from 499 species. A construction-preserving null retained species × coarse-state counts and coarse-state-specific palette distributions.
- Observer-disjoint validation recovered stable D rankings (median Spearman rho = 0.789), and fresh-image transport across 136 overlapping species was strong (rho = 0.968; Lin CCC = 0.972). In the prospective test, 158 species gave W = 0.517 versus null median 0.457 (p = 0.001).
- Greater D was associated with stronger within-species geographic colour organization. Thus cross-species generality was stronger in phenotype-space geometry than in a universal environmental coefficient: the BIO5–white association weakened under observer controls and failed transport, whereas bounded *Silene* molecular and PAL/WAL persistence evidence supported pigment-network accessibility as a mechanistic clue. White classification remained exposure-coupled (OR = 1.44).


---

## Introduction

Flower-colour research has developed along two complementary scales that rarely meet. Comparative macroecology asks how floral colour differs among species and environments, typically representing a species or community by mean, dominant or otherwise aggregated colour traits (Dalrymple et al. 2020). Flower-colour polymorphism research instead asks why distinct morphs coexist or replace one another within particular species and populations, where pollinators, abiotic selection, drift, gene flow and mating system can all matter (Narbona et al. 2018; Sapir et al. 2021). Genus-level work in *Protea* has shown that local ecological gradients can scale to broader differences in polymorphism incidence (Carlson & Holsinger 2015), while recent citizen-science studies have mapped range-wide colour variation at high throughput within individual species (Surmacz 2023; McKenzie et al. 2026). A less-tested comparative question lies between these traditions: **does the distribution of within-species flower-colour variation itself show repeatable structure across many species?**

That gap matters because intraspecific trait variation has at least three separable properties: **how much** variation a species contains, **which phenotypic directions** that variation occupies, and **how the variants are arranged geographically**. These need not covary: equal diversity can be geographically intermixed or partitioned into mosaics and clines, and species with similar diversity can vary along different phenotype-space contrasts. Distributional ITV frameworks already represent species as probability densities or individual-level trait distributions rather than fixed mean values (Carmona et al. 2019; Palacio et al. 2025). What remains less tested is whether distinct properties of the same within-species distribution—its **amount, direction and geographic organization**—carry independent, repeatable cross-species structure. Flower-colour polymorphism is unusually suited to that test because discrete morph frequencies, continuous colour displacement and geographic sorting can be measured within one trait system, while pigment pathways provide a bounded route to developmental interpretation.

Scaling that comparison first requires a measurement problem to be solved. Community-science photographs offer repeated individuals across broad spatial extents, but image conditions and observer behaviour can create apparent within-species diversity even when the underlying biological distribution is stable (Laitly et al. 2021; Di Cecco et al. 2021). Our methodological contribution is therefore not a new diversity index or clustering algorithm. It is an inferential architecture for treating a repeatedly sampled within-species distribution as a comparative species trait: fixed high-depth sampling, observer-disjoint validation, location-blind measurement, explicit technical missingness, construction-preserving nulls, and separation of target discovery from a species-disjoint prospective confirmation cohort. For four frozen biological colour states—white, yellow/orange, red/pink and blue/purple—we summarize the amount of within-species diversity as (D = 1-sum_k p_k^2), while retaining continuous nine-colour information to study the geometry of variation.

This design allows two ecological alternatives to be distinguished. If flower-colour variation is largely idiosyncratic across species, continuous displacement should not repeatedly concentrate along the same phenotype-space direction. If common developmental, genetic or functional constraints make some contrasts more accessible than others, a recurrent axis may appear even when the ecological causes differ among species. Independently, if larger D mainly reflects sampling breadth or unstructured noise, it need not be associated with geographic organization; if ecological or demographic sorting contributes to species differences in D, more diverse species may also show stronger within-species spatial structure without requiring the same geographic boundary or environmental driver across species.

The design is deliberately sequential because each stage solves a different inferential problem. First, an outcome-blind global metadata frame defines which species could have entered the study before flower colour is examined; this prevents the comparative sample from being assembled around visually interesting outcomes. Second, a 1,000-species high-depth resource is divided into discovery and species-disjoint validation cohorts so that a within-species diversity phenotype can be tested for reproducibility and recurrent colour-space structure can be discovered without treating discovery as confirmation. Third, because the white-versus-nonwhite axis was identified only after the original colour geometry had been opened, that axis, statistic, thresholds and structured null are frozen and tested in a newly sampled species- and photo-disjoint prospective confirmation cohort; a fresh cohort is necessary for the resulting test to have confirmatory rather than retrospective status. Finally, phylogeny, sampled photographic span and climate are added as targeted annotations or follow-up tests, not as new flower-colour cohorts, to ask whether broad ancestry, sampling extent or a simple environmental rule can explain the observed structure. The study therefore has two distinct aims: **methodologically**, to establish when repeated photographs support a validated species-level distributional trait; and **ecologically**, to ask whether cross-species generality resides in what phenotypic contrast varies and in how strongly that variation is geographically sorted, rather than necessarily in one universal driver or one shared map.

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

D is the Gini–Simpson diversity of the four frozen biological colour states (Simpson 1949). It combines richness and evenness and has a direct probability interpretation: under independent draws from the observed state distribution, D is the probability that two observations belong to different colour states. Because every species is represented on the same four-state scale, D provides a common bounded measure of the amount of within-species colour diversity. It is not converted to a global binary polymorphic/monomorphic outcome for the primary analyses.

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

### Secondary *Silene littorea* molecular-anchor extraction

For a molecular PAL anchor, we structured reported Casimiro-Soriguer et al. (2016) results from nine morph × stage mRNA-seq samples and 29 anthocyanin-pathway loci, retaining significant bud-stage contrasts involving white petals, regulatory contrasts, sequence follow-up and petal HPLC. The frozen extraction is source-derived, not raw-read reanalysis or replication.

### Secondary PAL/WAL persistence reanalysis

To ask whether tissue restriction could plausibly influence the persistence of achromatic pigment-loss phenotypes, we performed a descriptive reanalysis of Del Valle et al. (2019) Supplementary Tables S2 and S1. This analysis was secondary and post-publication; it was not an untouched confirmatory test.

For *Silene littorea* Table S2, frequencies of petal anthocyanin-loss (PAL) and whole-plant anthocyanin-loss (WAL) phenotypes were parsed across 21 populations and five survey years. For each phenotype we summarized positive population-years, the number of populations ever positive, positive-frequency range and median, and persistence across repeated years. Zero values remained zero and missing cells were not imputed.

For the cross-system Table S1 comparison, all **13 PAL** and **13 WAL** entries were retained in the source registry. Exact numeric frequencies were represented as [x,x], reported ranges as [a,b], and WAL entries reported as <x as censored intervals [0,x]. Qualitative “rare” or “extremely rare” WAL entries were retained without numeric imputation, and an explicitly greenhouse WAL entry was excluded from the natural numeric panel while remaining in the source registry. We report PAL lower- and upper-bound summaries, WAL numeric upper-bound summaries and counts of conservative interval separations. No cross-system p-value is used because study ascertainment, survey effort and reporting format are heterogeneous. The analysis tests a **differential-persistence plausibility** prediction; it does not estimate mutation rates or establish that tissue restriction causally raises natural frequency.

We also resolved PAL/WAL taxa against FCP cohorts to assess whether a direct D/H2 bridge was estimable. Source strings and taxonomic resolutions were retained, greenhouse-only records were excluded from the natural-frequency bridge, and sparse overlap was treated as non-estimable rather than tested.

### Complementary test: species-level D and within-species geographic organization

To distinguish structured geographic variation from unstructured colour noise or broader photographic coverage, **all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity**. Species-level organization was

`rho_i = Spearman(d_geo_ij, d_colour_ij)`.

Within each species, complete colour vectors were permuted among fixed coordinates, preserving geographic geometry and colour-vector composition; each species had **999 matched within-species null values**. Across species, D was correlated with observed `rho_i`; each frozen null realization was analysed identically, with upper-tail probability

`(1 + # {rho_null >= rho_obs}) / 1000`.

For the primary robustness analysis, **Rank(D) and rank(`rho_i`) are separately residualized** on ranked sampled span and clear technical-failure rate before correlation. Validation additionally used paired flower/background measurements as `Spearman(d_geo_ij, d_flower_ij - d_background_ij)`; **It is not the difference between separate flower and background Spearman coefficients**. Exact diversity-minimizing/maximizing ambiguity completions provide endpoint stress tests. These analyses test a structural correlate, not a causal maintenance mechanism.

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
| Published *S. littorea* transcriptome/biochemistry | Casimiro-Soriguer et al. (2016); structured extraction of reported Table S4 / main-text RNA-seq, sequence follow-up and HPLC results | 29 ABP-related loci across nine morph × stage samples; expanded sequencing survey of 38 individuals | Anchor one high-frequency PAL system to a quantified pigment-pathway signature without treating source data as new sequencing | Secondary source-derived molecular evidence |
| Published pigment-loss frequency tables | Del Valle et al. (2019) Supplementary Tables S2 and S1; deterministic descriptive parsing | 21 *S. littorea* populations; 13 PAL + 13 WAL literature systems | Ask whether flower-restricted pigment loss can persist at higher natural frequencies than whole-plant loss | Secondary descriptive mechanism evidence |

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

### Published molecular data anchor the *Silene* PAL phenotype near F3h1/Myb1a

The published *S. littorea* bud-stage transcriptome identified **F3h1** as the only locus significantly higher in both pigmented-versus-white contrasts among 29 ABP-related loci. F3h1 expression was **49.0×** higher in dark-pink than white buds (p = **0.039**) and **42.2×** higher in light-pink than white buds (p = **0.049**). **Myb1a** was **5.1×** higher in dark-pink than white (p = **0.009**) and **4.2×** higher in dark- than light-pink buds (p = **0.021**).

Sequence evidence did not identify a simple coding lesion: the source reported 622 SNPs across 29 ABP loci, but **F3h1 had zero SNPs** in the reported UTR/CDS table. Nine colour-associated *Ans* SNPs were synonymous, and expanded sequencing of **38 individuals** found no consistently colour-differentiating SNP. Petal HPLC showed cyanidin derivatives and white-versus-pigmented differences in rutin, quercetin and isovitexin. Together these measurements are consistent with a regulatory blockage near F3h1, potentially involving Myb1a, but identify no causal mutation. This is a single-species PAL anchor, not molecular validation of H2 (Fig. S9b).

### Biochemically anchored flower-restricted anthocyanin loss reaches higher reported natural frequencies

To test whether tissue restriction could plausibly influence the persistence of achromatic phenotypes, we reanalysed the frequency tables of Del Valle et al. (2019). Within *Silene littorea*, PAL and WAL are not visual labels alone in that source study: HPLC-DAD-MS^n tissue profiling showed petal anthocyanin-loss (PAL), with anthocyanins absent from petals but retained in photosynthetic tissues. Whole-plant anthocyanin-loss (WAL) lacked anthocyanins in both petals and photosynthetic tissues; flavone production was retained across phenotypes. Within *Silene littorea*, PAL whites were recorded at **8–21%** when present (median **15.5%**) across six positive population-years in two populations, both observed positive for at least three years. WAL whites occurred in more populations but remained at **0.05–0.86%** in every positive population-year (median **0.21%**). Every positive PAL frequency exceeded the maximum positive WAL frequency.

The broader literature table gave the same qualitative contrast. Across **13 PAL** and **13 WAL** systems, the median PAL lower bound was **5%**, whereas the median numeric WAL upper bound was **0.1%** and the largest quantified WAL upper bound was **1.4%**. Seven of 13 PAL lower bounds exceeded that 1.4% maximum, and every PAL upper endpoint exceeded it. These are descriptive reanalyses of an ascertained, heterogeneous literature sample with censored and qualitative frequency reporting; they do not estimate an unbiased cross-species effect or prove that tissue restriction itself causes higher frequency. They do, however, provide empirical support for differential persistence of flower-restricted versus whole-plant pigment-loss phenotypes (Fig. S9b).

A taxonomy-resolved overlap audit found only **two PAL and two natural-context WAL species** in the high-depth FCP cohorts, all in validation and none in discovery/prospective. No PAL-versus-WAL test of D or H2 geometry was therefore estimable.

### Greater D is associated with stronger within-species geographic colour organization

The positive D–spatial association reproduced across the two species-disjoint high-depth cohorts. The raw association was rho = **0.0892133** (p = **0.034**) in discovery and rho = **0.1016008** (p = **0.025**) in validation.

After controlling for sampled geographic span and clear ROI/flip technical-failure rate, the geometry-preserving analysis remained positive in discovery (partial rho = **0.1266367**, p = **0.007**) and validation (partial rho = **0.0992877**, p = **0.025**). In validation, the matched flower-minus-background response was also positive (partial rho = **0.1162411**, p = **0.010**).

The validation result remained supported under exact uniform ambiguity-endpoint completions: primary D_min4 rho = **0.0970781** (p = **0.029**) and D_max4 rho = **0.1252858** (p = **0.008**); flower-minus-background p-values were **0.009** and **0.006**. Species with greater measured flower-colour diversity therefore tend to show stronger internal geographic organization, although this association does not identify the causal process that creates or maintains that organization.

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

### Cross-species regularity in phenotype space

The strongest positive biological result is geometric. The original discovery/validation analyses showed that the recurrent construction-controlled component of within-species colour variation was overwhelmingly associated with a white-versus-nonwhite direction. White-versus-pigmented flower-colour combinations have historical precedent in floristic and experimental work, including observations that some anthocyanin-associated white/pigmented combinations are disproportionately represented in particular floras (Warren & Mackenzie 2001), but that precedent does not specify the mechanism of the present axis. Removing that axis eliminated the excess directional concentration, and the remaining non-white geometry did not support a shared hue direction.

The prospective-confirmation result changes the evidential status of this finding, but in a specific way. The white-axis target was no longer chosen after looking at the new cohort: species selection, measurement support, q_white, W, thresholds and the structured null were fixed before biological opening. At the primary tier, observed W = 0.517 exceeded a structured-null median of 0.457 (p = 0.001). Because that null already preserves coarse-state composition and is itself far above the isotropic expectation of 0.125, the prospectively confirmed quantity is **excess achromatic–chromatic alignment conditional on the measured coarse colour states**, not the existence of white-versus-nonwhite coarse combinations per se.

### What the achromatic–chromatic axis does not identify

The white-versus-nonwhite geometry is descriptive, not mechanistic. Reviews of flower-colour polymorphism emphasize that pollinator-mediated selection, abiotic selection, drift, gene flow, mating system and pigment genetics can all contribute in different systems (Sapir et al. 2021; Narbona et al. 2018). These analyses do not identify pigment chemistry, whether white states arise by pigment loss or non-white states by pigment gain, or the evolutionary direction of transitions. They also do not distinguish among developmental, genetic, pollinator-mediated or abiotic mechanisms.

The construction-preserving null asks whether continuous within-species geometry adds alignment after the frozen coarse-state composition and its palette mapping are held fixed. It therefore controls a construction baseline but does not validate the origin of the coarse white state itself. The direct highlight control shows that this distinction matters empirically: near-clipping is positively associated with frozen white classification within species (OR 1.444, 95% CI 1.389–1.502), with the complete interval above the predeclared negligible-coupling bound. The measured coarse white state therefore cannot be treated as free of image-exposure effects.

At the same time, this coupling does not account for the full H2 result. Removing all 2,205 response-blind high-clip rows retained excess alignment (W = 0.503 versus structured-null median 0.454, p = 0.001). The sensitivity lost 16 primary vectors and retained 89.9%, narrowly below the prespecified 90% threshold, so the frozen executable validity gate remained formally INDETERMINATE. Taken together, the evidence supports a bounded statement: an achromatic–chromatic excess remains after aggressive highlight exclusion, but the biological interpretation of the measured white state is exposure-coupled rather than artifact-cleared.

### Why an achromatic–chromatic axis may recur

One plausible source of recurrence is genetic and developmental accessibility. Anthocyanin-based pigmentation can be reduced through multiple structural and regulatory changes, allowing distinct molecular routes to converge on pale or white petals; the realized routes are shaped by pleiotropic costs (Wessinger & Rausher 2012). The recurrent achromatic–chromatic direction is therefore consistent with a many-to-one accessibility bias, but H2 neither establishes anthocyanin deficiency nor identifies transition direction.

Two linked evidence layers sharpen the pigment-pathway interpretation. In *S. littorea*, F3h1 is >42-fold lower in white buds than both pigmented morphs, Myb1a also differs, no causal colour-differentiating SNP survives expanded sequencing, and petal HPLC is concordant with an F3h-region blockage. Our PAL/WAL reanalysis adds persistence: flower-restricted loss reaches much higher reported natural frequencies than whole-plant loss. This fits an accessibility-and-maintenance filter, but the molecular experiment is single-species, the frequency systems are ascertained and only two PAL plus two natural WAL systems overlap FCP; no direct bridge to D/H2 is estimable.

Second, the environmental follow-up gives a context-dependent abiotic clue rather than a robust temperature rule. White records occupied warmer BIO5 environments in the frozen primary prospective analysis (median contrast +0.069 SD; Holm-adjusted p = 0.0354; conditional OR = 1.073, p = 0.000919), but the signal weakened when observer identity was controlled more aggressively: the observer-paired sensitivity was null (median 0.000 SD; p = 0.485; OR = 0.787, p = 0.202), while the observer-balanced contrast remained positive but marginal by Wilcoxon (median +0.054 SD; p = 0.075). The association also failed the frozen species-disjoint transport rule: discovery was null (p = 0.743) and validation missed the species-level criterion (p = 0.0541). Temperature is therefore not supported as a universal cross-species driver, and some of the within-cohort signal may reflect observer-associated geographic sampling. This fits evidence that thermal effects on floral pigmentation are context dependent (Lacey 2026) and that *Moricandia arvensis* can shift reversibly from lilac to white while losing detectable anthocyanins under summer conditions (Gómez et al. 2020; Narbona et al. 2026). Generality may therefore lie in shared pigment-network architecture acted on by different genetic and environmental perturbations, rather than in one universal BIO5 coefficient.

### Phenotype-space generality without a universal geographic map

The present data do not resolve whether geographic colour boundaries themselves are shared among species. Instead, the clearest cross-species regularity occurs in **phenotype space**: within-species colour variation repeatedly contains an achromatic–chromatic component, while its spatial realization may be shared, partly shared or species-specific. This distinction matters because a recurrent direction of variation does not require a recurrent map or a common environmental coefficient. Distinguishing shared from species-specific spatial sorting will require designs in which many species repeatedly sample both sides of the same candidate environmental or geographic contrasts.

### More colour diversity is more geographically organized

The replicated D–spatial association shows that species with greater measured colour diversity also tend to have stronger internal geographic organization. This is not explained by sampled span in validation, persists after sampled-span and technical-failure adjustment, remains positive in a matched flower-minus-background contrast, and survives ambiguity-endpoint stress tests.

The result is structural rather than causal. Spatially varying abiotic selection, pollinator turnover, restricted dispersal or gene flow, demographic history, drift and mating-system differences could all contribute. Together with H2, it supports a **two-layer ecological question**: **what varies** shows recurrent cross-species structure, whereas **where that variation is sorted** remains unresolved and potentially context dependent.

These results also show what is lost when species are represented by a mean or by a single scalar measure of variability. Distributional ITV approaches already establish the value of retaining within-species trait distributions; the present result adds an empirical decomposition of the **same distribution** into **amount, direction and spatial organization**. Those properties are nonredundant here: direction carries recurrent cross-species phenotype-space geometry, whereas amount is positively associated with geographic organization. A mean colour cannot distinguish a broad but unstructured distribution from a geographically partitioned polymorphism, and D alone cannot identify which phenotype-space contrast is varying.

### Two simple explanations fail fresh-data tests

The H3 tests sharpen what the species-level phenotype is not trivially reducible to. Validation-cohort D showed no detectable broad tree-wide phylogenetic conservation under any of the three frozen tree placements, while the apparent discovery association with sampled photographic span collapsed essentially to zero in the species-disjoint validation cohort. Together, these out-of-sample results show that reproducible between-species differences in D are not accounted for by either broad shared ancestry as detectable here or the geographic extent over which photographs happened to be sampled.

This inference is deliberately bounded. H3a is a non-support result rather than an equivalence test, so it does not establish a zero phylogenetic effect and does not exclude finer-scale lineage effects, particular clades or repeated evolutionary origins. H3b concerns photographic sampled span, not true biological range size.

### Scope, representativeness and source dependence

The 42,111-species frame gives the analysis broad taxonomic opportunity, but the high-depth cohorts are selected for repeated-observation support and are not a probability sample of global plant diversity. The paper therefore does not estimate the prevalence of flower-colour polymorphism.

Likewise, the prospective confirmation cohort is species-disjoint and prospectively tested, but it comes from the same iNaturalist source/opportunity universe and uses the same measurement system as the earlier cohorts. The fresh-image D transport check also remains within that same source and measurement system. Validation structure should match the intended generalization claim rather than being treated as generically independent (Roberts et al. 2017). The strongest current wording is fresh-image/same-system transport for D and prospective species-disjoint confirmation for the frozen H2 axis, not independent-source replication.

A stronger external validation would apply the same frozen q_white/W estimand and support rules to an independently generated image source, curated field dataset or another measurement system without retuning the axis.

### Conclusion

The methodological advance is a validated route from repeated observations to a comparative within-species trait distribution: measurement reproducibility, location-blind phenotype construction, target discovery and prospective confirmation are treated as separate inferential stages. The ecological advance is that cross-species generality appears in the **structure of intraspecific variation**. Continuous within-species colour displacement shows excess achromatic–chromatic alignment, while species with greater four-state diversity also show stronger geographic organization of that variation even though sampled photographic span does not explain the pattern in validation. Thus **what varies** can be recurrent across species even when **where it is sorted** and the environmental coefficient that sorts it are context dependent. Bounded *S. littorea* molecular and PAL/WAL persistence evidence is consistent with a many-to-one pigment-network accessibility and maintenance filter, whereas BIO5 does not transport as a universal rule and broad tree-wide phylogenetic conservation is not detected. The measured coarse white state remains exposure-coupled, although high-clip exclusion retains excess H2 alignment. More generally, **treating species as distributions rather than mean phenotypes** is necessary but not sufficient: a within-species distribution has empirically separable ecological properties, and **amount, direction and geographic organization need not carry the same information**. In this system, the strongest cross-species regularity lies in phenotype-space direction and in the coupling between diversity amount and spatial organization, not in a universal environmental coefficient.

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

Carlson, J. E., & Holsinger, K. E. (2015). Extrapolating from local ecological processes to genus-wide patterns in colour polymorphism in South African *Protea*. *Proceedings of the Royal Society B: Biological Sciences*, 282, 20150583. https://doi.org/10.1098/rspb.2015.0583

Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2019). Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. *Ecology*, 100(12), e02876. https://doi.org/10.1002/ecy.2876

Casimiro-Soriguer, I., Narbona, E., Buide, M. L., del Valle, J. C., & Whittall, J. B. (2016). Transcriptome and biochemical analysis of a flower color polymorphism in *Silene littorea* (Caryophyllaceae). *Frontiers in Plant Science*, 7, 204. https://doi.org/10.3389/fpls.2016.00204

Dalrymple, R. L., Kemp, D. J., Flores-Moreno, H., Laffan, S. W., White, T. E., Hemmings, F. A., & Moles, A. T. (2020). Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*, 228, 1972–1985. https://doi.org/10.1111/nph.16737

Di Cecco, G. J., Barve, V., Belitz, M. W., Stucky, B. J., Guralnick, R. P., & Hurlbert, A. H. (2021). Observing the Observers: How Participants Contribute Data to iNaturalist and Implications for Biodiversity Science. *BioScience*, 71(11), 1179–1188. https://doi.org/10.1093/biosci/biab093

Laitly, A., Callaghan, C. T., Delhey, K., & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution*, 11, 4071–4083. https://doi.org/10.1002/ece3.7307

Luong, Y., Gasca-Herrera, A., Misiewicz, T. M., & Carter, B. E. (2023). A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences*, 11(5), e11546. https://doi.org/10.1002/aps3.11546

McKenzie, P. F., Church, S. H., & Hopkins, R. (2026). High-Throughput iNaturalist Image Analysis Reveals Flower Color Divergence in *Monarda fistulosa*. *The American Naturalist*, 208(1), 101–109. https://doi.org/10.1086/739413

Narbona, E., Wang, H., Ortiz, P. L., Arista, M., & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology*, 20(Suppl. 1), 8–20. https://doi.org/10.1111/plb.12575

Palacio, F. X., Graco-Roza, C., de Bello, F., & Carmona, C. P. (2025). Integrating intraspecific trait variability in functional diversity: An overview of methods and a guide for ecologists. *Ecological Monographs*, 95(2), e70024. https://doi.org/10.1002/ecm.70024

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, 9, Article 39. https://doi.org/10.2202/1544-6115.1585

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40, 913–929. https://doi.org/10.1111/ecog.02881

Simpson, E. H. (1949). Measurement of Diversity. *Nature*, 163, 688. https://doi.org/10.1038/163688a0

Sapir, Y., Gallagher, M. K., & Senden, E. (2021). What Maintains Flower Colour Variation within Populations? *Trends in Ecology & Evolution*, 36(6), 507–519. https://doi.org/10.1016/j.tree.2021.01.011

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

**Figure 5. Spatial organization accompanies species-level polymorphism while two simple explanations fail fresh-data tests.** (a) Observed D–spatial-organization partial correlations (diamonds) against the mean and 95% interval of the frozen geometry-preserving null for discovery, validation and the validation matched flower-minus-background response. Observed partial rho = 0.1266367 (p = 0.007), 0.0992877 (p = 0.025) and 0.1162411 (p = 0.010), respectively. (b) Validation-cohort Blomberg K under the three frozen phylogenetic placement scenarios: S1 K = 0.0710190, p = 0.2716; S2 K = 0.0601476, p = 0.4134; S3 K = 0.0707577, p = 0.2674. No frozen placement supported detectable broad tree-wide conservation at p < 0.05; this is not an equivalence test. (c) The discovery association with sampled span (rho = 0.1798786, p = 0.00089996) collapsed to essentially zero in the species-disjoint validation cohort (rho = -0.0025855, p = 0.9586021). Sampled photographic span is not true biological range size.

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

