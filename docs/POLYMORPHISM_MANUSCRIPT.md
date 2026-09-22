# Within-species flower-colour variation shows achromatic–chromatic alignment beyond coarse colour-state composition

**Working manuscript draft — 2026-09-18**

**Status:** repository-grounded draft. Numerical claims are restricted to frozen result/claim files in this repository. A bounded literature audit has been integrated; journal-specific formatting remains pending.

## Abstract

Flower colour is commonly summarized as a species mean or categorical state, removing within-species diversity from macroecological analysis. We quantified a continuous four-state diversity phenotype, D, and tested whether continuous within-species colour displacement shows recurrent geometry across plant species. Observer-disjoint reserve partitions recovered stable between-species rankings in D (median Spearman rho = 0.789). Discovery/reserve colour-space analyses localized recurrent geometry to a white-versus-nonwhite contrast, which was then frozen and tested prospectively in a species-disjoint third cohort from the same iNaturalist opportunity universe. The prospective primary test included 158 species and yielded W = 0.517 versus a coarse-state-preserving structured-null median of 0.457 (1.13-fold; p = 0.001). Thus the confirmatory quantity is excess alignment conditional on measured coarse colour-state composition, not the full difference from isotropy. A post-confirmatory gate-reapplied null retained support, whereas a background-white proxy identified unresolved exposure/background-context confounding of the coarse white state. Across the original high-depth cohorts, greater D also remained associated with stronger within-species geographic colour organization after sampled-span, technical-failure, background and ambiguity checks. The evidence supports bounded geometric regularity and species-specific spatial organization while leaving both mechanism and image-formation contributions to the measured white state unresolved.

**Keywords:** flower colour; polymorphism; intraspecific variation; citizen science; reproducibility; colour space; phylogenetic signal

---

## Introduction

Macroecological analyses usually represent species by a single trait value or by a dominant categorical state. That compression is often practical, but it removes the distribution of phenotypes within species. Recent syntheses increasingly treat intraspecific trait variation as an ecological object in its own right rather than only residual variation around a species mean (Palacio et al. 2025). Intraspecific variation can itself be biologically informative, yet it is difficult to scale because apparent diversity may reflect uneven sampling, observer identity or measurement error rather than stable differences among species. Flower colour provides a tractable test case because within-population colour variation is a long-standing evolutionary problem with multiple possible maintaining processes (Sapir et al. 2021; Narbona et al. 2018). Repeated community-science photographs can sample many individuals over broad spatial extents, and prior studies show both the promise and the measurement challenges of deriving floral or organismal colour from such images (Laitly et al. 2021; Luong et al. 2023; McKenzie et al. 2026). The resulting within-species diversity should therefore only be interpreted as a species phenotype if it is reproducible under explicit validation.

We therefore treat flower-colour polymorphism as a continuous quantity rather than forcing species into a binary polymorphic/monomorphic classification. This choice separates our macroecological estimand from the narrower classical definition of within-population discrete colour polymorphism while retaining the broader problem of how within-species colour variation is structured (Sapir et al. 2021; Narbona et al. 2018). For four frozen biological colour states—white, yellow/orange, red/pink and blue/purple—we define species-level diversity as (D = 1-sum_k p_k^2). This formulation preserves gradation in the relative frequencies of colour states while avoiding an arbitrary binary threshold. Our first question is a measurement-validity question: do independent sets of observers recover similar between-species rankings in D?

Amount of variation is not the same as geometry of variation. Two species can have similar D while differing in which colours separate their within-species modes. We therefore ask whether continuous within-species colour displacement is directionless or repeatedly concentrated along a particular colour-space axis. The original discovery and reserve analyses first detected label-free directional concentration and then, through a construction-preserving audit, localized that concentration to a fixed white-versus-equal-nonwhite contrast. Because that named contrast was identified after the original H2 geometry had been opened, those cohorts provide discovery and target-localization evidence rather than an untouched confirmatory test of the named axis.

The decisive H2 test is therefore prospective. Before opening a new biological cohort, we froze the white-versus-nonwhite axis, the W statistic, the admissibility thresholds, the structured null, support gates and one-shot execution contract. We then evaluated that fixed target in a species-disjoint third cohort selected from the same iNaturalist opportunity universe. Separately, we retained a previously frozen species-disjoint test of whether species with greater D also show stronger within-species geographic colour organization. That analysis provides a structural link between the amount of variation and how variation is arranged in space, without assigning a causal mechanism. Finally, we used two frozen H3 tests as alternative-explanation filters: whether species differences in D show broad tree-wide phylogenetic conservation, and whether an apparent association with sampled photographic span survives species-disjoint replication. Together these tests distinguish recurrent geometry and spatial organization from two simple features of ancestry or sampling opportunity.

---

## Methods

### Global sampling frame and inferential cohorts

The project began from the **Repeated Global Flower-Colour Atlas (RGFCA)**, an upstream image-first framework originally designed to test whether within-species flower-colour discontinuities recur in the same broad geographic regions across species. The present paper inherits RGFCA's outcome-blind species discovery, high-depth sampling and measurement infrastructure but changes the estimand from shared global boundary geography to species-level polymorphism amount, colour-space geometry and species-specific spatial organization.

The outcome-blind RGFCA metadata-discovery frame contained **42,111 unique iNaturalist species**. The frame was built before the present polymorphism analyses by combining a one-pass baseline with a cache-resistant repeated global discovery census over an 18 × 9 equal-area grid. Twenty metadata-only V2 rounds made 3,240 fixed cell-level request attempts with zero request errors; V1 observation IDs were excluded from V2, and the deduplicated V1 + V2 species union yielded 42,111 species. No candidate image pixels or flower-colour outcomes were used to define this union. The complete lineage is recorded in `docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md`.

A separate metadata-only capacity scan retained all 42,111 species in the universe while quantifying high-depth photo availability after observer capping. Exactly **4,730 species** had capacity for at least 100 retained photographs (`U100`). This capacity is an observation-process property, not a plant trait or evidence of polymorphism. The 42,111-species frame therefore defines the broad opportunity universe rather than a probability sample for estimating global polymorphism prevalence.

The original high-depth programme measured 100 photographs for each of 1,000 species divided into discovery and reserve source cohorts. After the frozen classifiability rule and a minimum of 40 classifiable photographs per species, the discovery inferential frame contained 369 species and the reserve inferential frame contained 363 species. These cohorts were kept separate for validation and replication.

### Acquisition of the original discovery and reserve high-depth cohorts

The discovery and reserve cohorts inherit the frozen RGFCA iNaturalist acquisition contract rather than a separate flower-polymorphism sampling campaign. Candidate records were required to be iNaturalist Research Grade, species-rank observations with photographs and georeferences, the frozen flowering annotation (term 12, value 13), positional accuracy no worse than 5 km, unobscured/open coordinates and one of five allowed photo licences (CC0, CC BY, CC BY-SA, CC BY-NC or CC BY-NC-SA). Candidate-page and candidate-species selection did not use flower colour, and candidate image pixels remained unopened during acquisition.

Within each species, any one observer contributed at most two retained photographs. Final selection used deterministic geographic maximin sampling to obtain exactly 100 raw photographs per species. The frozen candidate pool therefore contained 1,000 species × 100 photographs. A fixed hash-ranked 500-species measurement budget defined the discovery cohort; the reserve cohort was the species-disjoint complement of the other 500 species from the same frozen candidate pool. Thus discovery and reserve differ in species identity but share the same acquisition rules and raw 100-photo denominator.

The acquisition query did not impose a native-range restriction or an explicit `captive=false` / `wild=true` parameter. Research Grade is therefore treated only as the iNaturalist quality-grade criterion and not as proof that every record represents a native or exclusively wild population. Spatial estimands in this paper refer to the observed community-photograph records.

A later third-cohort H2 test was constructed from an outcome-blind candidate frame that excluded the legacy high-depth species and all P500-selected species. The candidate frame contained 3,230 species. Deterministic hash-based selection froze 500 species before fresh metadata retrieval. Fresh retrieval yielded 499 species with 100 authorized rows each, giving 49,900 rows for the prospective biological execution. No failed species or row was replaced.

The high-depth cohort sizes are therefore hypothesis-specific validation denominators. They are not used as estimates of polymorphism prevalence among the 42,111-species frame.

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

D is interpreted as a continuous within-species colour-diversity phenotype. It is not converted to a global binary polymorphic/monomorphic outcome for the primary analyses.

### H1: observer-disjoint reproducibility

The primary H1 protocol tests whether D is reproducible when observers, rather than photographs, are separated between estimates. Observer identities and observation counts determine the split; morph labels and D do not. This observer-level separation is motivated by evidence that iNaturalist observations carry an observer process, including specialization and heterogeneous contribution patterns (Di Cecco et al. 2021), rather than behaving as exchangeable photographs from a fully specified sampling design.

The first-frozen primary protocol generated 200 observer-disjoint partitions. The reserve decision rule required adequate paired-species support, a median split Spearman correlation of at least 2/3, and a 5th-percentile correlation of at least 0.5. Lin's concordance correlation coefficient (CCC), absolute differences and Spearman-Brown projected reliability were retained as agreement diagnostics.

A later deterministic single-split analysis imposed a stronger rho >= 0.80 criterion. Because that stricter protocol was frozen after the first repeated-partition result had already been opened, it is treated as a deliberately harder stress test rather than as a replacement primary analysis.

### H2: continuous colour geometry in the original cohorts

H2 uses normalized nine-colour palette coordinates

[
[mathrm{white},mathrm{yellow},mathrm{orange},mathrm{red},mathrm{pink},mathrm{magenta},mathrm{purple},mathrm{blue},mathrm{bronze}].
]

Within each eligible species, continuous palette rows are transformed into Hellinger space and partitioned by deterministic unlabeled two-means. The resulting two-mode displacement is converted to a unit direction (u_i). Species must also pass a coarse-state second-mode gate and a continuous minor-cluster gate.

Two admissibility tiers are fixed:

- primary tier: both relevant minor-mode fractions >= 0.10;
- strict sensitivity tier: both fractions >= 0.20.

The original broad geometry was evaluated under both isotropic and construction-preserving references. The decisive construction-preserving structured null permutes normalized nine-colour rows across already selected H2 species within each frozen coarse morph while preserving every species × coarse-morph row count, then refits the same unlabeled two-means construction.

Audit of the original discovery/reserve geometry localized the supported recurrent component to the fixed white-versus-equal-nonwhite unit contrast

[
q_{mathrm{white}} = operatorname{normalize}(1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8).
]

For a set of eligible species, the targeted alignment statistic is

[
W = mean_i (u_i^T q_white)^2.
]

The original cohorts establish discovery/audit evidence because q_white was fixed only after the broad geometry had been opened. Projection of q_white out of the species displacement vectors is used to test whether a residual recurrent hue direction remains.

### Prospective third-cohort H2 confirmation

The third-cohort test was designed specifically to separate target discovery from confirmation. Before biological opening, the following were fixed:

- the 500-species deterministic selection rule;
- fresh metadata retrieval and no-replacement rule;
- the location-blind measurement pipeline;
- minimum 40 classifiable rows per measurement-evaluable species;
- minimum 250 measurement-evaluable species for the support gate;
- q_white;
- W;
- the 0.10 primary and 0.20 strict tiers;
- the construction-preserving structured null;
- 999 null repetitions;
- the no-axis-refit rule;
- one-shot/no-rerun execution;
- durable serialization and read-back validation.

All 256 terminal measurement partitions had to complete before the metadata-colour join and H2 stage could open. The support gate was evaluated before W was calculated.

For each tier, the upper-tail Monte Carlo probability is

[
p = (1 + #(W_null >= W_obs)) / 1000.
]

With 999 randomly generated null worlds, p = 0.001 is the minimum attainable value under this plus-one rule; it should be interpreted as the Monte Carlo resolution of the frozen test rather than an exact exhaustive tail probability (Phipson & Smyth 2010).

The prospective H2 target is supported when the primary tier is evaluable and p < 0.05. The strict 0.20 tier is a pre-specified sensitivity test.

### Post-confirmatory H2 validity diagnostics

After the prospective result had been terminalized, we performed three explicitly post-confirmatory diagnostics that cannot replace or redefine the frozen H2 decision. First, we decomposed the primary W relative to the isotropic expectation and the coarse-state-preserving structured-null baseline, and summarized per-species W contributions according to whether white occurred among the two leading coarse states. Second, to assess the observed-set gate asymmetry, we started from the 185 species passing the coarse-state gate and generated 299 additional structured-null worlds with the continuous minor-cluster threshold reapplied after every refit; the frozen primary seed was retained. Third, we assessed a separate background-bearing sample as an indirect exposure/context proxy by comparing background white fraction between white- and nonwhite-classified images within species. This proxy cannot distinguish digital exposure from biological background differences, framing or mask leakage. The successful third-cohort run did not retain image pixels, masks or background palette fractions, so the previously specified clip_fraction, near_clip_fraction and luminance_q99 highlight diagnostics cannot be reconstructed from its sealed artifact without image reacquisition.

### Complementary test: species-level D and within-species geographic organization

We additionally retain a previously frozen, species-disjoint analysis asking whether species with greater D also show stronger internal geographic colour organization. For each species, all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity. Species-specific geographic organization was then summarized as

`rho_i = Spearman(d_geo_ij, d_colour_ij)`,

where positive values indicate that geographically more distant photographs tend to be more colour-dissimilar. The original RGFCA spatial randomization preserved species identity, coordinates, the complete pairwise geographic geometry and the multiset of colour vectors, while permuting complete colour vectors among photographs within species. Each species therefore had one observed `rho_i` and 999 matched within-species null values.

The threshold-free D–spatial statistic is the Spearman correlation across species between fixed species-level D and observed `rho_i`. For each of the same 999 frozen spatial permutations, D is correlated with the corresponding permuted species-level `rho_i` values; the upper-tail randomization probability is `(1 + # {rho_null >= rho_obs}) / 1000`. No spatial distances, colour distances or permutation assignments were regenerated for the present manuscript.

The primary robustness analysis uses a partial-rank version of the same statistic. Rank(D) and rank(`rho_i`) are separately residualized on an intercept, ranked sampled geographic span and ranked clear ROI/flip technical-failure rate; the reported partial Spearman value is the Pearson correlation of those residuals. The statistic is recalculated for each of the 999 frozen within-species spatial null realizations.

Reserve also provides a matched flower-minus-background response. This is calculated within species as `Spearman(d_geo_ij, d_flower_ij - d_background_ij)`, using paired flower and background colour measurements from the same photographs and joint same-photo permutations. It is not the difference between separate flower and background Spearman coefficients.

Exact D_min4 and D_max4 completions provide uniform endpoint stress tests for ambiguous-palette observations under the frozen four-state model. The reporting-only receipt for the current manuscript reproduces these already frozen Step-8/Step-9 outputs; no model, threshold, distance metric or null was refitted for this submission. These analyses test a structural correlate of D, not the causal process maintaining polymorphism.

### H3a: broad phylogenetic signal

H3a asks whether reserve-cohort D shows broad tree-wide phylogenetic structure. The primary statistic is Blomberg's K, calibrated by 9,999 tip-label permutations on each of three frozen V.PhyloMaker2 placement scenarios (S1, S2 and S3). Reserve is the fresh primary cohort; each scenario retains 341 reserve species.

Support requires raw-D K to have p < 0.05 on all three scenarios. Pagel's lambda and an opportunity-adjusted residual trait are secondary/sensitivity analyses and cannot rescue a failed primary reserve test.

### H3b: sampled photographic span

H3b tests replication of a discovery association between D and sampled photographic span. The primary association is Spearman correlation between D and log-transformed sampled span with a 20,000-permutation two-sided probability. Reserve is the species-disjoint replication cohort. Observer/classifiability-adjusted partial-rank estimates and S1-S3 rank-PGLS models are retained as sensitivities.

The predictor is sampled photographic span under the fixed measurement design; it is not interpreted as true biological range size.

### Reproducibility and frozen decisions

The analysis is governed by frozen protocols, machine-readable result files, input hashes, workflow receipts and no-rescue rules. The third-cohort prospective execution used one authorized run and forbade species replacement, threshold changes, axis refitting and rerun-based selection after biological opening. Image pixels were not persisted.

The P500 prospective expansion is retained separately as a measurement-transport result. Its H2 calculation failed during post-calculation serialization before a durable biological H2 result was written. Under the one-shot contract, P500 supplies neither confirmation nor refutation of H2 and is not replayed for prospective status.

---

## Results

### H1: species-level polymorphism is reproducible across observer-disjoint photo sets

The first-frozen repeated-partition H1 test supported observer-disjoint reproducibility in reserve. Across 200 partitions, the median number of paired species was 329. Median split Spearman rho was **0.7891**, with a 5th percentile of **0.7652** and a 95th percentile of 0.8109. Median Lin CCC was **0.8548**, and median Spearman-Brown projected reliability was **0.8821**.

The later deterministic stress test retained 363 reserve species with zero observer leakage and yielded rho = **0.7927** (bootstrap 95% interval **0.7419–0.8324**) and CCC = **0.8474**. This split missed its deliberately stricter prespecified rho = 0.80 floor by 0.0073. We therefore retain the first-frozen H1 support while explicitly rejecting a claim of near-perfect or split-invariant reliability.

These results admit D as a reproducible high-depth species phenotype for the subsequent geometry analyses, but they do not estimate global polymorphism prevalence.

### H2 discovery and audit: the recurrent component localizes to white versus nonwhite

In the original discovery and reserve cohorts, the fixed white-axis statistic exceeded the construction-preserving structured null at both admissibility tiers.

At the primary 0.10 tier, discovery contained 152 eligible species with W = **0.514625** (null median 0.430808, p = **0.001**) and reserve contained 129 species with W = **0.514586** (null median 0.466546, p = **0.001**).

At the strict 0.20 tier, discovery contained 75 species with W = **0.542355** (null median 0.443626, p = **0.001**) and reserve contained 65 species with W = **0.510517** (null median 0.469943, p = **0.008**).

When the white contrast was projected out, residual directional concentration no longer exceeded the same construction-preserving null. A non-white-only diagnostic likewise did not support a recurrent chromatic hue axis. Thus the legacy H2 signal was narrow: the repeatable component was concentrated along an achromatic–chromatic white-versus-nonwhite direction rather than a general shared hue direction.

Because q_white was identified after the original broad H2 geometry had been opened, these results establish target localization but not untouched prospective confirmation.

### Prospective third-cohort measurement gate

The prospective third-cohort run completed all **256 / 256** terminal partitions and produced **49,900** terminal rows from **499** species. All measurement IDs were unique; no species or rows were replaced; image pixels were not persisted.

Of the 49,900 rows, **25,788** were classifiable into the four biological states and 24,112 were nonclassifiable. Using the prespecified minimum of 40 classifiable rows per species, **377** species were measurement-evaluable, exceeding the required minimum of 250. The measurement-support gate therefore passed before H2 was opened.

### Prospective third-cohort H2: the frozen white axis is confirmed

At the primary 0.10 tier, **158** species produced admissible nonzero displacement vectors. Observed W was **0.517**. Across 999 structured-null worlds, the null median was **0.457** and the 95% interval was **0.436–0.475**. Observed W was 1.13 times the null median, with upper-tail p = **0.001**.

At the strict 0.20 tier, **86** species were eligible. Observed W was **0.533**, compared with a null median of **0.459** and a 95% interval of **0.433–0.487**. Observed W was 1.16 times the null median, again with p = **0.001**.

The frozen terminal verdict was therefore

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

This is an untouched prospective test of a previously frozen axis in a species-disjoint cohort. Because the third cohort was drawn from the same iNaturalist source/opportunity universe and processed with the same measurement system, it is not described as an independent-source replication.

### What the prospective H2 contrast confirms, and post hoc validity diagnostics

The isotropic expectation is 0.125, whereas the frozen primary structured-null median was **0.457** and observed W was **0.517**. The confirmatory result therefore concerns the increment above a coarse-state-preserving construction baseline. Of the 158 primary vector species, **137 (86.7%)** had white as one of their two leading coarse morphs; mean per-species white-axis contribution was 0.562 for those species versus 0.224 for the other 21.

A post-confirmatory diagnostic re-applied the continuous minor-cluster gate in each of 299 structured-null worlds, starting from the 185 coarse-gate species. The null median was **0.447** (95% interval **0.429–0.470**; plus-one p = **0.0033**), so the frozen conditional null was slightly more conservative.

A previously frozen digital-highlight protocol was not executed on the third cohort. The terminal artifact retained neither image pixels nor background palette fractions. In a separate post hoc background-available sample of 461 species, white-classified images had a higher median background white fraction than nonwhite-classified images (0.0692 versus 0.0553); 66.4% of species showed the same within-species direction (paired Wilcoxon p = 4.6 × 10^-12). This is compatible with exposure/scene-brightness or ROI effects but also with biological background differences and does not prove overexposure. Because the structured null conditions on measured coarse-state composition, it cannot falsify an artifact that acts upstream by creating or enriching the white state. Digital exposure/background-context confounding therefore remains unresolved.

### Greater D is associated with stronger within-species geographic colour organization

The positive D–spatial association reproduced across the two species-disjoint high-depth cohorts. The raw association was rho = **0.0892133** (p = **0.034**) in discovery and rho = **0.1016008** (p = **0.025**) in reserve.

After controlling for sampled geographic span and clear ROI/flip technical-failure rate, the geometry-preserving analysis remained positive in discovery (partial rho = **0.1266367**, p = **0.007**) and reserve (partial rho = **0.0992877**, p = **0.025**). In reserve, the matched flower-minus-background response was also positive (partial rho = **0.1162411**, p = **0.010**).

The reserve result remained supported when ambiguous-palette observations were assigned uniformly to the exact diversity-minimizing or diversity-maximizing four-state endpoints. For the primary spatial response, D_min4 gave rho = **0.0970781** (p = **0.029**) and D_max4 gave rho = **0.1252858** (p = **0.008**); for the flower-minus-background response the corresponding p-values were **0.009** and **0.006**. Thus species with greater measured flower-colour diversity tend to show stronger internal geographic organization, but this association does not identify the process that creates or maintains that organization.

### H3a: no detectable broad tree-wide conservation of D

Reserve Blomberg-K tests were nonsignificant under all three frozen tree-placement scenarios:

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Pagel's lambda was small in each scenario and its tests against lambda = 0 were also unsupported. Opportunity-adjusted K sensitivities remained nonsignificant. The frozen verdict was `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

This result closes the tested claim of broad tree-wide signal under the frozen reserve design; it does not imply that phylogeny is irrelevant to flower-colour polymorphism at all evolutionary scales.

### H3b: discovery span effect collapses in reserve

Discovery showed a positive association between D and sampled photographic span (n = 369, rho = **0.1798786**, p = **0.00089996**). The species-disjoint reserve did not reproduce that effect (n = 363, rho = **-0.0025855**, p = **0.9586021**). The observer/classifiability-adjusted partial-rank result was likewise near zero (rho = 0.0055187, p = 0.9162042), and S1-S3 rank-PGLS sensitivities were unsupported.

The frozen verdict was `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

---

## Discussion

### A measurable phenotype before a mechanism

The first result is methodological but biologically consequential: within-species flower-colour diversity can be summarized as a continuous species-level phenotype that is reproducible across completely disjoint observer sets under the first-frozen high-depth validation design. This does not mean that D is measured without error. The stricter deterministic split deliberately exposes that limitation: an observed rho of 0.793 is strong enough to preserve broad species ordering but not strong enough to justify claims of near-perfect or split-invariant reliability.

This distinction matters for macroecological work with citizen-science photographs. Repeated observations can recover more than a modal species colour, but the reliability of the derived distribution should be tested directly rather than assumed from sample size alone. That caution is consistent with direct evaluations of colour information in citizen-science photographs and with broader evidence that observer behaviour is part of the iNaturalist observation process (Laitly et al. 2021; Di Cecco et al. 2021).

The methodological contribution is architectural rather than a claim to a new standalone statistic. Gini–Simpson diversity, rank correlations, Hellinger transformation, two-means clustering, Jensen–Shannon divergence, permutation tests and phylogenetic signal statistics are established tools. What is specific to this study is their assembly around fixed high-depth species sampling, observer-disjoint validation, location-blind image measurement, explicit technical-versus-ambiguous missingness, construction-preserving nulls and a new species/photo-disjoint prospective confirmation cohort. This design treats a within-species phenotype distribution as a species-level comparative trait while keeping measurement validity, target discovery and confirmation as separate inferential stages.

### A recurrent axis, not arbitrary colour-space variation

The strongest positive biological result is geometric. The original discovery/reserve analyses showed that the recurrent construction-controlled component of within-species colour variation was overwhelmingly associated with a white-versus-nonwhite direction. White-versus-pigmented flower-colour combinations have historical precedent in floristic and experimental work, including observations that some anthocyanin-associated white/pigmented combinations are disproportionately represented in particular floras (Warren & Mackenzie 2001), but that precedent does not specify the mechanism of the present axis. Removing that axis eliminated the excess directional concentration, and the remaining non-white geometry did not support a shared hue direction.

The third-cohort result changes the evidential status of this finding, but in a specific way. Species selection, measurement support, q_white, W, thresholds and the structured null were fixed before biological opening. At the primary tier, observed W = 0.517 exceeded a structured-null median of 0.457 (p = 0.001). Because that null already preserves coarse-state composition and is itself far above the isotropic expectation of 0.125, the prospectively confirmed quantity is **excess achromatic–chromatic alignment conditional on the measured coarse colour states**, not the existence of white-versus-nonwhite coarse combinations per se.

### What the achromatic–chromatic axis does not identify

The white-versus-nonwhite geometry is descriptive, not mechanistic. Reviews of flower-colour polymorphism emphasize that pollinator-mediated selection, abiotic selection, drift, gene flow, mating system and pigment genetics can all contribute in different systems (Sapir et al. 2021; Narbona et al. 2018). These analyses do not identify pigment chemistry, whether white states arise by pigment loss or non-white states by pigment gain, or the evolutionary direction of transitions. They also do not distinguish among developmental, genetic, pollinator-mediated or abiotic mechanisms.

The construction-preserving null asks whether continuous within-species geometry adds alignment after the frozen coarse-state composition and its palette mapping are held fixed. It therefore controls a construction baseline but does not validate the origin of the coarse white state itself. An exposure, scene-brightness or segmentation artifact acting upstream on white classification is conditioned on rather than falsified by this null. The post hoc background-white proxy makes that limitation explicit without proving that exposure caused the observed classifications.

### From a repeated global atlas to species-level generality

The present analysis changes the level at which generality is sought. RGFCA was originally constructed to ask whether different species place strong flower-colour discontinuities in the same broad geographic regions. That shared-geography estimand did not provide the positive biological spine retained here. Instead, the same high-depth framework showed that species-specific spatial organization was measurable even when cross-species boundary geography was weak or unsupported. The current paper therefore treats that heterogeneity itself as the object of study: the common rule is sought in the amount and colour-space geometry of polymorphism, while the geographic realization is allowed to remain species-specific. In this sense, the strongest cross-species generality is more evident in **phenotype space than in geographic space**: species repeatedly use a similar achromatic–chromatic direction of variation, but they need not express or maintain that variation in the same places.

### Spatial organization is the strongest ecological clue

The replicated association between D and within-species geographic colour organization provides a positive clue about why species differ in polymorphism. Species with greater D are not merely those sampled across larger geographic extents: the reserve sampled-span association with D collapses to zero, whereas the D–spatial-organization relationship persists after sampled-span and clear technical-failure adjustment and remains positive in a matched flower-minus-background contrast. The ambiguity-endpoint analysis further shows that the association is not tied to one arbitrary treatment of unresolved palette compositions.

This result is structural rather than causal. Stronger geographic organization could arise from spatially varying abiotic selection, turnover in pollinator communities, restricted dispersal or gene flow, demographic history, drift, mating-system differences, or combinations of these processes. The current photographs and occurrence geometry cannot distinguish among them. Together with the prospective H2 result, the evidence suggests a two-stage working model: measured within-species colour variation contains an achromatic–chromatic component stronger than expected from coarse-state composition alone, while species-specific ecological and demographic processes may determine how strongly measured variants are maintained or sorted across space. The first statement is the bounded geometric inference supported by H2; the second remains a mechanistic hypothesis motivated by the replicated D–spatial association. Neither establishes that the measured coarse white state is free of image-exposure or background-context effects.

### Two simple explanations fail fresh-data tests

The H3 tests sharpen what the species-level phenotype is not trivially reducible to. Reserve D showed no detectable broad tree-wide phylogenetic conservation under any of the three frozen tree placements, while the apparent discovery association with sampled photographic span collapsed essentially to zero in the species-disjoint reserve. Together, these out-of-sample results show that reproducible between-species differences in D are not accounted for by either broad shared ancestry as detectable here or the geographic extent over which photographs happened to be sampled.

This inference is deliberately bounded. H3a is a non-support result rather than an equivalence test, so it does not establish a zero phylogenetic effect and does not exclude finer-scale lineage effects, particular clades or repeated evolutionary origins. H3b concerns photographic sampled span, not true biological range size.

### Scope, representativeness and source dependence

The 42,111-species frame gives the analysis broad taxonomic opportunity, but the high-depth cohorts are selected for repeated-observation support and are not a probability sample of global plant diversity. The paper therefore does not estimate the prevalence of flower-colour polymorphism.

Likewise, the third cohort is species-disjoint and prospectively tested, but it comes from the same iNaturalist source/opportunity universe and uses the same measurement system as the earlier cohorts. Validation structure should match the intended generalization claim rather than being treated as generically independent (Roberts et al. 2017). The strongest current wording is prospective species-disjoint confirmation or transport of the frozen axis, not independent-source replication.

A stronger external validation would apply the same frozen q_white/W estimand and support rules to an independently generated image source, curated field dataset or another measurement system without retuning the axis.

### Conclusion

Within-species flower-colour diversity can be measured reproducibly as a continuous species phenotype under high-depth photographic sampling. In a pre-frozen species-disjoint third cohort, continuous colour displacement showed excess alignment with the fixed white-versus-nonwhite axis relative to a coarse-state-preserving structured null. That increment is robust to reapplying the continuous-cluster gate within null worlds, but digital exposure/background-context confounding of the coarse white state remains unresolved. Species with greater D also show stronger within-species geographic colour organization across the original high-depth cohorts after sampled-span, technical-failure, background and ambiguity checks. The evidence supports bounded geometric regularity plus species-specific spatial organization without identifying the ecological, evolutionary or image-formation processes responsible for the measured white states.

---

## Frozen claim boundaries for submission

The manuscript may claim:

1. observer-disjoint reproducibility of continuous four-state D under the first-frozen high-depth validation rule;
2. targeted discovery/audit localization of legacy H2 geometry to white versus nonwhite;
3. untouched prospective species-disjoint confirmation of excess alignment with the already frozen white axis relative to the coarse-state-preserving structured null;
4. replicated positive association between species-level D and within-species geographic colour organization, with background, technical-failure and ambiguity stress tests;
5. non-support for the tested broad phylogenetic-signal claim;
6. non-replication of the sampled-photographic-span association.

The manuscript must not claim:

- global polymorphism prevalence from the high-depth cohorts;
- independent-source replication of H2;
- pigment-loss/gain mechanism or evolutionary direction;
- adaptive causation by pollinators, climate or habitat;
- a recurrent non-white hue axis;
- absence of all phylogenetic structure;
- irrelevance of true geographic range size;
- near-perfect or split-invariant H1 reliability;
- a durable P500 H2 result.

---

## Repository evidence map

- Current claim ledger: `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- Paper architecture: `docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260918.md`
- Figure plan: `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`
- Supporting Information map: `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`
- Literature audit: `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`
- H1 reconciliation: `docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md`
- H2 target freeze: `docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`
- Third-cohort protocol: `docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md`
- Third-cohort result/claim freeze: `docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`
- Post-confirmatory validity diagnostics: `docs/POLYMORPHISM_H2_POSTHOC_VALIDITY_DIAGNOSTICS_20260922.md`
- Third-cohort measurement result: `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`
- Third-cohort H2 result: `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`
- H3a protocol/result manifest: `docs/POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md`, `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`
- H3b result freeze: `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`
- P500 terminal status: `docs/P500_PROSPECTIVE_TERMINAL_POSTMORTEM_20260916.md`

## References

Di Cecco, G. J., Barve, V., Belitz, M. W., Stucky, B. J., Guralnick, R. P., & Hurlbert, A. H. (2021). Observing the Observers: How Participants Contribute Data to iNaturalist and Implications for Biodiversity Science. *BioScience*, 71(11), 1179–1188. https://doi.org/10.1093/biosci/biab093

Laitly, A., Callaghan, C. T., Delhey, K., & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution*, 11, 4071–4083. https://doi.org/10.1002/ece3.7307

Luong, Y., Gasca-Herrera, A., Misiewicz, T. M., & Carter, B. E. (2023). A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences*, 11(5), e11546. https://doi.org/10.1002/aps3.11546

McKenzie, P. F., Church, S. H., & Hopkins, R. (2026). High-Throughput iNaturalist Image Analysis Reveals Flower Color Divergence in *Monarda fistulosa*. *The American Naturalist*, 208(1), 101–109. https://doi.org/10.1086/739413

Narbona, E., Wang, H., Ortiz, P. L., Arista, M., & Imbert, E. (2018). Flower colour polymorphism in the Mediterranean Basin: occurrence, maintenance and implications for speciation. *Plant Biology*, 20(Suppl. 1), 8–20. https://doi.org/10.1111/plb.12575

Palacio, F. X., Graco-Roza, C., de Bello, F., & Carmona, C. P. (2025). Integrating intraspecific trait variability in functional diversity: An overview of methods and a guide for ecologists. *Ecological Monographs*, 95(2), e70024. https://doi.org/10.1002/ecm.70024

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, 9, Article 39. https://doi.org/10.2202/1544-6115.1585

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40, 913–929. https://doi.org/10.1111/ecog.02881

Sapir, Y., Gallagher, M. K., & Senden, E. (2021). What Maintains Flower Colour Variation within Populations? *Trends in Ecology & Evolution*, 36(6), 507–519. https://doi.org/10.1016/j.tree.2021.01.011

Warren, J., & Mackenzie, S. (2001). Why are all colour combinations not equally represented as flower-colour polymorphisms? *New Phytologist*, 151, 237–241. https://doi.org/10.1046/j.1469-8137.2001.00159.x

Literature-use boundaries are frozen in `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`. These references support background and interpretation; they do not alter the repository's machine-readable empirical results or claim ceiling.
