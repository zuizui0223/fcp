# A recurrent achromatic–chromatic axis structures within-species flower-colour polymorphism across plant species

**New Phytologist Full Paper — submission-format working draft**

**Authors:** [AUTHOR LIST TO CONFIRM]

**Affiliations:** [AFFILIATIONS TO INSERT]

**Corresponding author:** [NAME / EMAIL TO INSERT]

**Word counts (current working draft):**
- Summary: 169 words
- Introduction: 507 words
- Materials and Methods: 1,675 words
- Results: 866 words
- Discussion: 846 words
- Main text (Introduction through Discussion): 3,894 words
- Figures: 5
- Tables: 1
- Supporting Information: evidence map + planned supplementary figures/tables

**Keywords (alphabetical):** achromatic–chromatic axis; citizen science; flower colour; intraspecific variation; polymorphism; prospective confirmation

## Summary

- Flower colour is usually reduced to one species state, obscuring within-species diversity. We asked whether flower-colour polymorphism can be measured reproducibly as a species phenotype and whether its variation repeatedly follows a common colour-space direction.
- A continuous four-state diversity score was validated with observer-disjoint photo sets. Continuous nine-colour geometry was audited in discovery/reserve cohorts, then an already frozen white-versus-nonwhite axis was tested prospectively in a pre-frozen species-disjoint third cohort.
- Reserve observer-disjoint partitions gave median Spearman rho = 0.789. The third cohort completed 49,900 rows from 499 species; 377 passed the support gate. The frozen axis was supported at the primary tier (158 species, W = 0.51725, p = 0.001) and strict tier (86 species, W = 0.53293, p = 0.001).
- Within-species flower-colour diversity is reproducibly measurable under high-depth photographic sampling, and its strongest prospectively confirmed recurrent geometry is achromatic–chromatic rather than a general shared hue direction. Fresh reserve tests also make broad ancestry and sampled photographic extent insufficient as simple explanations of species differences in D.

---

## Introduction

Macroecological analyses usually represent species by a single trait value or by a dominant categorical state. That compression is often practical, but it removes the distribution of phenotypes within species. Recent syntheses increasingly treat intraspecific trait variation as an ecological object in its own right rather than only residual variation around a species mean (Palacio et al. 2025). Intraspecific variation can itself be biologically informative, yet it is difficult to scale because apparent diversity may reflect uneven sampling, observer identity or measurement error rather than stable differences among species. Flower colour provides a tractable test case because within-population colour variation is a long-standing evolutionary problem with multiple possible maintaining processes (Sapir et al. 2021; Narbona et al. 2018). Repeated community-science photographs can sample many individuals over broad spatial extents, and prior studies show both the promise and the measurement challenges of deriving floral or organismal colour from such images (Laitly et al. 2021; Luong et al. 2023; McKenzie et al. 2026). The resulting within-species diversity should therefore only be interpreted as a species phenotype if it is reproducible under explicit validation.

We therefore treat flower-colour polymorphism as a continuous quantity rather than forcing species into a binary polymorphic/monomorphic classification. This choice separates our macroecological estimand from the narrower classical definition of within-population discrete colour polymorphism while retaining the broader problem of how within-species colour variation is structured (Sapir et al. 2021; Narbona et al. 2018). For four frozen biological colour states—white, yellow/orange, red/pink and blue/purple—we define species-level diversity as (D = 1-sum_k p_k^2). This formulation preserves gradation in the relative frequencies of colour states while avoiding an arbitrary binary threshold. Our first question is a measurement-validity question: do independent sets of observers recover similar between-species rankings in D?

Amount of variation is not the same as geometry of variation. Two species can have similar D while differing in which colours separate their within-species modes. We therefore ask whether continuous within-species colour displacement is directionless or repeatedly concentrated along a particular colour-space axis. The original discovery and reserve analyses first detected label-free directional concentration and then, through a construction-preserving audit, localized that concentration to a fixed white-versus-equal-nonwhite contrast. Because that named contrast was identified after the original H2 geometry had been opened, those cohorts provide discovery and target-localization evidence rather than an untouched confirmatory test of the named axis.

The decisive H2 test is therefore prospective. Before opening a new biological cohort, we froze the white-versus-nonwhite axis, the W statistic, the admissibility thresholds, the structured null, support gates and one-shot execution contract. We then evaluated that fixed target in a species-disjoint third cohort selected from the same iNaturalist opportunity universe. Finally, we used two frozen H3 tests as alternative-explanation filters: whether species differences in D show broad tree-wide phylogenetic conservation, and whether an apparent association with sampled photographic span survives species-disjoint replication. These tests ask whether the measurable phenotype can be reduced to two simple features of ancestry or sampling opportunity; they do not search post hoc for replacement predictors.

---

---

## Materials and Methods

### Global sampling frame and inferential cohorts

The project began from a global flower-colour sampling frame containing 42,111 species. This frame defines the broad opportunity universe, not a probability sample for estimating the global prevalence of polymorphism.

The original high-depth programme measured 100 photographs for each of 1,000 species divided into discovery and reserve source cohorts. After the frozen classifiability rule and a minimum of 40 classifiable photographs per species, the discovery inferential frame contained 369 species and the reserve inferential frame contained 363 species. These cohorts were kept separate for validation and replication.

A later third-cohort H2 test was constructed from an outcome-blind candidate frame that excluded the legacy high-depth species and all P500-selected species. The candidate frame contained 3,230 species. Deterministic hash-based selection froze 500 species before fresh metadata retrieval. Fresh retrieval yielded 499 species with 100 authorized rows each, giving 49,900 rows for the prospective biological execution. No failed species or row was replaced.

The high-depth cohort sizes are therefore hypothesis-specific validation denominators. They are not used as estimates of polymorphism prevalence among the 42,111-species frame.

### Photographic measurement and outcome firewall

All high-depth image rows were processed with a frozen location-blind measurement machine rather than by manual selection after outcome inspection. For the prospective third cohort, the inherited machine was pinned to source commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`. It used the unchanged ROI-v4 flower-detection step, EfficientSAM segmentation, a normalized nine-colour flower palette (white, yellow, orange, red, pink, magenta, purple, blue and bronze), the four coarse biological states, and the same ROI, flip and palette-quality gates used in the qualified prospective design. Rows failing these gates were assigned terminal nonclassifiable states; failure never triggered a replacement image or species.

The third-cohort metadata draw was row-disjoint from 178,462 previously used observation IDs and 178,462 previously used photo IDs. Before any candidate pixel was opened, the selected-species manifest and authorized metadata were checksum-verified. Each photo ID was converted deterministically to a blind measurement ID. Measurement workers could access only the measurement ID, image filename and photo licence; species identity, taxon ID, coordinates, observer ID, observation ID and prospective rank remained sealed until all terminal measurements were complete. The 49,900 rows were distributed across two blind batches, 32 semantic shards per batch and four compute partitions per shard, giving 256 terminal partitions with no early stopping. Image pixels and flower masks were not persisted after partition sealing.

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

The first-frozen primary protocol generated 200 observer-disjoint partitions. Within each species and partition, observers were kept intact and ordered by the number of measured rows they contributed. Ties were resolved by a frozen SHA256 ordering of seed, species and observer ID. Observers were then greedily assigned to the currently smaller half by all-row count, so balancing used observer identity and sampling effort but never colour labels, classifiability, D, geography or H2/H3 outcomes. A species contributed to a partition-level estimate only when each half retained at least 20 classifiable rows; a >=15-per-half analysis was prespecified as sensitivity only.

For each partition, D was calculated separately in the two observer-disjoint halves and species were compared using Spearman correlation. The reserve decision rule required a median of at least 100 paired species, median split Spearman rho >=2/3 and 5th-percentile rho >=0.5. Lin's concordance correlation coefficient (CCC), absolute differences and Spearman-Brown projected reliability were retained as agreement diagnostics rather than substitute decision statistics.

A later deterministic single-split analysis imposed a stronger rho >= 0.80 criterion. Because that stricter protocol was frozen after the first repeated-partition result had already been opened, it is treated as a deliberately harder stress test rather than as a replacement primary analysis.

### H2: continuous colour geometry in the original cohorts

H2 uses normalized nine-colour palette coordinates

[
[mathrm{white},mathrm{yellow},mathrm{orange},mathrm{red},mathrm{pink},mathrm{magenta},mathrm{purple},mathrm{blue},mathrm{bronze}].
]

Within each eligible species, each nine-colour composition was first normalized to unit row sum and transformed to Hellinger coordinates by element-wise square root. Deterministic unlabeled two-means was fitted without using the four coarse morph labels. Mean normalized compositions were then calculated for the larger and smaller continuous clusters, and their difference defined the species displacement vector Delta_i. Because the cluster labels are arbitrary, H2 uses the sign-invariant unit axis u_i = Delta_i / ||Delta_i||. Species must also pass a coarse-state second-mode gate and a continuous minor-cluster gate.

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

The fresh metadata stage made exactly one request per frozen selected species under the prespecified acquisition rule. Of 500 selected species, 499 supplied 100 fresh authorized rows and one supplied 99; the latter was excluded before pixel opening with no replacement and no biological interpretation. The resulting pre-pixel denominator was therefore fixed at 499 species and 49,900 rows.

All 256 terminal measurement partitions had to complete before the metadata-colour join and H2 stage could open. Every frozen row had to receive exactly one terminal state. Only rows classified into the four biological coarse states counted toward measurement support; a species was measurement-evaluable at n_classifiable >=40, and at least 250 measurement-evaluable species were required before the H2 statistic could be opened. A failure of this gate was prespecified as underidentification, not evidence against polymorphism or against the white axis.

The prospective run was also required to survive durable terminalization. Calculation of W in memory was not sufficient: the result had to be serialized, read back from disk, validated as an H2_COMPLETE stage, packaged as an artifact and committed immutably. This requirement was qualified synthetically before biological opening to prevent a post-calculation serialization failure from being misclassified as a biological result.

For each tier, the upper-tail Monte Carlo probability is

[
p = (1 + #(W_null >= W_obs)) / 1000.
]

With 999 randomly generated null worlds, p = 0.001 is the minimum attainable value under this plus-one rule; it should be interpreted as the Monte Carlo resolution of the frozen test rather than an exact exhaustive tail probability (Phipson & Smyth 2010).

The prospective H2 target is supported when the primary tier is evaluable and p < 0.05. The strict 0.20 tier is a pre-specified sensitivity test.

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

---

## Table 1. Inferential roles of the high-depth cohorts

| Cohort / execution | Primary role | Key denominator | H2 status |
|---|---|---:|---|
| Discovery | Measurement calibration and H2 target discovery/audit | 369 D-eligible species | Legacy post-audit localization |
| Reserve | Fresh H1/H3 replication and legacy H2 validation | 363 D-eligible species | Legacy post-audit validation |
| P500 | Prospective measurement-pipeline transport | 499 species; 49,900 rows; 373 measurement-evaluable | No durable biological H2 verdict |
| Third cohort | Untouched prospective test of already frozen q_white/W | 499 species; 49,900 rows; 377 measurement-evaluable | H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED |

The third cohort is species-disjoint but remains within the same iNaturalist source/opportunity universe and measurement system; it is not an independent-source replication.

---

## Results

### H1: species-level polymorphism is reproducible across observer-disjoint photo sets

The first-frozen repeated-partition H1 test supported observer-disjoint reproducibility in reserve. Across 200 partitions, the median number of paired species was 329. Median split Spearman rho was **0.7891**, with a 5th percentile of **0.7652** and a 95th percentile of 0.8109. Median Lin CCC was **0.8548**, and median Spearman-Brown projected reliability was **0.8821**.

The later deterministic stress test retained 363 reserve species with zero observer leakage and yielded rho = **0.7927** (bootstrap 95% interval **0.7419–0.8324**) and CCC = **0.8474**. This split missed its deliberately stricter prespecified rho = 0.80 floor by 0.0073. We therefore retain the first-frozen H1 support while explicitly rejecting a claim of near-perfect or split-invariant reliability.

These results admit D as a reproducible high-depth species phenotype for the subsequent geometry analyses, but they do not estimate global polymorphism prevalence.

### H2 discovery and audit: recurrent geometry localizes to white versus nonwhite

The original H2 analysis first established recurrent label-free geometry before naming its biological direction. At the primary 0.10 tier, discovery had 152 vector species with leading-axis concentration lambda1 = **0.541412**; reserve had 129 vector species with lambda1 = **0.535045**. Reserve mean squared projection onto the frozen discovery axis was **0.524157**, and the independently fitted discovery and reserve leading axes had absolute alignment of approximately **0.986**. Discovery concentration and reserve frozen-axis transport each exceeded both the isotropic reference and the subsequent coarse-state-preserving structured null (structured-null p = **0.001** in both cases). The strict 0.20 tier showed the same direction of support.

Audit then localized that recurrent component to white versus nonwhite. In the original discovery and reserve cohorts, the fixed white-axis statistic exceeded the construction-preserving structured null at both admissibility tiers.

At the primary 0.10 tier, discovery contained 152 eligible species with W = **0.514625** (null median 0.430808, p = **0.001**) and reserve contained 129 species with W = **0.514586** (null median 0.466546, p = **0.001**).

At the strict 0.20 tier, discovery contained 75 species with W = **0.542355** (null median 0.443626, p = **0.001**) and reserve contained 65 species with W = **0.510517** (null median 0.469943, p = **0.008**).

When the white contrast was projected out, residual directional concentration no longer exceeded the same construction-preserving null. A non-white-only diagnostic likewise did not support a recurrent chromatic hue axis. Thus the legacy H2 signal was narrow: the repeatable component was concentrated along an achromatic–chromatic white-versus-nonwhite direction rather than a general shared hue direction.

Because q_white was identified after the original broad H2 geometry had been opened, these results establish target localization but not untouched prospective confirmation.

### Prospective third-cohort measurement gate

The prospective third-cohort run completed all **256 / 256** terminal partitions and produced **49,900** terminal rows from **499** species. All measurement IDs were unique; no species or rows were replaced; image pixels were not persisted.

Of the 49,900 rows, **25,788** were classifiable into the four biological states and 24,112 were nonclassifiable. Using the prespecified minimum of 40 classifiable rows per species, **377** species were measurement-evaluable, exceeding the required minimum of 250. The measurement-support gate therefore passed before H2 was opened.

### Prospective third-cohort H2: the frozen white axis is confirmed

At the primary 0.10 tier, **158** species produced admissible nonzero displacement vectors. Observed W was **0.5172457461**. Across 999 structured-null worlds, the null median was **0.4571428150** and the 95% interval was **0.4358491120–0.4752987776**. Observed W was 1.1315 times the null median, with upper-tail p = **0.001**.

At the strict 0.20 tier, **86** species were eligible. Observed W was **0.5329282123**, compared with a null median of **0.4593196659** and a 95% interval of **0.4328679570–0.4867224043**. Observed W was 1.1603 times the null median, again with p = **0.001**.

The frozen terminal verdict was therefore

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

This is an untouched prospective test of a previously frozen axis in a species-disjoint cohort. Because the third cohort was drawn from the same iNaturalist source/opportunity universe and processed with the same measurement system, it is not described as an independent-source replication.

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

---

## Discussion

### A measurable phenotype before a mechanism

The first result is methodological but biologically consequential: within-species flower-colour diversity can be summarized as a continuous species-level phenotype that is reproducible across completely disjoint observer sets under the first-frozen high-depth validation design. This does not mean that D is measured without error. The stricter deterministic split deliberately exposes that limitation: an observed rho of 0.793 is strong enough to preserve broad species ordering but not strong enough to justify claims of near-perfect or split-invariant reliability.

This distinction matters for macroecological work with citizen-science photographs. Repeated observations can recover more than a modal species colour, but the reliability of the derived distribution should be tested directly rather than assumed from sample size alone. That caution is consistent with direct evaluations of colour information in citizen-science photographs and with broader evidence that observer behaviour is part of the iNaturalist observation process (Laitly et al. 2021; Di Cecco et al. 2021).

### A recurrent axis, not arbitrary colour-space variation

The strongest positive biological result is geometric. The original discovery/reserve analyses showed that the recurrent construction-controlled component of within-species colour variation was overwhelmingly associated with a white-versus-nonwhite direction. White-versus-pigmented flower-colour combinations have historical precedent in floristic and experimental work, including observations that some anthocyanin-associated white/pigmented combinations are disproportionately represented in particular floras (Warren & Mackenzie 2001), but that precedent does not specify the mechanism of the present axis. Removing that axis eliminated the excess directional concentration, and the remaining non-white geometry did not support a shared hue direction.

The third-cohort result changes the evidential status of this finding. The white-axis target was no longer chosen after looking at the new cohort: the species selection, measurement support rule, q_white, W, thresholds and structured null were fixed before biological opening. The primary and strict tests both returned p = 0.001, with observed W well above their null distributions. The axis is therefore not merely a retrospective description of the first two cohorts; it transported prospectively to a new species-disjoint cohort from the same opportunity universe.

### What the achromatic–chromatic axis does not identify

The white-versus-nonwhite geometry is descriptive, not mechanistic. Reviews of flower-colour polymorphism emphasize that pollinator-mediated selection, abiotic selection, drift, gene flow, mating system and pigment genetics can all contribute in different systems (Sapir et al. 2021; Narbona et al. 2018). These analyses do not identify pigment chemistry, whether white states arise by pigment loss or non-white states by pigment gain, or the evolutionary direction of transitions. They also do not distinguish among developmental, genetic, pollinator-mediated or abiotic mechanisms.

The construction-preserving null strengthens the claim that the observed alignment is not explained simply by the frozen coarse-state composition and global mapping from coarse states to the nine-colour palette. It does not convert geometric alignment into a causal mechanism.

### Two simple explanations fail fresh-data tests

The H3 tests sharpen what the species-level phenotype is not trivially reducible to. Reserve D showed no detectable broad tree-wide phylogenetic conservation under any of the three frozen tree placements, while the apparent discovery association with sampled photographic span collapsed essentially to zero in the species-disjoint reserve. Together, these out-of-sample results show that reproducible between-species differences in D are not accounted for by either broad shared ancestry as detectable here or the geographic extent over which photographs happened to be sampled.

This inference is deliberately bounded. H3a is a non-support result rather than an equivalence test, so it does not establish a zero phylogenetic effect and does not exclude finer-scale lineage effects, particular clades or repeated evolutionary origins. H3b concerns photographic sampled span, not true biological range size.

### Scope, representativeness and source dependence

The 42,111-species frame gives the analysis broad taxonomic opportunity, but the high-depth cohorts are selected for repeated-observation support and are not a probability sample of global plant diversity. The paper therefore does not estimate the prevalence of flower-colour polymorphism.

Likewise, the third cohort is species-disjoint and prospectively tested, but it comes from the same iNaturalist source/opportunity universe and uses the same measurement system as the earlier cohorts. Validation structure should match the intended generalization claim rather than being treated as generically independent (Roberts et al. 2017). The strongest current wording is prospective species-disjoint confirmation or transport of the frozen axis, not independent-source replication.

A stronger external validation would apply the same frozen q_white/W estimand and support rules to an independently generated image source, curated field dataset or another measurement system without retuning the axis.

### Conclusion

Within-species flower-colour diversity can be measured reproducibly as a continuous species phenotype under high-depth photographic sampling. Across the original discovery/reserve analyses, the strongest recurrent colour-space component localized to a white-versus-nonwhite axis. A pre-frozen species-disjoint third cohort then prospectively confirmed that same axis at both primary and strict admissibility tiers. Fresh reserve tests further showed no detectable broad tree-wide conservation of D and reduced the apparent discovery association with sampled photographic span to essentially zero. The current evidence therefore isolates a robust geometric regularity while making broad ancestry and sampled extent insufficient as simple explanations under the tested designs; its evolutionary and ecological mechanisms remain open.

---

---

## Acknowledgements

[TO COMPLETE BEFORE SUBMISSION: funding, institutional support, data-provider acknowledgements, and individual contributions that do not meet authorship criteria.]

## Competing interests

[TO CONFIRM BEFORE SUBMISSION.]

## Author contributions

[TO COMPLETE AFTER FINAL AUTHOR LIST.]

## Data availability

Frozen protocols, analysis code, machine-readable results, claim ledgers and canonical figures are versioned in the `zuizui0223/fcp` GitHub repository. The third-cohort prospective result is preserved as an immutable repository commit and GitHub Actions artifact. A permanent archival DOI/version should be added before submission.

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

## Figure legends

**Figure 1. From a global sampling frame to a measurable species-level flower-colour polymorphism phenotype.** (a) Descriptive distributions of D = 1 - sum_k p_k^2 in the discovery (n = 369) and reserve (n = 363) high-depth inferential cohorts; dashed lines mark cohort medians. These validation cohorts are not a prevalence sample. (b) Sampling architecture from the 42,111-species opportunity frame. The original validation lane proceeds to the 500-species discovery and reserve source cohorts and their >=40-classifiable D inferential sets, whereas the prospective-confirmation lane separately proceeds through pre-frozen third-cohort selection/fresh metadata to 499 species × 100 authorized rows and 377 measurement-evaluable species with zero replacements.

**Figure 2. Observer-disjoint reproducibility of the continuous polymorphism score D.** (a) Distributional summary of Spearman split-half correlations across 200 first-frozen observer-disjoint partitions in discovery and reserve, shown as 5th percentile–median–95th percentile intervals. Reserve median rho = 0.7891 and q05 = 0.7652; the dashed line marks the primary median floor of 2/3. (b) Later deterministic observer-disjoint stress tests in discovery and reserve with bootstrap 95% intervals. The reserve estimate was rho = 0.7927 (95% CI 0.7419–0.8324; CCC = 0.8474), narrowly below the prespecified 0.80 floor; this later stress test constrains but does not overwrite the chronologically earlier primary H1 result.

**Figure 3. Discovery and audit of the recurrent white-versus-nonwhite colour-space target in the original cohorts.** (a) Loadings of the fixed zero-sum q_white contrast: white is opposed to the equal mean of the eight non-white palette coordinates. The panel explicitly records that this named axis was isolated only after the original broad H2 geometry had been opened. (b) Legacy targeted W values (diamonds) against the median (points) and 95% interval (bars) of the construction-preserving structured null. Primary 0.10: discovery N = 152, W = 0.514625, p = 0.001; reserve N = 129, W = 0.514586, p = 0.001. Strict 0.20: discovery N = 75, W = 0.542355, p = 0.001; reserve N = 65, W = 0.510517, p = 0.008. Projection-removal and non-white-only falsification diagnostics are reported in Supporting Information.

**Figure 4. Prospective species-disjoint confirmation of the frozen white-versus-nonwhite axis.** Structured-null W distributions from 999 frozen null worlds; dashed lines show null medians and solid vertical lines show observed W. (a) Primary 0.10 tier: N = 158, observed W = 0.5172457461, null median = 0.4571428150, 95% interval 0.4358491120–0.4752987776, upper-tail p = 0.001. (b) Strict 0.20 sensitivity: N = 86, observed W = 0.5329282123, null median = 0.4593196659, 95% interval 0.4328679570–0.4867224043, p = 0.001. The prospective cohort completed 49,900 terminal rows from 499 species, with 377 measurement-evaluable species and zero replacements before H2 opening; it is species-disjoint within the same iNaturalist opportunity universe, not an independent-source replication.

**Figure 5. Broad ancestry and sampled photographic extent fail as simple explanations of species differences in D.** (a) Reserve Blomberg K under the three frozen phylogenetic placement scenarios: S1 K = 0.0710190, p = 0.2716; S2 K = 0.0601476, p = 0.4134; S3 K = 0.0707577, p = 0.2674. No frozen placement supported detectable broad tree-wide conservation at the prespecified p < 0.05 criterion. (b) The discovery association with sampled span (rho = 0.1798786, p = 0.00089996) collapsed to essentially zero in the species-disjoint reserve (rho = -0.0025855, p = 0.9586021). H3a is not an equivalence test, and sampled photographic span is not true biological range size; the tests therefore reject neither all phylogenetic effects nor all geographic effects, but they make these two broad explanations insufficient under the frozen designs.

## Supporting Information

The complete evidence map is provided in `docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`.

**Fig. S1.** Full H1 reserve partition diagnostics for all 200 observer-disjoint partitions, including paired-species counts, Spearman rho, CCC, signed bias and absolute D differences.

**Fig. S2.** H1 discovery-cohort concordance together with the later deterministic reserve stress-test diagnostics and its prespecified rho = 0.80 floor.

**Fig. S3.** Broad pre-target H2 geometry in the original cohorts, including leading-axis concentration, discovery-axis reserve transport and discovery–reserve axis alignment.

**Fig. S4.** Construction-preserving H2 null audit showing the quantities preserved when normalized nine-colour rows are permuted within coarse morph while species × coarse-morph row counts remain fixed.

**Fig. S5.** Residual H2 tests after projecting out q_white together with the low-sample non-white-only diagnostics.

**Fig. S6.** Third-cohort chain of custody from deterministic species selection and fresh-metadata freeze through 256 terminal measurement receipts, support-gate completion and the durable H2_COMPLETE result.

**Fig. S7.** P500 execution chronology showing successful measurement support but no durable biological H2 verdict because post-calculation serialization failed before terminal result persistence.

**Fig. S8.** H3a sensitivity analyses across S1–S3 phylogenies for raw D, finite-sample sensitivity, opportunity-adjusted residuals, Blomberg K and Pagel lambda.

**Fig. S9.** H3b sampled-span sensitivities including raw, finite-sample, observer/classifiability-adjusted partial-rank and rank-PGLS reserve analyses.

