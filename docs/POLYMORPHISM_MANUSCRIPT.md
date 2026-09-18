# A recurrent achromatic–chromatic axis structures within-species flower-colour polymorphism across plant species

**Working manuscript draft — 2026-09-18**

**Status:** repository-grounded draft. Numerical claims are restricted to frozen result/claim files in this repository. Literature citations and journal formatting are not yet integrated.

## Abstract

Flower colour is commonly summarized as a species mean or categorical state, which removes within-species colour diversity from macroecological analyses. We asked whether flower-colour polymorphism can instead be treated as a reproducible species-level phenotype and, if so, whether within-species variation repeatedly follows a common direction in colour space. We began from a 42,111-species flower-colour sampling frame and used species-level high-depth photographic cohorts for hypothesis-specific validation rather than prevalence estimation. Polymorphism was quantified as a continuous four-state diversity score, (D = 1-sum_k p_k^2), over white, yellow/orange, red/pink and blue/purple states. In the species-disjoint reserve cohort, repeated observer-disjoint partitions recovered stable species rankings in D (median Spearman rho = 0.789, 5th percentile = 0.765; median Lin CCC = 0.855), although a later deliberately stricter deterministic split yielded rho = 0.793 and missed a prespecified 0.80 floor. Continuous colour-space analyses of the original discovery and reserve cohorts identified a recurrent white-versus-nonwhite axis beyond a coarse-state-preserving construction null, but that named axis was isolated after the original geometry had been opened. We therefore tested the already frozen axis prospectively in a pre-frozen species-disjoint third cohort drawn from the same iNaturalist opportunity universe. The measurement pipeline completed 49,900 rows from 499 species without replacement; 377 species passed the predeclared measurement-support gate. The fixed white-axis statistic was supported at both the primary 0.10 tier (158 species, W = 0.51725, structured-null p = 0.001) and the stricter 0.20 tier (86 species, W = 0.53293, p = 0.001). By contrast, broad phylogenetic signal in D was unsupported across three frozen reserve phylogenies, and a discovery association with sampled photographic span failed species-disjoint reserve replication. Within-species flower-colour diversity is therefore measurable as a species phenotype under high-information photographic sampling, and its strongest replicated geometric regularity in these data is a recurrent achromatic–chromatic axis rather than a general shared hue direction.

**Keywords:** flower colour; polymorphism; intraspecific variation; citizen science; reproducibility; colour space; phylogenetic signal

---

## Introduction

Macroecological analyses usually represent species by a single trait value or by a dominant categorical state. That compression is often practical, but it removes the distribution of phenotypes within species. Intraspecific variation can itself be biologically informative, yet it is difficult to scale because apparent diversity may reflect uneven sampling, observer identity or measurement error rather than stable differences among species. Flower colour provides a tractable test case: repeated photographs can sample many individuals over broad spatial extents, but the resulting within-species diversity should only be interpreted as a species phenotype if it is reproducible under observer-disjoint resampling.

We therefore treat flower-colour polymorphism as a continuous quantity rather than forcing species into a binary polymorphic/monomorphic classification. For four frozen biological colour states—white, yellow/orange, red/pink and blue/purple—we define species-level diversity as (D = 1-sum_k p_k^2). This formulation preserves gradation in the relative frequencies of colour states while avoiding an arbitrary binary threshold. Our first question is a measurement-validity question: do independent sets of observers recover similar between-species rankings in D?

Amount of variation is not the same as geometry of variation. Two species can have similar D while differing in which colours separate their within-species modes. We therefore ask whether continuous within-species colour displacement is directionless or repeatedly concentrated along a particular colour-space axis. The original discovery and reserve analyses first detected label-free directional concentration and then, through a construction-preserving audit, localized that concentration to a fixed white-versus-equal-nonwhite contrast. Because that named contrast was identified after the original H2 geometry had been opened, those cohorts provide discovery and target-localization evidence rather than an untouched confirmatory test of the named axis.

The decisive H2 test is therefore prospective. Before opening a new biological cohort, we froze the white-versus-nonwhite axis, the W statistic, the admissibility thresholds, the structured null, support gates and one-shot execution contract. We then evaluated that fixed target in a species-disjoint third cohort selected from the same iNaturalist opportunity universe. Finally, we asked whether the resulting species-level polymorphism phenotype shows broad phylogenetic structure or a replicable association with sampled photographic span. These latter analyses are bounded explanatory tests: failure to support them does not negate the measurement or geometry results, and they are not used to search post hoc for alternative predictors.

---

## Methods

### Global sampling frame and inferential cohorts

The project began from a global flower-colour sampling frame containing 42,111 species. This frame defines the broad opportunity universe, not a probability sample for estimating the global prevalence of polymorphism.

The original high-depth programme measured 100 photographs for each of 1,000 species divided into discovery and reserve source cohorts. After the frozen classifiability rule and a minimum of 40 classifiable photographs per species, the discovery inferential frame contained 369 species and the reserve inferential frame contained 363 species. These cohorts were kept separate for validation and replication.

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

The primary H1 protocol tests whether D is reproducible when observers, rather than photographs, are separated between estimates. Observer identities and observation counts determine the split; morph labels and D do not.

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
W = operatorname{mean}_i (u_i^	op q_{mathrm{white}})^2.
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
p = rac{1+#(W_{mathrm{null}}ge W_{mathrm{obs}})}{1000}.
]

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

At the primary 0.10 tier, **158** species produced admissible nonzero displacement vectors. Observed W was **0.5172457461**. Across 999 structured-null worlds, the null median was **0.4571428150** and the 95% interval was **0.4358491120–0.4752987776**. Observed W was 1.1315 times the null median, with upper-tail p = **0.001**.

At the strict 0.20 tier, **86** species were eligible. Observed W was **0.5329282123**, compared with a null median of **0.4593196659** and a 95% interval of **0.4328679570–0.4867224043**. Observed W was 1.1603 times the null median, again with p = **0.001**.

The frozen terminal verdict was therefore

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

This is an untouched prospective test of a previously frozen axis in a species-disjoint cohort. Because the third cohort was drawn from the same iNaturalist source/opportunity universe and processed with the same measurement system, it is not described as an independent-source replication.

### H3a: broad phylogenetic signal is not supported

Reserve Blomberg-K tests were nonsignificant under all three frozen tree-placement scenarios:

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Pagel's lambda was small in each scenario and its tests against lambda = 0 were also unsupported. Opportunity-adjusted K sensitivities remained nonsignificant. The frozen verdict was `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

This result closes the tested claim of broad tree-wide signal under the frozen reserve design; it does not imply that phylogeny is irrelevant to flower-colour polymorphism at all evolutionary scales.

### H3b: sampled-span association fails reserve replication

Discovery showed a positive association between D and sampled photographic span (n = 369, rho = **0.1798786**, p = **0.00089996**). The species-disjoint reserve did not reproduce that effect (n = 363, rho = **-0.0025855**, p = **0.9586021**). The observer/classifiability-adjusted partial-rank result was likewise near zero (rho = 0.0055187, p = 0.9162042), and S1-S3 rank-PGLS sensitivities were unsupported.

The frozen verdict was `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

---

## Discussion

### A measurable phenotype before a mechanism

The first result is methodological but biologically consequential: within-species flower-colour diversity can be summarized as a continuous species-level phenotype that is reproducible across completely disjoint observer sets under the first-frozen high-depth validation design. This does not mean that D is measured without error. The stricter deterministic split deliberately exposes that limitation: an observed rho of 0.793 is strong enough to preserve broad species ordering but not strong enough to justify claims of near-perfect or split-invariant reliability.

This distinction matters for macroecological work with citizen-science photographs. Repeated observations can recover more than a modal species colour, but the reliability of the derived distribution should be tested directly rather than assumed from sample size alone.

### A recurrent axis, not arbitrary colour-space variation

The strongest positive biological result is geometric. The original discovery/reserve analyses showed that the recurrent construction-controlled component of within-species colour variation was overwhelmingly associated with a white-versus-nonwhite direction. Removing that axis eliminated the excess directional concentration, and the remaining non-white geometry did not support a shared hue direction.

The third-cohort result changes the evidential status of this finding. The white-axis target was no longer chosen after looking at the new cohort: the species selection, measurement support rule, q_white, W, thresholds and structured null were fixed before biological opening. The primary and strict tests both returned p = 0.001, with observed W well above their null distributions. The axis is therefore not merely a retrospective description of the first two cohorts; it transported prospectively to a new species-disjoint cohort from the same opportunity universe.

### What the achromatic–chromatic axis does not identify

The white-versus-nonwhite geometry is descriptive, not mechanistic. These analyses do not identify pigment chemistry, whether white states arise by pigment loss or non-white states by pigment gain, or the evolutionary direction of transitions. They also do not distinguish among developmental, genetic, pollinator-mediated or abiotic mechanisms.

The construction-preserving null strengthens the claim that the observed alignment is not explained simply by the frozen coarse-state composition and global mapping from coarse states to the nine-colour palette. It does not convert geometric alignment into a causal mechanism.

### Stable geometry without a broad predictor

The H3 results provide useful limits. Broad reserve phylogenetic signal was unsupported under all three frozen tree placements, and the sampled-span association selected in discovery collapsed essentially to zero in reserve. A stable measurable phenotype and a prospectively confirmed recurrent geometric axis therefore do not require the broad predictors tested here to show corresponding replication.

These negative tests should not be generalized beyond their estimands. The phylogenetic result does not exclude finer-scale lineage effects, particular clades or repeated evolutionary origins. The span result concerns photographic sampled span, not true biological range size.

### Scope, representativeness and source dependence

The 42,111-species frame gives the analysis broad taxonomic opportunity, but the high-depth cohorts are selected for repeated-observation support and are not a probability sample of global plant diversity. The paper therefore does not estimate the prevalence of flower-colour polymorphism.

Likewise, the third cohort is species-disjoint and prospectively tested, but it comes from the same iNaturalist source/opportunity universe and uses the same measurement system as the earlier cohorts. The strongest current wording is prospective species-disjoint confirmation or transport of the frozen axis, not independent-source replication.

A stronger external validation would apply the same frozen q_white/W estimand and support rules to an independently generated image source, curated field dataset or another measurement system without retuning the axis.

### Conclusion

Within-species flower-colour diversity can be measured reproducibly as a continuous species phenotype under high-depth photographic sampling. Across the original discovery/reserve analyses, the strongest recurrent colour-space component localized to a white-versus-nonwhite axis. A pre-frozen species-disjoint third cohort then prospectively confirmed that same axis at both primary and strict admissibility tiers. In contrast, broad phylogenetic signal and a discovery association with sampled photographic span were not supported in fresh reserve tests. The current evidence therefore points to a robust geometric regularity in how flower-colour polymorphism is expressed, while leaving its evolutionary and ecological mechanisms open.

---

## Frozen claim boundaries for submission

The manuscript may claim:

1. observer-disjoint reproducibility of continuous four-state D under the first-frozen high-depth validation rule;
2. targeted discovery/audit localization of legacy H2 geometry to white versus nonwhite;
3. untouched prospective species-disjoint confirmation of the already frozen white axis in the third cohort;
4. non-support for the tested broad phylogenetic-signal claim;
5. non-replication of the sampled-photographic-span association.

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
- Paper architecture: `docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260913.md`
- H1 reconciliation: `docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md`
- H2 target freeze: `docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`
- Third-cohort protocol: `docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md`
- Third-cohort result/claim freeze: `docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`
- Third-cohort measurement result: `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`
- Third-cohort H2 result: `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`
- H3a protocol/result manifest: `docs/POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md`, `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`
- H3b result freeze: `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`
- P500 terminal status: `docs/P500_PROSPECTIVE_TERMINAL_POSTMORTEM_20260916.md`

## References

Literature citations are intentionally not populated in this repository-grounded draft. They should be integrated only after a separate source audit so that background claims remain distinguishable from frozen empirical results.
