# Flower-colour polymorphism paper architecture

This document defines the **current manuscript architecture**. It is a reader-facing map of the active New Phytologist paper, not a development history.

Authoritative surfaces:

- submission manuscript: `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- claim ledger: `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- data-lineage map: `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`
- figure plan: `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`
- reproducibility contract: `CURRENT_PAPER_REPRODUCIBILITY.md`

## 1. Paper in one sentence

Across high-depth community-science image cohorts, within-species flower-colour variation behaves as a reproducible comparative trait, shows recurrent achromatic–chromatic geometry in a prospective species-disjoint test, and is more strongly geographically organized in species with greater colour diversity.

## 2. Questions

The paper asks three linked questions:

1. **Is within-species flower-colour variation reproducible as a species-level comparative trait?**
2. **Does within-species colour displacement repeatedly occupy a common direction in phenotype space?**
3. **Is greater colour diversity associated with stronger within-species geographic organization?**

The first establishes the inferential unit. The second tests cross-species regularity in phenotype space. The third asks whether the amount of variation is spatially structured rather than merely broad or noisy.

## 3. Data architecture and why each stage exists

| Stage | Data source | Scale | Why required | Role |
|---|---|---:|---|---|
| Outcome-blind opportunity frame | iNaturalist metadata | 42,111 species; 4,730 high-depth-capable | Define candidate species before flower colour is examined | Sampling frame only |
| Discovery | iNaturalist high-depth images | 500 species × 100 photographs; 369 D-eligible | Estimate within-species distributions and identify candidate structure | Discovery/calibration |
| Species-disjoint validation | Same acquisition contract, different species | 500 species × 100 photographs; 363 D-eligible | Test reproducibility and spatial associations outside discovery species | Validation |
| Fresh-image D transport | New photo IDs, same frozen measurement system | 136 overlapping evaluable species | Test transport beyond the original image sample | Same-system transport |
| Prospective confirmation | Newly sampled species and photo IDs after target freeze | 499 species × 100 photographs; 377 measurement-evaluable; 158 primary H2 species | Test a frozen target without reusing discovery species/photos | Primary confirmatory H2 |
| Post-confirmatory annotations | Highlight metrics, WorldClim 2.1, V.PhyloMaker2, sampled span | Existing image cohorts with targeted annotations | Bound image-formation effects and simple alternative explanations | Interpretation filters |

The central inferential sequence is:

`outcome-blind frame → discovery/validation → target freeze → prospective confirmation → bounded explanatory tests`

The prospective confirmation cohort remains within the same iNaturalist source and measurement system, so it is species- and photo-disjoint confirmation rather than independent-source replication.

## 4. Measurement layer

Flower colour is represented in two complementary ways.

### Four-state diversity

The four frozen biological states are:

- white;
- yellow/orange;
- red/pink;
- blue/purple.

For species-level state frequencies (p_k),

[
D = 1 - sum_k p_k^2.
]

D measures the amount of within-species colour diversity.

### Continuous colour geometry

Nine-colour compositions retain continuous phenotype-space information. Unlabelled within-species two-mode displacement defines a sign-invariant direction (u_i). The confirmed target is the fixed white-versus-equal-nonwhite direction (q_{white}), summarized by

[
W = operatorname{mean}_i (u_i^T q_{white})^2.
]

The confirmatory contrast is observed W relative to a construction-preserving structured null, not relative to isotropy alone.

## 5. Results spine

### Result 1 — D is reproducible and finite-sample robust

Observer-disjoint validation supports repeatable between-species D rankings, and fresh-image transport is strong across 136 overlapping species. Finite-sample correction leaves species rankings essentially unchanged (raw-versus-corrected Spearman rho > 0.99997 in both high-depth cohorts) and preserves the D–spatial, sampled-span and phylogenetic sensitivity conclusions.

Interpretation: repeated photographs can support a comparative within-species distributional trait, while measurement remains imperfect and same-system. The small plug-in bias of raw Gini–Simpson D is not driving the paper's D-based results.

### Result 2 — recurrent geometry is achromatic–chromatic

The original cohorts localize recurrent continuous displacement to a white-versus-nonwhite direction. After that target is frozen, a new species- and photo-disjoint prospective cohort confirms excess alignment relative to the coarse-state-preserving structured null.

Primary prospective result:

- n = 158 species;
- W = 0.517;
- structured-null median = 0.457;
- p = 0.001.

Interpretation: the confirmed quantity is **excess achromatic–chromatic alignment beyond the measured coarse-state construction baseline**.

### Result 3 — greater D accompanies stronger geographic organization

The association is positive in discovery and species-disjoint validation and remains supported after sampled-span, technical-failure, background and ambiguity checks.

Interpretation: greater measured diversity is spatially organized rather than explained simply by broader sampled photographic extent.

### Secondary empirical mechanism evidence and alternative-explanation boundaries

Two secondary results provide bounded biological clues about why an achromatic endpoint may recur:

- a descriptive reanalysis of Del Valle et al. (2019) shows markedly higher natural frequencies for petal anthocyanin-loss than whole-plant anthocyanin-loss phenotypes, consistent with differential persistence when pigment loss is flower-restricted;
- white records occupy warmer BIO5 environments in the frozen primary prospective analysis, but the signal weakens under observer controls and does not transport as a common cross-cohort rule; BIO5 is therefore an observer-sensitive contextual clue, not robust evidence for temperature sorting.

Additional boundaries are:
- no detectable broad tree-wide phylogenetic conservation under the tested validation design;
- the discovery sampled-span association collapses in validation;
- direct highlight analysis shows that coarse white classification is exposure-coupled.

These secondary results constrain mechanism without establishing a universal cause of the prospective H2 geometry.

## 6. Discussion logic

The Discussion follows one hierarchy.

1. **Validated distributional trait:** repeated photographs recover species differences in within-species colour diversity.
2. **What varies:** cross-species regularity is strongest along an achromatic–chromatic phenotype-space direction.
3. **Why it may recur:** multiple genetic/developmental routes can converge on achromatic phenotypes; PAL/WAL frequency contrasts support differential persistence as one filter, while the observer-sensitive, nontransporting BIO5 result remains only a contextual environmental clue.
4. **Where it is sorted:** greater D is associated with stronger within-species geographic organization, but the present data do not identify one common map or causal maintenance mechanism.
5. **What simple explanations do not suffice:** broad phylogenetic conservation, sampled photographic span and one universal BIO5 effect are unsupported under the tested designs.
6. **Scope:** source dependence and exposure coupling bound the biological interpretation.

The ecological synthesis is:

**what varies shows recurrent cross-species structure; where that variation is sorted remains context dependent and unresolved as a universal map.**

## 7. Figure sequence

1. **Figure 1:** data provenance, sampling scale and why each inferential stage exists.
2. **Figure 2:** observer-disjoint reproducibility of D.
3. **Figure 3:** localization of the white-versus-nonwhite target in the original cohorts.
4. **Figure 4:** prospective confirmation of the frozen achromatic–chromatic axis.
5. **Figure 5:** positive D–spatial association plus bounded phylogenetic and sampled-span alternatives.

This sequence mirrors the manuscript logic:

`define → validate → discover → confirm → interpret`

## 8. Current title

**Within-species flower-colour variation shows recurrent achromatic–chromatic geometry across plant species**

The title emphasizes the positive biological result without implying a universal mechanism, a shared geographic map, or independent-source replication. The structured-null boundary is stated explicitly in the Summary, Results and Discussion.

## 9. Hard nonclaims

The manuscript does not claim:

- global prevalence of flower-colour polymorphism;
- independent-source replication;
- artifact-free measurement of white;
- pigment chemistry or evolutionary transition direction;
- universal pollinator or climate causation;
- a universal BIO5–white relationship;
- a shared or species-specific geographic map across all species;
- absence of finer-scale phylogenetic effects;
- that sampled photographic span is true biological range size;
- that the D–spatial association identifies a causal maintenance mechanism.

## 10. Submission surface

The active submission surface consists of:

- the New Phytologist manuscript;
- five main figures;
- one data-provenance/inferential-necessity table;
- Supporting Information;
- cover letter;
- claim ledger;
- data-lineage map;
- frozen machine-readable results and reproducibility package.

Development-only execution history is preserved for provenance but is not part of the manuscript's reader-facing architecture.
