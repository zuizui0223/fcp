# Flower-colour polymorphism paper architecture

This document defines the **current manuscript architecture**. It is a reader-facing map of the active New Phytologist paper, not a development history.

Authoritative surfaces:

- submission manuscript: `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`
- claim ledger: `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- data-lineage map: `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`
- figure plan: `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`
- reproducibility contract: `CURRENT_PAPER_REPRODUCIBILITY.md`

## 1. Paper in one sentence

Species-wide flower-colour ITV is repeatedly **geographically partitioned** across plant species: nearby conspecific observations are more colour-homogeneous than expected from each species' own colour composition, while geographic turnover contains a dominant IBD-like component and a smaller technically bounded environmental component; a separate prospective test shows that phenotype-space displacement also recurrently contains an achromatic–chromatic component.

## 2. Questions

The paper asks four linked questions:

1. **Can species-wide flower-colour distributions be measured reproducibly?**
2. **Does range-wide colour diversity reside as local coexistence, or is it distributed among geographic localities?**
3. **Is geographic colour turnover explained only by distance/history, or do environmental differences add nonredundant IBE-like structure?**
4. **Does within-species colour displacement repeatedly occupy a common direction in phenotype space?**

Together these questions distinguish **overall sampled diversity**, **spatial allocation**, **spatial process** and **phenotype-space direction**. D is not a within-population polymorphism measure: it can combine local coexistence with differentiation among sampled locations.

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

D measures species-wide sampled colour-state diversity across retained observations; it can reflect local coexistence, differentiation among sampled locations, or both.

### Continuous colour geometry

Nine-colour compositions retain continuous phenotype-space information. Unlabelled within-species two-mode displacement defines a sign-invariant direction (u_i). The confirmed target is the fixed white-versus-equal-nonwhite direction (q_{white}), summarized by

[
W = operatorname{mean}_i (u_i^T q_{white})^2.
]

The confirmatory contrast is observed W relative to a construction-preserving structured null, not relative to isotropy alone.

## 5. Results spine

### Result 1 — species-wide D is measurable and reproducible

Observer-disjoint validation and fresh-image transport show that repeated photographs recover stable between-species differences in species-wide sampled four-state diversity. This establishes the comparative distribution before asking where it is spatially allocated.

### Result 2 — species-wide flower-colour ITV is geographically partitioned

At the fixed post hoc 50-km scale, local pairwise diversity is lower than expected from each species' exact overall colour composition in discovery, validation and the third cohort:

- discovery: depletion = **0.02053**, p = **0.005**;
- validation: depletion = **0.01867**, p = **0.005**;
- third: depletion = **0.01468**, p = **0.005**.

The result is also supported at 25, 100 and 250 km in all three cohorts.

Direct falsification tests remain positive across all three cohorts after:
- excluding same-observer pairs;
- removing all white records;
- replacing coarse states with continuous nine-colour distances.

Interpretation: range-wide FCP is often a **distributed polymorphism**, not unrestricted local mixing of all species-wide variants.

### Result 3 — geographic turnover contains both IBD-like and bounded IBE-like components

Continuous nine-colour turnover shows positive unique geographic and BIO5-associated components in discovery, validation and the third cohort. The IBD-like component is larger in every cohort.

The BIO5 component is not promoted to adaptive causation because:
- validation flower-minus-background response is unsupported;
- strict same-observer-pair sensitivities are unsupported.

Interpretation: geography is the dominant structured component; a smaller environmental residual is a candidate, technically bounded sorting signal.

### Result 4 — recurrent phenotype-space geometry is achromatic–chromatic

The original cohorts localize recurrent continuous displacement to a white-versus-nonwhite direction. After that target is frozen, a new species- and photo-disjoint prospective cohort confirms excess alignment relative to the coarse-state-preserving structured null:

- n = 158;
- W = **0.517**;
- null median = **0.457**;
- p = **0.001**.

Interpretation: the confirmed quantity is excess achromatic–chromatic alignment beyond the measured coarse-state construction baseline.

### Result 5 — simple universal explanations fail

- Broad environmental heterogeneity does not transport as a general explanation for species-wide D.
- A universal signed warm-to-white BIO5 rule is unsupported across cohorts.
- Broad tree-wide phylogenetic conservation is not detected under the tested validation design.
- The discovery sampled-span association collapses in validation.
- Coarse white classification remains exposure-coupled.

Mechanistic molecular context remains in Discussion/SI and does not upgrade the ecological results.

## 6. Discussion logic

The Discussion follows this hierarchy.

1. **Validated distributional trait:** species-wide colour distributions can be measured reproducibly.
2. **Where ITV resides:** range-wide flower-colour diversity is locally depleted and therefore geographically distributed.
3. **Evolutionary process boundary:** IBD-like structure is stronger than the bounded BIO5-associated residual; selection, gene flow, drift, history and plasticity remain competing mechanisms.
4. **What varies:** a separate prospective test confirms recurrent achromatic–chromatic displacement geometry.
5. **What current data cannot establish:** local adaptation requires fitness/genetic evidence; reproductive isolation requires mating/gene-flow evidence.
6. **Mechanistic context:** published pigment-network work supplies developmental plausibility only.

The ecological synthesis is:

**species-wide flower-colour ITV is not merely individual noise or ubiquitous local coexistence; a reproducible component is distributed among geographic localities.**

The evolutionary interpretation is deliberately bounded:

**spatially varying selection is one biologically plausible generator of distributed FCP, but the present comparative photographs identify the pattern rather than proving adaptation.**

## 7. Figure sequence

1. **Figure 1:** data provenance, sampling scale and inferential chronology.
2. **Figure 2:** observer-disjoint reproducibility of D.
3. **Figure 3:** localization of the achromatic–chromatic target in the original cohorts.
4. **Figure 4:** prospective confirmation of the frozen achromatic–chromatic axis.
5. **Figure 5:** distributed polymorphism and spatial process — local depletion across cohorts, robustness diagnostics, and IBD-like versus IBE-like turnover.

Phylogenetic and sampled-span negative tests move to Supporting Information.

This sequence mirrors the manuscript logic:

`measure → locate ITV geographically → decompose spatial process → retain complementary phenotype geometry`

## 8. Current title

**Species-wide flower-colour variation is geographically partitioned across plant species**

The title centers the most replicated ecological pattern: local colour diversity is depleted relative to each species' own range-wide composition in three species-disjoint cohorts and survives observer, white-state and continuous-colour falsification tests. It does not imply that the partitioning is adaptive.

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
