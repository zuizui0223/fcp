# Moricandia mechanistic chain linking heat to an achromatic floral endpoint — 2026-09-25

## Evidence chain

This document integrates only evidence that is directly supported by the public Moricandia experiment or its public source files. It is post-publication mechanistic triangulation, not a new confirmatory study.

### 1. Environmental manipulation

The published experiment separates two summer temperature regimes that share the stated 16/8 h photoperiod:

- mild summer: 30/20 C day/night, weighted mean 23.75 C;
- hot summer: 35/25 C day/night, weighted mean 28.75 C.

This contrast therefore isolates an increase in thermal severity while holding the stated summer photoperiod constant.

### 2. Floral anthocyanin responds dose-dependently to hotter summer conditions

In the 15 individuals followed from mild summer to hot summer:

- mean cyanidin-equivalent concentration: **1.893 -> 0.531**;
- hot/mild mean ratio: **0.281**;
- mean reduction: **71.9%**;
- **13/15** individuals decreased;
- paired Wilcoxon one-sided p = **0.000427**;
- two-sided p = **0.000854**.

Reversibility strengthens the interpretation:

- mild summer -> spring return: 1.236 -> 4.708, 12/14 individuals increased, one-sided p = 0.000305;
- spring-control -> hot summer: 5.690 -> 0.179, 14/14 individuals decreased, one-sided p = 0.000061.

Thus the anthocyanin response tracks thermal treatment direction rather than a generic monotonic time decline.

### 3. The flavonol comparator is largely maintained

In the same 15 mild -> hot individuals:

- mean kaempferol-equivalent flavonol: **58.44 -> 57.18**;
- hot/mild mean ratio: **0.979**;
- 8/15 individuals decreased and 7/15 increased;
- paired two-sided Wilcoxon p = **0.679**.

The preregistered log-ratio branch-specificity statistic is not estimable because two individuals have measured zero cyanidin means. No pseudocount is introduced and the frozen primary verdict remains `ANTHOCYANIN_BRANCH_SPECIFICITY_NOT_SUPPORTED_UNDER_THIS_TEST`.

A zero-safe post-open descriptive robustness using `(hot-mild)/(hot+mild)` gives:

- median symmetric cyanidin change = **-0.673**;
- median symmetric kaempferol change = **-0.030**;
- cyanidin-minus-kaempferol response <0 in **12/15** individuals;
- paired one-sided Wilcoxon p = **0.00513**;
- sign-test p = **0.0176**.

This robustness cannot rescue the frozen primary, but it is consistent with stronger suppression of the anthocyanin branch than of the kaempferol/flavonol branch.

### 4. The pigment-network nodes overlap the natural-FCP white-state network

Public spring-versus-mild-summer RNA-seq shows significant lower expression in summer-white flowers for every frozen FCP white-state node class represented in the significant-DEG table:

- **MYB90**: logFC -4.056, FDR 4.81e-11;
- **CHS**: logFC -1.293, FDR 0.034;
- **U75C1/U78D2**: four significant transcripts, median logFC -1.5315;
- **TTG1/WD40**: logFC -1.411, FDR 0.0279.

Additional pathway context genes are lower:
- PAL;
- 4CL;
- DFR (logFC -3.164).

The transcript comparison changes both temperature and photoperiod, so it cannot by itself identify temperature as the sole transcriptional cue. The independent mild-versus-hot anthocyanin dose test supplies the temperature-specific physiological bridge.

### 5. The endpoint is a white, UV-absorbing flower

The peer-reviewed study reports that the same individuals produce bright lilac, UV-reflecting flowers in spring and small white, UV-absorbing flowers in summer. Anthocyanin derivatives generate the spring lilac colour, while flavonols remain important UV absorbers.

The public Source Data workbook does not include the Figure 2b raw reflectance spectra, so this final phenotype contrast cannot be independently recomputed numerically from Source Data and is retained at the published evidence level.

## Integrated mechanistic model

The directly supported within-species chain is:

`
higher summer thermal severity
        ->
strong reduction in floral anthocyanin
        +
relative maintenance of kaempferol/flavonol
        ->
overlapping downregulation of pigment-network structural/regulatory nodes
        ->
white / UV-absorbing summer floral phenotype
`

The ordering of transcript and pigment evidence should not be overinterpreted as a formal mediation analysis. What is supported is that thermal treatment, pigment abundance, pigment-network expression and floral phenotype move coherently within the same biological system.

## Connection to FCP

This system gives a concrete mechanistic realization of the FCP many-to-one accessibility model:

- natural genetic/regulatory variants in multiple FCP species perturb CHS, F3H, ANS/UFGT and MYB/bHLH/WD40-related nodes;
- environmental heat can perturb overlapping parts of the same network;
- distinct perturbations can therefore converge on a low-chromatic / achromatic endpoint.

The general prediction is **genotype/morph x temperature interaction**, not a universal positive temperature main effect.

## Hard boundaries

This chain does not:
- rescue the failed discovery/reserve BIO5 replication;
- establish that heat is the universal ecological cause of the FCP white/nonwhite axis;
- establish repeated nonwhite -> white evolutionary direction;
- prove isotope-traced metabolic flux rerouting;
- generalize anthocyanin mechanisms to betalain or carotenoid white systems.
