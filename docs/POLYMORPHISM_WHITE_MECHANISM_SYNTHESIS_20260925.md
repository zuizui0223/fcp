# Recurrent white/nonwhite flower-colour geometry: mechanism synthesis — 2026-09-25

## Executive statement

The current evidence does **not** support one universal ecological cause of recurrent white/nonwhite flower-colour variation.

It supports a different, more mechanistic model:

> **Many genetic/regulatory changes and some thermal environments can perturb overlapping parts of the floral pigment network and converge on a pale/achromatic endpoint; the ecological conditions that expose, maintain or sort that endpoint remain species- and context-dependent.**

This synthesis is separate from the frozen New Phytologist manuscript.

---

## 1. Phenotype-space generality — supported

The prospective third cohort confirmed excess within-species alignment with the frozen white-versus-nonwhite axis relative to the coarse-state-preserving structured null.

This establishes recurrent achromatic–chromatic geometry.

It does **not** establish:
- transition direction;
- pigment chemistry;
- environmental causation;
- pollinator causation.

---

## 2. Recurrent nonwhite -> white evolutionary direction — not supported

OpenTree endpoint analysis:
- 196 species;
- 116 nonwhite-dominant;
- 80 white-dominant;
- 100/100 topology resolutions completed;
- median q(nonwhite->white) / q(white->nonwhite) = **0.7007**;
- fraction with ratio > 1 = **0**.

Frozen verdict:

`DIRECTIONAL_ASYMMETRY_NOT_SUPPORTED_UNDER_THIS_TEST`.

Therefore “recurrent white/nonwhite axis” must not be rewritten as “repeated evolution from pigmented to white”.

---

## 3. Universal environmental cause — not supported

### BIO5

Third cohort alone showed a prospective positive within-species BIO5–white association.

But the fixed replication in the original species-disjoint cohorts failed:

**Discovery**
- 271 eligible species;
- median white-minus-nonwhite BIO5 contrast = -0.0068 SD;
- Wilcoxon p = 0.743;
- conditional OR = 1.033, p = 0.131.

**Reserve**
- 260 eligible species;
- median contrast = +0.0662 SD;
- Wilcoxon p = 0.0541;
- conditional OR = 1.048, p = 0.0319.

Frozen replication verdict:

`LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST`.

BIO5 is therefore a cohort-dependent clue, not a common/global FCP driver.

### Seasonal heat

The prospective seasonal-heat follow-up showed strong broad/species-stratified signals, but the same-species × same-climate-cell criterion was p = 0.054995 and the full gate failed.

Thus heat is plausible in particular systems, but not confirmed as a universal cross-species sorter.

---

## 4. Universal pollinator cause — not supported

Independent tests did not support:
- Sphingidae as a general whitening predictor;
- Lepidoptera fraction;
- effort-adjusted bee partner breadth;
- local bee richness;
- local Bombus fraction.

Pollinators can be decisive in individual systems, but no common global pollinator explanation is supported by the tested datasets.

---

## 5. Natural proximal molecular convergence — supported as a mechanistic synthesis

At least five independent natural white/pigmented systems converge on reduced, blocked or rerouted anthocyanin/flavonoid-network output:

| Species | Natural white/pigmented mechanism |
|---|---|
| *Gymnadenia rhellicani* | R2R3-MYB stop variant -> ANS regulation -> cyanidin reduction |
| *Ipomoea purpurea* | CHS loss-of-function and bHLH regulatory variants |
| *Parrya nudicaulis* | petal-specific CHS downregulation / cis-regulatory control |
| *Pleione limprichtii* | FLS, ANS, UFGT and candidate MBW-network regulation |
| *Silene littorea* | strong F3H downregulation with candidate Myb1a regulation |

This is **pathway-level convergence**, not one universal white-flower gene.

Important bound:
- *Abronia fragrans* is a betalain-system white/pink polymorphism, so anthocyanin loss cannot be universal across angiosperms.

---

## 6. Heat acts on the same molecular network — supported mechanistically

Normalized FCP natural-white nodes:

`MYB, CHS, F3H, ANS, UFGT, FLS, bHLH, WD40`

Direct floral heat-response literature includes:

`MYB, CHS, CHI, F3H, DFR, ANS, UFGT`

Explicit overlap:

`MYB, CHS, F3H, ANS, UFGT`

= **5/8 normalized FCP node classes**.

This is descriptive mechanistic overlap, not an enrichment test.

---

## 7. Moricandia gives a direct within-species thermal bridge

### 7a. Temperature dose under the same summer photoperiod

Published Source Data contain 15 individuals followed through:

`14.17 C spring -> 23.75 C mild summer -> 28.75 C hotter summer`.

The mild and hotter summer regimes share the stated 16/8 h summer photoperiod.

After averaging two flowers within individual-period:

- mean cyanidin: **1.893 -> 0.531**;
- hot/mild ratio = **0.281**;
- mean reduction = **71.9%**;
- **13/15** individuals decreased;
- paired Wilcoxon one-sided p = **0.000427**;
- sign-test p = **0.00369**.

### 7b. Reversibility

Separate 14-individual sequence:

`23.75 C mild summer -> 14.17 C spring return`

- cyanidin: **1.236 -> 4.708**;
- **3.81x increase**;
- 12/14 individuals increased;
- one-sided Wilcoxon p = **0.000305**.

Spring-like -> hot sequence:
- 14/14 individuals decrease;
- cyanidin **5.690 -> 0.179**;
- p = **6.10e-5**.

The period-2-to-period-3 change differs between mild->hot and mild->spring-return sequences:
- Mann-Whitney p = **3.72e-5**.

Thus a generic “third measurement period causes fading” explanation is inconsistent with the reversal data.

### 7c. Same system, same pigment network

Public transcriptome reanalysis of the summer-white state shows significant downregulation of every normalized FCP node represented in the published significant-DEG table:

- MYB90: logFC **-4.056**, FDR 4.81e-11;
- CHS: **-1.293**, FDR 0.034;
- U75C1/U78D2: four significant rows, median logFC **-1.5315**;
- TTG1/WD40: **-1.411**, FDR 0.0279.

Additional pathway genes PAL, 4CL and DFR are also lower.

Therefore Moricandia provides the current strongest experimental chain:

`
greater thermal severity
   -> lower floral anthocyanin
   -> lower expression at FCP-overlapping pigment-network nodes
   -> white / achromatic floral phenotype
`

The RNA-seq spring-versus-mild-summer comparison changes both temperature and photoperiod; only the mild-versus-hot anthocyanin contrast holds the stated summer photoperiod constant.

---

## 8. Natural same-species ecological bridge: Parrya

In *Parrya nudicaulis*:
- natural white-flower frequency increases with growing-season temperature across Alaskan populations;
- the natural white phenotype maps to reduced flux near the entry of the anthocyanin pathway, with CHS-associated regulatory downregulation.

This is strong ecological × molecular compatibility, but the temperature-frequency association remains observational.

---

## 9. Context dependence is expected, not an exception

*Ipomoea purpurea* provides an important counterexample: the anthocyanin network is temperature-sensitive, but the sign is not a universal whitening effect.

Published heat responses in other species also depend on:
- genotype;
- flower developmental stage;
- light environment;
- regulatory background.

This explains why:
- heat can be a real molecular route to an achromatic endpoint;
- yet BIO5 can fail as a universal cross-species main effect.

---

## 10. Current working model

`
         genetic structural lesions
        / 
regulatory variants -----------\
                                \
                                 -> pigment-network flux -> chromatic endpoint
                                /
thermal / environmental input -/
`

More explicitly:

`
(genetic/regulatory state) × (thermal environment)
        -> anthocyanin/flavonoid network output
        -> white/pale versus pigmented flower
`

### What is general

The **network architecture and accessible phenotypic endpoint** may be recurrent.

### What is contingent

The sign and magnitude of ecological selection or plastic response are species-, genotype- and environment-specific.

---

## 11. Decisive next experiment

The strongest next experiment is not another global BIO5 regression.

It is:

`natural colour genotype/morph × temperature`

with temperature manipulated while photoperiod, light and water are controlled, measuring on the same flowers:

1. spectral colour / CIELAB;
2. anthocyanin concentration;
3. flavonol/flavone allocation;
4. CHS, F3H, DFR, ANS, UFGT;
5. MYB, bHLH and WD40/TTG1 regulators.

The primary estimand should be the **morph/genotype × temperature interaction**.

The key prediction is not simply “hotter = whiter”. It is:

> genotypes already carrying reduced pigment-network activity should differ predictably in the magnitude or threshold of thermal pigment suppression.

---

## Claim ceiling

Supported working mechanism:

> Distinct genetic and thermal perturbations can converge on overlapping pigment-network nodes and thereby generate similar achromatic floral phenotypes.

Not supported:
- universal high-temperature whitening;
- universal pollinator-driven whitening;
- repeated nonwhite->white macroevolutionary asymmetry;
- one universal anthocyanin gene;
- heat as the explanation of the global FCP H2 result.
