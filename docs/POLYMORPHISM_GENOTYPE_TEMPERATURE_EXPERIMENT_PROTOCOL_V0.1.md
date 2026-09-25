# FCP colour genotype/morph × temperature experiment protocol — v0.1

## Purpose

Directly test the mechanism suggested by the current FCP evidence:

`
pigment-network genotype / regulatory state × temperature
    -> pigment-network output
    -> chromatic versus achromatic floral phenotype
`

This protocol is prospective design work. It does not alter the frozen New Phytologist manuscript and it is not a preregistration until a focal species/population panel and execution date are frozen separately.

## Primary biological hypothesis

Natural white and pigmented morphs differ in their thermal reaction norms because they occupy different states of the floral pigment network.

The primary prediction is **not** a universal hotter -> white main effect.

The primary prediction is:

`
morph/genotype × temperature interaction != 0
`

for at least one pigment-network phenotype.

## Experimental structure

### Biological material

Use naturally segregating white and pigmented material from the same species.

Preferred design:
- >= 12 independent maternal families / genotypes per morph;
- >= 2 offspring / ramets per genotype per temperature treatment;
- genotype identity retained throughout the analysis;
- source population recorded and balanced across morphs where possible.

If clonal replication is unavailable, maternal family is the genetic grouping unit and the interpretation is family-level rather than genotype-level.

### Temperature treatments

Minimum:
- control temperature;
- warm treatment approximately +5 C in daytime and nighttime mean relative to control.

Photoperiod, irradiance, soil moisture, nutrients, pot size, chamber humidity and developmental stage must be held as constant as practicable.

If two growth chambers are used, treatment must be crossed with chamber over time or replicated across chambers. Temperature must not be permanently confounded with chamber.

### Flower sampling

Sample flowers at the same developmental stage, preferably first day of anthesis.

For each biological replicate, the same flower or immediately adjacent flowers from the same inflorescence should supply:
1. reflectance / image colour;
2. pigment chemistry;
3. RNA if feasible.

Record flowering date because thermal treatment may shift phenology.

## Frozen primary outcomes

### P1 — chromatic phenotype

Primary colour metric:
- continuous achromatic–chromatic coordinate derived from calibrated petal reflectance or image colour.

Preferred:
- spectral reflectance 300–700 nm plus a human-independent chromatic summary;
- CIELAB chroma or an equivalent continuous chromatic coordinate.

White/nonwhite manual class is secondary.

### P2 — anthocyanin output

Primary biochemical outcome:
- total anthocyanin concentration per fresh or dry petal mass.

Record flavonol concentration separately.

### P3 — pigment-network expression

Target panel:
- CHS
- F3H
- DFR
- ANS
- UFGT / relevant anthocyanin UGT
- MYB activator/repressor candidate
- bHLH
- WD40 / TTG1-class regulator where identifiable

Expression is normalized to prevalidated reference genes.

## Primary models

For continuous outcome y:

`
y ~ morph * temperature + (1 | genotype_or_family) + block/chamber terms
`

The coefficient of interest is:

`
morph:temperature
`

Run separately for:
- chromatic phenotype;
- anthocyanin concentration.

Expression outcomes are secondary mechanistic outcomes unless the focal regulatory genes are frozen before opening treatment results.

## Key mechanistic predictions

### H1 — pigmented morph thermal fading

If the pigmented morph uses a heat-sensitive anthocyanin configuration:

- pigmented flowers should lose chroma and anthocyanin under warming;
- white flowers near the pigment floor should change less;
- interaction should therefore shrink the white–pigmented phenotypic distance under heat.

### H2 — pathway-node concordance

In the pigmented morph under heat, one or more of the FCP-convergent nodes should decrease:
- MYB
- CHS
- F3H
- ANS
- UFGT / UGT
- WD40 / TTG1

The direction should be concordant with the observed reduction in anthocyanin.

### H3 — branch specificity

If heat selectively suppresses the anthocyanin branch rather than all flavonoid metabolism:

- anthocyanin should decrease more strongly than UV-absorbing flavonols;
- the anthocyanin/flavonol allocation ratio should decline.

This ratio must be analysed only when both quantities are strictly positive or when a measurement model for zero / below-detection values is specified before outcome opening. No post hoc pseudocount is allowed.

## Controls that matter

Mandatory:
- calibrated white/grey reference for images or spectroscopy;
- fixed developmental stage;
- fixed photoperiod and irradiance;
- randomized plant position;
- chamber replication/crossover;
- soil-water monitoring;
- genotype/family tracking.

Strongly recommended:
- leaf anthocyanin/flavonoid measurement to distinguish petal-specific from whole-plant pathway effects;
- pollen viability / seed set to test pleiotropic thermal cost;
- recovery treatment returning warmed plants to control temperature to test reversibility.

## Decision logic

The mechanism receives strong support only if:

1. morph × temperature interaction is detected for anthocyanin or chromatic phenotype;
2. colour and anthocyanin change in concordant directions;
3. at least one frozen pigment-network node changes in the predicted direction;
4. the effect is not explained by chamber or developmental-stage imbalance.

A temperature main effect without a morph/genotype interaction is evidence for plasticity but does not establish genotype-specific FCP mechanism.

A morph main effect without a temperature interaction establishes constitutive colour differentiation but not thermal sensitivity.

## Interpretation ceiling

A positive experiment would support:

> naturally segregating floral pigment states differ in thermal reaction norm through the same pigment network that repeatedly generates white/pigmented variation.

It would not by itself establish:
- the historical direction of evolutionary transitions;
- that temperature maintains the polymorphism in nature;
- a universal heat response across angiosperms;
- a universal anthocyanin mechanism in betalain/carotenoid systems.

## Link to current FCP evidence

This experiment is motivated by:
- recurrent FCP white/nonwhite phenotype geometry;
- five natural systems with pigment-network molecular convergence;
- Moricandia within-individual temperature-dose suppression of cyanidin;
- Moricandia summer-white downregulation of MYB90, CHS, U75C1/U78D2 and TTG1;
- Ipomoea CHS-null genotype × heat fitness interaction;
- failure of a universal BIO5 effect across FCP cohorts.

The experiment is designed specifically to test the interaction model implied by those mixed comparative results.
