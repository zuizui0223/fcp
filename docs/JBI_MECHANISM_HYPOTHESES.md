# JBI mechanism extension: breeding system, dispersal and phenology

Status: **prospective mechanism registry v1 (2026-08-27)**, written before inspection of the GIFT coverage output. It extends the display-core-v6 space–time analysis and is not a preregistered independent study.

## Response

The response remains the documented spatial organization of intraspecific flower-colour variation:

- C-only = local coexistence documented without retained geographic-structure evidence;
- S-only = geographic segregation documented without retained local-coexistence evidence;
- C+S = both positive axes documented;
- unresolved = not used as a biological negative.

## H6 — mating system / effective gene flow

### Basis

Self-fertilization can reduce effective pollen-mediated gene flow and increase spatial genetic structure. Classical population-genetic work explicitly predicts and observes stronger microgeographic differentiation in predominantly selfing plants, and comparative fine-scale spatial-genetic analyses report stronger spatial structure in selfers. Reviews of flower-colour polymorphism also identify autonomous selfing as a process capable of modifying FCP maintenance and geographic morph mosaics.

Key literature: Allard (1975, *Genetics* 79 Suppl:115–126); Vekemans & Hardy (2004, *Molecular Ecology* 13:921–935); Narbona et al. (2018, *Plant Biology* 20:S1, 8–20; doi:10.1111/plb.12575); Pannell (2010, Encyclopedia of Life Sciences, doi:10.1002/9780470015902.a0021909).

### Direction fixed before trait inspection

**Prediction H6:** greater self-fertilization propensity should shift documented organization toward S relative to C, conditional on geographic extent and documentation propensity.

This is not a claim that selfing creates colour divergence. The interpretation is reduced effective gene flow / increased population structure as a permissive modifier of spatial fixation or differentiation.

## H7 — dispersal × geographic fragmentation

Seed dispersal shapes population connectivity, range expansion and genetic structure, but a categorical dispersal syndrome is only a proxy for realized gene dispersal. The global FSGS literature also shows that seed-dispersal vectors matter while a simple main effect of anthropogenic fragmentation is not universal.

Key literature: Levin et al. (2003, *Annual Review of Ecology, Evolution, and Systematics* 34:575–604); Beckman et al. (2023, *Annual Review of Ecology, Evolution, and Systematics* 54); de Kort et al. / global FSGS meta-analysis (2023, PMCID: PMC10198778).

### Direction fixed before trait inspection

**Prediction H7:** dispersal syndrome is not assigned a universal monotonic main-effect direction. The primary test is an interaction:

`fragmentation × low-connectivity dispersal syndrome -> S relative to C`.

A fragmentation main effect remains a negative-control component because the v1 range-fragmentation screen was not supported.

## H9 — flowering duration × temporal environmental variability

Flowering phenology changes which pollinators and abiotic conditions a morph experiences. FCP case studies show morph-specific flowering peaks and spatiotemporal pollinator turnover, but the direction of a species-level flowering-duration main effect is not uniquely predicted: a long season may increase exposure to temporally varying selection, or average over that variation.

### Direction fixed before trait inspection

No standalone directional test is assigned to flowering duration. The primary test is:

`flowering duration × dynamic temporal environmental variability -> C relative to S`.

Calendar start/end month alone is not compared globally because hemispheres and climatic regimes make ordinal calendar month non-comparable. A cyclic flowering-duration variable is used where start and end are available.

## Trait source and missingness

Primary trait source: GIFT 3.2 aggregated species-level traits, using `bias_ref = FALSE`, `bias_deriv = FALSE`, and categorical agreement >= 0.66. Candidate columns are self-fertilization, lifecycle, dispersal syndrome, flowering start and flowering end.

No trait record = missing documentation, never biological absence or a default category.

## Coverage gates fixed before results

A candidate enters inferential modelling only if all of the following hold among the 32 informative C/S species:

1. >= 18 species have a non-missing trait value;
2. each of C-only, S-only and C+S retains >= 4 species;
3. a categorical trait has >= 2 represented categories and no primary contrast is supported by fewer than 4 species in a category;
4. an interaction is only fit when its complete-case sample also satisfies the class-count gate.

Traits failing these gates remain descriptive coverage results and are not searched for a favorable recoding.

## Inference boundary

- H6 may be tested as a main mechanism modifier if coverage passes.
- H7 is interaction-first and must not be converted post hoc into a favorable categorical main effect.
- H9 is interaction-first; flowering calendar position itself is not a global directional hypothesis.
- Family-cluster / family-bootstrap inference and geographic-extent control remain aligned with the existing v6 conditional model.
- Multiple tests within this mechanism family are adjusted; raw exploratory values may be shown but are not promoted after adjustment.
