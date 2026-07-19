# Ecological hypotheses and full analysis plan for flower-colour polymorphism

## Scope

The current evidence does not support treating all flower-colour polymorphism (FCP) as one biological phenomenon. The analysis must distinguish:

1. within-population coexistence of discrete morphs;
2. among-population geographic colour differentiation;
3. continuous colour variation;
4. plastic colour shifts;
5. artificial horticultural variation.

Only the first category is strict within-population FCP. The second is a spatial differentiation problem and may arise through local adaptation, colonisation history, drift, or gene flow.

## Competing mechanisms

### H1 Pollinator-mediated balancing selection

Predictions:
- strongest in obligately outcrossing and deceptive species;
- rare morph advantage, temporal reversals, or morph-specific visitation;
- stronger within-population coexistence than among-population replacement;
- interaction with pollinator diversity, visual-system diversity, and reward status.

### H2 Pollinator-mediated divergent selection

Predictions:
- morph frequencies track local pollinator assemblages;
- reciprocal transplant or common-garden experiments show local-morph fitness advantage;
- colour divergence exceeds neutral population structure;
- more common among populations than within populations.

### H3 Abiotic stress and pigment pleiotropy

Predictions:
- anthocyanin-presence/absence polymorphisms dominate;
- pigmented morphs increase under drought, high UV, cold, or stressful soils;
- unpigmented morphs may outperform under benign conditions;
- effects depend on pigment pathway and whether floral and vegetative pigmentation are genetically coupled.

### H4 Pigment-production or allocation trade-off

Predictions:
- pigmentation covaries negatively with flower number, flower size, nectar, seed production, or growth under benign conditions;
- trade-off is strongest for costly pathway activation or nitrogen-containing betalains;
- the same morph may gain stress tolerance but lose display investment.

This hypothesis cannot be tested from flower-colour labels alone. Pigment class and allocation traits are required.

### H5 Herbivore, florivore, pathogen, and reproductive-interference selection

Predictions:
- morph fitness differences remain after controlling for pollinator visitation;
- colour correlates with damage, pathogen load, hybridisation, or heterospecific pollen transfer;
- opposing agents can maintain polymorphism even when pollinators favour one morph.

### H6 Neutral drift, founder effects, and range-expansion history

Predictions:
- morph-frequency differentiation follows neutral genomic structure;
- no consistent environment or pollinator association after controlling for ancestry;
- strongest in small, isolated, recently colonised, or range-edge populations;
- geographic colour sectors or clines coincide with colonisation routes.

### H7 Polymorphism facilitates niche breadth or range expansion

Predictions:
- within-population FCP precedes or predicts broader realised climatic niche;
- polymorphic populations show lower temporal variance in recruitment or reproduction;
- morphs occupy partially complementary microclimates;
- species-level association remains after phylogeny, sampling effort, range age, and mating system are controlled.

A cross-sectional correlation between current range size and recorded FCP cannot establish this direction. Reverse causation is equally plausible: widespread species experience more environmental heterogeneity and are more intensively studied.

## Review of the current analysis logic

### What is valid

- Separating within-population and among-population evidence is essential.
- Controlling GBIF occurrence count addresses a major observation-effort bias.
- Same-genus matched controls are preferable to unrelated controls.
- Treating unknown controls as not-confirmed rather than monomorphic is correct.
- The current results justify a directional association statement only, not a causal mechanism.

### What remains invalid or underidentified

- Species-level binary FCP status collapses frequency, spatial scale, persistence, pigment chemistry, and plasticity.
- Range size is both a biological variable and a detection variable.
- Current controls are not verified monomorphic species.
- Taxonomic matching does not replace phylogenetic covariance modelling.
- GBIF latitude range is an incomplete proxy for geographic range and a poor proxy for niche breadth.
- A within-versus-among logistic model asks which evidence type is recorded, not what mechanism maintains colour variation.
- Literature-derived labels are affected by research traditions: orchids are more likely to be studied for pollinator NFDS, while anthocyanin systems are more likely to be studied for abiotic stress.

## Required data model

Create one row per species-population-study combination, not only one row per species.

Minimum fields:
- accepted species name and taxonomic backbone ID;
- population coordinates and year;
- within/among/continuous/plastic classification;
- morph identities and frequencies;
- pigment class: anthocyanin, carotenoid, betalain, structural, unknown;
- floral-only versus whole-plant pigmentation;
- reward status and pollination system;
- mating system and self-incompatibility;
- pollinator assemblage and visitation by morph;
- male and female fitness by morph;
- herbivory, florivory, pathogen and heterospecific-pollen measures;
- climate, soil, UV, elevation and habitat variables;
- neutral genomic structure or at least marker-based population structure;
- experimental design: observational, common garden, reciprocal transplant, manipulation;
- evidence quality and sample size.

## Full analysis sequence

### Analysis A: evidence architecture

Fit a multilevel multinomial model for evidence state:

`within / among / mixed / unclear`

Predictors:
- pigment class;
- reward/deception;
- mating system;
- life history;
- study effort;
- taxonomic random effects;
- publication and study random effects.

Purpose: identify ascertainment bias before ecological interpretation.

### Analysis B: within-population maintenance

Response:
- morph frequency through time or morph-specific relative fitness.

Models:
- hierarchical selection-gradient model;
- morph-frequency change model;
- negative frequency-dependence term;
- interactions with reward status, pollinator diversity and environmental stress.

Key test:

`relative fitness ~ morph frequency + environment + pollinator context + interactions`

A negative coefficient of morph frequency is direct evidence for rare-morph advantage.

### Analysis C: among-population differentiation

Response:
- population morph frequency or quantitative colour coordinates.

Model:
- spatial GLMM or Bayesian multinomial model;
- environmental predictors;
- pollinator predictors;
- neutral genetic covariance or ancestry matrix;
- spatial random field.

Partition explanatory power into:
- environment;
- pollinator community;
- neutral structure/history;
- shared spatial component.

### Analysis D: adaptation versus drift

Preferred tests:
- reciprocal transplant/common garden where available;
- QST-FST or animal-model equivalent for quantitative colour;
- genotype-environment association for pigment loci;
- comparison of colour clines with neutral genomic clines;
- isolation-by-environment versus isolation-by-distance.

Interpretation:
- colour differentiation exceeding neutral expectation and predicting fitness supports selection;
- colour tracking neutral structure without fitness or environment association supports drift/history;
- both may act together.

### Analysis E: pigment cost and allocation trade-offs

Restrict to studies with pigment chemistry or validated pigment-presence/absence.

Responses:
- flower number, flower size, nectar, biomass, seed set, drought survival.

Model:

`trait or fitness ~ pigment morph * stress treatment + phylogeny + study random effect`

The critical evidence is a morph-by-environment crossover, not a universal cost of colour production.

### Analysis F: niche breadth and range expansion

Use climatic niche breadth based on cleaned occurrences and environmental space, not latitude range alone.

Metrics:
- multivariate climatic hypervolume;
- within-species environmental variance;
- occupancy-modelled area;
- marginality and specialisation indices.

Comparative models:
- phylogenetic logistic/ordinal model for FCP occurrence;
- PGLS or phylogenetic Bayesian model for niche breadth;
- sister-species contrasts;
- sensitivity to study effort and range age.

Competing causal models:
1. environmental heterogeneity -> FCP;
2. FCP -> broader niche;
3. range expansion/history -> both FCP and niche breadth.

Compare them using explicit causal diagrams and structural-equation models only after temporal and phylogenetic variables are available.

## Priority order

1. Finish high-impact manual classification for unresolved species.
2. Add pigment class, reward status, mating system and evidence design to confirmed species.
3. Build a population-study evidence table for the best 30-50 species.
4. Run environment/pollinator/neutral variance partitioning on species with population data.
5. Treat the current global range-size analysis as a screening result, not the flagship causal result.
6. Make the flagship result a mechanism-aware synthesis distinguishing balancing selection, divergent selection and neutral history.

## Current defensible conclusion

Flower-colour variation is maintained by multiple interacting mechanisms. The present global data show no evidence that narrow ranges generally promote strict within-population FCP. The earlier negative range-size association is more plausibly a mixture of geographic colour differentiation, observation effort and lineage-specific research bias. Pollinator attraction, abiotic stress, allocation trade-offs and drift are competing, non-exclusive hypotheses that require population-level fitness, pigment and neutral-genetic data to distinguish.