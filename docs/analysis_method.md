# Research analysis protocol

## Central question
Which ecological and geographic conditions predict natural flower-colour polymorphism after controlling for literature and observation effort?

## Two-stage design

### Stage 1: ascertainment model
The present repository contains literature-derived candidates rather than a random sample of angiosperms. We first model whether a candidate is validated as a natural polymorphism as a function of discovery effort, using log1p(n_works). Family-level standardized residuals identify clades with more or fewer validated cases than expected at equal effort. Leave-one-family-out refits test whether the effort effect is driven by a single large family.

This stage estimates validation conditional on candidature. It does not estimate biological prevalence.

### Stage 2: biological macroecology model
When data/species_macroecology_covariates.csv is present, the pipeline fits a binomial model with island status, absolute latitude, life history, pollination system and publication effort as predictors, using family-clustered robust standard errors.

The covariate contract also reserves growth form, breeding system, range size, GBIF occurrences, accepted names, taxon keys and phylogeny tip labels for later models.

## Required extension for prevalence inference
A defensible global prevalence model requires a sampled background of species that were eligible for detection, not only positive or ambiguous candidates. The target background should be stratified by family and region. Candidate or validated status can then be modelled jointly with publication effort and occurrence availability.

## Final model hierarchy
1. Detection model for entry into the literature candidate set.
2. State model for natural flower-colour polymorphism as a function of ecology and geography.
3. Strict, moderate and inclusive evidence sensitivity analyses.
4. Family- and phylogeny-aware random effects.
5. Interaction models only after main effects are stable.

## Interpretation rules
- Never call the candidate fraction a global prevalence estimate.
- Report effort-adjusted effects alongside raw patterns.
- Treat family residuals as hypothesis-generating until phylogenetic modelling is available.
- Keep acquisition code in PR1 and biological analysis code in PR2.
