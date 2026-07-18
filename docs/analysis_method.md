# Evidence analysis method

## Scope

The first analysis stage describes the evidence dataset as an ascertainment-biased literature sample. It does **not** estimate the global biological prevalence of flower-colour polymorphism.

## Analysis units

- Primary unit: canonical species.
- Secondary aggregation: plant family.
- Evidence variables: number of retained works, title matches, context matches, maximum evidence score and total evidence score.

## Operational evidence tiers

- **A — strong:** at least two title matches, two context matches, two works and maximum score at least 24.
- **B — moderate:** at least one title match, one context match and maximum score at least 18.
- **C — weak:** at least one title match or two context matches.
- **D — minimal:** all remaining candidates.

These are review-prioritisation tiers, not confirmed biological states.

## Required diagnostics

1. Evidence-tier counts.
2. Family representation and Wilson intervals for the fraction of candidates with strong or moderate evidence.
3. Family concentration to expose taxonomic sampling imbalance.
4. Threshold sensitivity across title, context and score cutoffs.
5. Explicit reporting of literature ascertainment and taxonomic coverage limitations.

## Later inferential stage

Biological comparisons should begin only after manual validation produces a defensible species-level response variable. The preferred sequence is:

1. Freeze inclusion and exclusion criteria.
2. Blind or duplicate-review a validation subset and calculate agreement.
3. Resolve synonyms and accepted names.
4. Add sampling-effort covariates, including publication count and database coverage.
5. Fit family- or phylogeny-aware models, with sensitivity analyses across evidence thresholds.
6. Separate natural polymorphism, ontogenetic colour change and artificial or horticultural variation.

The baseline script establishes the reproducible descriptive layer required before those models are attempted.
