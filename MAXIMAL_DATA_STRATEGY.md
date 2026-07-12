# Maximal data strategy

The project no longer treats a fixed sample size as the ceiling. The empirical universe is now defined as a reproducible all-angiosperm census from a taxonomic backbone, with evidence layers added independently.

## Core principle

Do not ask only:

> Which species are known to be flower-colour polymorphic?

That creates severe ascertainment bias. Instead construct:

1. a maximal taxonomic universe of accepted angiosperm species;
2. an evidence layer for pollination mode and floral display opportunity;
3. an outcome layer for flower-colour structure;
4. an effort layer describing how strongly each species has been searched.

Every species starts as `outcome_class = unknown`. Unknown is never silently converted to monomorphic.

## Outcome states

- `monomorphic_confirmed`
- `geographic_mosaic`
- `mixed`
- `local_coexistence`
- `unknown`
- `excluded_trait_definition`

The distinction between metapopulation persistence and local coexistence is retained throughout.

## Maximal census workflow

```bash
python build_angiosperm_census.py --out data/angiosperm_census.csv
python prioritize_census_review.py \
  --input data/angiosperm_census.csv \
  --out data/angiosperm_census_prioritized.csv
```

Use `--max-records 1000` only for smoke testing. Omit it for the full run.

## Evidence accumulation without bias

The full census is immutable except for taxonomic-version updates. Evidence can be accumulated asynchronously in separate columns or tables.

Recommended evidence tables:

- `species_pollination_evidence.csv`
- `species_display_evidence.csv`
- `species_colour_outcome_evidence.csv`
- `species_sampling_effort.csv`
- `species_range_environment.csv`

Each record should include source, evidence type, date, extraction method, confidence, and reviewer.

## Review priority is not inclusion probability

The prioritization script creates four review queues:

- A: families rich in documented floral-colour systems;
- B: wind-pollinated controls;
- C: specialized pollination systems;
- D: general background.

All taxa remain in the census. Priority controls review order only and must not be used as a response variable.

## Analysis layers

### Layer 1: occurrence of any documented flower-colour polymorphism

Model detection explicitly. `unknown` is not a negative.

### Layer 2: spatial organization conditional on polymorphism

Compare:

- geographic mosaic;
- mixed spatial structure;
- local coexistence.

### Layer 3: display-opportunity sensitivity analysis

Repeat analyses for:

1. all angiosperms;
2. species with a visible floral display;
3. species with comparable petal/tepal displays.

### Layer 4: pollination system

Retain abiotic and biotic systems rather than excluding wind-pollinated taxa. At minimum distinguish:

- wind;
- water;
- insect;
- bird;
- bat;
- other vertebrate;
- mixed/ambophilous;
- unknown.

## Why maximal coverage matters

Range size, research effort, horticultural interest, conspicuousness, and pollination mode all affect whether a polymorphism is reported. A maximal census allows these observation processes to be modeled instead of conditioning the dataset on prior discovery.

## Practical stopping rule

There is no fixed target `N`. Data collection continues until one of the following is reached:

1. the complete census has been taxonomically frozen and all high-priority evidence queues are exhausted;
2. additional review produces negligible change in posterior estimates;
3. a preregistered coverage threshold is reached for every major phylogenetic and pollination stratum.

The full census is the sampling frame; manual review depth is the resource-limited quantity.
