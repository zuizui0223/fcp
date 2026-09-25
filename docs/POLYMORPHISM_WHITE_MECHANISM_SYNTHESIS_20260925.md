# FCP white-state mechanism synthesis — 2026-09-25

Status: **mechanism triage after independent environment and bee tests; transition direction remains unresolved**.

This document does not modify the frozen New Phytologist manuscript.

## Current evidence hierarchy

### 1. Recurrent phenotype axis — established, but signless

The prospective third-cohort H2 result supports excess white-versus-nonwhite achromatic–chromatic alignment relative to the frozen structured null. Because H2 uses a sign-invariant squared projection, it does not establish whether evolution proceeds preferentially from pigmented/nonwhite states toward white or in the reverse direction.

### 2. Environment — one broad-scale candidate

The prospectively frozen three-variable environmental mechanism family used 281 species after response-blind high-clip exclusion.

- BIO5 maximum temperature: **PASS** under the frozen gate.
  - median within-species white-minus-nonwhite contrast = +0.069 SD;
  - 57.3% of species positive;
  - species-level Wilcoxon p = 0.0118, Holm p = 0.0354;
  - conditional-logistic OR = 1.073 per within-species SD, 95% CI 1.029–1.119, p = 0.000919 after near-clip adjustment.
- BIO14 precipitation of the driest month: not supported.
- mean solar radiation: not supported.

The BIO5 signal is not a single-species artifact: leave-one-species-out inference retains the same conclusion and the bootstrap median contrast is positive.

### 3. Scale/observer diagnostics — BIO5 is not a local causal result

Post-result diagnostics constrain the BIO5 interpretation.

Strict same-observer strata:
- 106 species;
- BIO5 species-level Wilcoxon p = 0.485;
- conditional-logistic OR = 0.787, p = 0.202.

Observer-balanced broad-range sensitivity:
- 352 species;
- median BIO5 contrast = +0.054 SD;
- Wilcoxon p = 0.075;
- sign p = 0.048.

Geographic white/nonwhite matching:
- 25 km, 246 species: no support;
- 50 km high-resolution WorldClim sensitivity: positive but modest (photo-level p = 0.0266; species-level Wilcoxon p = 0.0405);
- 100 km: no support;
- broader 250 km photo-level contrast remains positive but species-level evidence is weak.

Interpretation: the candidate is **broad-scale temperature-associated sorting**, not demonstrated local heat selection.

### 4. Independent bee breadth — formally negative

Outcome-blind coverage against Noori et al. (2026) curated GloBI v3.1 exact-matched 346/499 third-cohort species.

The prospectively frozen bee-breadth mechanism test retained:
- 113 eligible colour/bee species;
- 48 plant families.

Primary effort-adjusted bee partner breadth:
- Spearman rho = -0.0424;
- one-sided permutation p = 0.3219;
- family-clustered fractional GLM p = 0.847;
- gate = **BEE_BREADTH_NOT_SUPPORTED_UNDER_THIS_TEST**.

This rejects the specific cross-species prediction that broader bee partner sets generally suppress white flower states. It does not reject pollinator-mediated evolution generally.

### 5. Heat versus local bee environment — heat remains the stronger discriminator

A later joint discriminator used spatially localized same-plant bee-interaction information. Because white outcomes and the prior BIO5 result were already known, this is not untouched confirmation.

At the outcome-blind preflight, the smallest radius meeting the coverage rule was 500 km.

Joint conditional models:

- heat + local bee richness:
  - BIO5 OR = 1.138, p = 0.000300;
  - local bee richness OR = 0.946, p = 0.149.

- heat + local Bombus fraction + local bee richness:
  - BIO5 OR = 1.155, p = 0.000266;
  - Bombus fraction OR = 1.019, p = 0.647;
  - local bee richness OR = 0.928, p = 0.0663.

Species-level local-pollinator contrasts are also null.

Thus the current data do not support bee breadth, local bee richness or local Bombus dominance as the common cross-species explanation of the white state. Temperature remains the only tested ecological axis with a reproducible broad-scale signal.

## Current mechanistic model

The strongest model presently allowed by the evidence is:

`
recurrently accessible white/nonwhite phenotype axis
        +
broad-scale temperature-associated ecological sorting
        +
species-specific local ecology/demography
`

This is **not yet** a model of recurrent pigmented -> white evolution, because transition direction has not been established.

## Critical next gate

The active transition-direction analysis is isolated on:

`analysis/white-transition-direction-20260925`

with draft PR #70.

Primary endpoint states are frozen as:
- white fraction >= 0.90: white-dominant;
- white fraction <= 0.10: nonwhite-dominant;
- intermediate species excluded.

The test compares the OpenTree/Grafen CTMC rates:

`
q(nonwhite-dominant -> white-dominant)
`

versus

`
q(white-dominant -> nonwhite-dominant)
`.

Only if this direction gate supports asymmetry can the programme move from “white-associated sorting” toward a comparative claim about recurrent evolution toward white-dominant states.

## Hard boundaries

The current evidence does not establish:
- anthocyanin loss or any molecular loss-of-function mechanism;
- heat-induced whitening;
- nocturnal moth/bat causation;
- absence of morph-specific bee selection within populations;
- an ancestral chromatic state for all focal taxa;
- causal selection rather than correlated biogeographic sorting.
