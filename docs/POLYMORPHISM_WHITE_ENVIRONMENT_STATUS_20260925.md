# White-environment mechanism status — 2026-09-25

## Current result

A separate prospective mechanism test was executed on the frozen third-cohort flower-colour measurements. It does **not** modify the current New Phytologist manuscript.

### Primary pre-specified environmental family

Primary panel after response-blind high-clip exclusion and technical availability:
- 12,583 classifiable photographs;
- 281 species with >=5 white and >=5 non-white photographs;
- species are the primary comparative replicate;
- three pre-specified mechanisms with Holm correction.

| Mechanism | Prediction | Species-level result | Photo-level corroboration | Frozen gate |
|---|---|---|---|---|
| BIO5 heat | hotter -> white | median white-minus-nonwhite = +0.069 within-species SD; 57.3% species positive; Wilcoxon p=0.0118; Holm p=0.0354 | OR=1.073 per within-species SD, 95% CI 1.029–1.119, p=0.000919 after near-clip adjustment | **PASS** |
| BIO14 drought relief | wetter driest month -> white | median -0.0187 SD; Holm p=0.692 | OR=0.993, p=0.749 | not supported |
| mean solar radiation | greater radiation -> less white | median -0.00683 SD; Holm p=0.692 | OR=0.987, p=0.544 | not supported |

The BIO5 result is not driven by one or a few species: leave-one-species-out Wilcoxon p ranges from approximately 0.0087 to 0.0152, and a 20,000-bootstrap 95% interval for the median standardized contrast is approximately +0.017 to +0.116.

### Geographic heterogeneity

A post-result geographic diagnostic shows:
- northern-hemisphere species: n=233, median BIO5 contrast +0.078 SD, Wilcoxon p≈0.010;
- southern-hemisphere species: n=48, median +0.029 SD, p≈0.749;
- tropical species: n=17, no positive signal;
- mid-latitude species: n=212, median +0.084 SD, p≈0.012;
- high-latitude species: n=52, same positive median direction but low precision.

Therefore the current signal should be described as a **temperate/broad-scale heat-sorting candidate**, not a universal global law.

## Observer diagnostics

Two post-hoc diagnostics address observer structure.

### Strict paired-observer test

Only observers who photographed both white and non-white flowers of the same species contribute.

For BIO5:
- 144 paired species-observer strata;
- 106 species;
- species-level median paired difference = 0;
- Wilcoxon p=0.485;
- conditional-logistic OR=0.787, 95% CI 0.544–1.138, p=0.202.

This does not retain the primary heat signal.

However, this restriction changes the spatial estimand sharply. The paired-observer strata span a median of about 56 km, whereas the full eligible species panels span a median of about 3,320 km. Thus this test mainly asks about local within-observer sorting, not the broad geographic environmental sorting targeted by the WorldClim test.

### Observer-balanced broad-range test

Each observer is first averaged within species x colour, so an observer cannot gain extra weight from contributing more photographs. Between-observer geography is retained.

For BIO5:
- 352 species with >=3 observers per colour;
- median contrast +0.054 SD;
- 55.4% positive;
- Wilcoxon p=0.0750;
- two-sided sign p=0.0484.

The direction persists but statistical strength is weaker.

BIO14 and solar radiation remain null under observer balancing.

## Current ecological interpretation

The evidence now separates three propositions:

1. **Broad geographic temperature sorting is a real candidate.** White states occur modestly more often toward warmer parts of species' sampled ranges, and this was the only pre-specified environmental mechanism to clear the frozen gate.
2. **The effect is weak and scale-dependent.** It is not recovered when inference is restricted to the much smaller spatial ranges covered by observers who photographed both states.
3. **Drought relief and chronic solar-radiation filtering are not supported as cross-species common explanations under these tests.**

The appropriate claim is therefore:

> Across the third-cohort species, white flower states show a weak recurrent tendency to occupy warmer portions of species' sampled climatic ranges. The association survives exposure control and unequal observer weighting, but not a strict same-observer/local-scale restriction, so it is evidence for broad-scale environmental sorting rather than demonstrated causal heat-driven whitening.

## What remains unresolved

This result does not establish:
- pigmented -> white evolutionary direction;
- heat-induced pigment loss rather than colonization history, correlated habitat, or phenology;
- a genetic anthocyanin-loss mechanism;
- pollinator causation;
- external replication outside iNaturalist.

The next high-value discriminator is **pollinator turnover versus temperature** using an independent pollinator dataset: test whether the BIO5 contrast persists after a species-range pollinator-guild predictor is added, and whether white is specifically associated with nocturnal pollination rather than temperature alone.
