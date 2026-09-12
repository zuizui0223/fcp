# H3b reserve sampled-span replication — frozen result

Date frozen: 2026-09-12 JST

Protocol: `docs/POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md`

Workflow run: `34677468362`
Artifact: `10292399238` (`polymorphism-h3b-reserve-span-20260912`)
Artifact digest: `sha256:34c5e646725f6865e313b17f1b70f2471db8169f443fb0649da83044d654386c`

## Status

**Frozen verdict: `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.**

This is a species-disjoint fresh replication of the discovery-selected association between sampled geographic span and the continuous four-state polymorphism score D.

## Discovery calibration

- n = 369 species
- Spearman rho(D, log1p sampled span) = 0.1798786
- two-sided 20,000-permutation p = 0.00089996
- unbiased-D sensitivity: rho = 0.1796832, p = 0.00079996
- observer/classifiable-adjusted partial rank association: rho = 0.1585085, p = 0.00269987

These values are calibration/selection history and are not independent confirmation.

## Fresh reserve replication

- n = 363 species
- Spearman rho(D, log1p sampled span) = -0.0025855
- two-sided 20,000-permutation p = 0.9586021
- unbiased-D sensitivity: rho = -0.0023742, p = 0.9629019
- observer/classifiable-adjusted partial rank association: rho = 0.0055187, p = 0.9162042

The direction and magnitude do not replicate.

## Phylogenetic sensitivity

The prespecified rank-PGLS sensitivity uses the frozen S1-S3 V.PhyloMaker2 trees. Each reserve tree retains 341 tips.

- S1: beta_span_rank = -0.0055931, p = 0.9136754, lambda ~= 1e-7
- S2: beta_span_rank = -0.0055931, p = 0.9136754, lambda ~= 1e-7
- S3: beta_span_rank = -0.0055931, p = 0.9136754, lambda ~= 1e-7

No placement scenario supports the sampled-span association.

## Interpretation

The discovery association between sampled geographic span and D does not generalize to the species-disjoint reserve cohort. It must not be used as a general predictor of flower-colour polymorphism in the main paper.

This result does **not** imply that true geographic range size is unrelated to polymorphism. The frozen predictor is sampled span from the fixed measured-photo design, not a biological range-size estimate.

## Consequence for manuscript architecture

H3b sampled-span is closed as a non-replicated discovery result. The mainline should not be rescued by threshold changes, alternative span definitions, post-hoc trait selection, or additional H3b predictor hunting. The paper should return to the measurement-validity and colour-space geometry claims, with H3a/H3b reported as bounded negative results where useful.
