# FCP white-transition direction protocol — 2026-09-25

Status: **prospectively frozen before opening phylogenetic transition-rate results**.

This analysis is separate from the frozen New Phytologist manuscript and from the white-environment mechanism analysis. It addresses a prerequisite for mechanistic language: whether the recurrent white/non-white axis shows macroevolutionary asymmetry consistent with repeated transitions toward white-dominant states.

## Biological question

The frozen H2 result is sign-invariant and therefore establishes a recurrent white-versus-nonwhite axis, not evolutionary direction. This analysis asks whether, across the third-cohort species tree, transitions between endpoint-dominant colour states are asymmetric.

The states are defined before phylogenetic inference:

- `nonwhite_dominant`: white fraction <= 0.10;
- `white_dominant`: white fraction >= 0.90;
- intermediate species are excluded from the directional endpoint analysis.

Only species with >=40 globally classifiable rows in the four frozen biological colour states are eligible.

## Data source

Biological state source:
- immutable third-cohort prospective measurement run `35177668182`;
- same frozen white/non-white coarse-state construction used by the polymorphism programme.

Phylogeny source:
- Open Tree of Life taxonomy and induced topology queried after state definitions are frozen;
- only exact, non-approximate TNRS matches are used;
- duplicated OTT IDs are excluded rather than merged.

## Tree uncertainty

The induced OpenTree topology is not time-scaled and may contain polytomies.

For 100 pre-specified replicates:
1. unresolved polytomies are randomly resolved with `ape::multi2di(random=TRUE)`;
2. branch lengths are assigned with Grafen's method, power = 1;
3. the same endpoint states are analyzed on each resolved tree.

Seed: `20260925`.

## Primary model

For each resolved tree, fit binary continuous-time Markov models with `ape::ace(type="discrete")`:

- ER: one shared transition rate;
- ARD: separate rates
  - `q_NW = q(nonwhite_dominant -> white_dominant)`;
  - `q_WN = q(white_dominant -> nonwhite_dominant)`.

Primary directional statistic:

`
R = q_NW / q_WN
`

and equivalently `log(R)`.

The likelihood-ratio statistic comparing ARD to ER has 1 degree of freedom.

## Estimability gate

The analysis is `not_estimable` unless:
- at least 150 endpoint-dominant species remain on the induced tree;
- both endpoint states contain at least 50 tree tips;
- at least 80 of 100 tree-resolution replicates complete both ER and ARD fits.

## Direction gate

The result is labeled `NONWHITE_TO_WHITE_ASYMMETRY_SUPPORTED` only if all hold:

1. median `R > 1`;
2. at least 90% of completed tree resolutions have `R > 1`;
3. median ARD-vs-ER likelihood-ratio p < 0.05;
4. at least 80% of completed tree resolutions have likelihood-ratio p < 0.05.

If estimable but any criterion fails, the result is `DIRECTIONAL_ASYMMETRY_NOT_SUPPORTED_UNDER_THIS_TEST`.

## Sensitivity

A threshold sensitivity is descriptive only and cannot rescue the primary gate:
- 0.05 / 0.95 endpoint thresholds, if each state retains >=30 tree tips.

No alternative threshold is searched after the result is opened.

## Hard nonclaims

Even a positive direction gate does not establish:
- the ancestral flower colour of angiosperms or of all focal lineages;
- within-population mutation direction;
- anthocyanin loss as the molecular mechanism;
- heat, pollinator, or other ecological causation;
- absolute evolutionary rates, because OpenTree + Grafen branch lengths are not a dated phylogeny.

A positive result supports only a comparative macroevolutionary asymmetry among endpoint-dominant states under this topology-based design.
