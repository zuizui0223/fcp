# JBI Chapter 1 — spatial state-distribution protocol

## Central question

The spatial arm no longer asks which environment each flower colour adapted to.
The primary question is:

> When all evaluable present-day flower-colour observations are superimposed globally, is their spatial arrangement random after conditioning on species distributions, or do independent species repeatedly place colour transitions along shared geographic boundaries?

The global map may suppress species identity for visualization, but species identity is never removed from the inferential null model.

## Analysis layers

### Layer 1 — global descriptive map

Plot every evaluable photograph at its observation location, using flower-colour state as the point colour. Species is not foregrounded in the display layer.

This map is descriptive only. It must not be interpreted as evidence of global colour structure because species composition, observation effort, and geographic range boundaries are all spatially heterogeneous.

### Layer 2 — species-conditioned random-labelling test

Observation locations are fixed. Flower-colour labels are permuted strictly within species.

For species i, the null is:

`flower-colour label independent of location | species i`.

This preserves, by construction:

- observation effort and point geometry;
- regions containing many species;
- each species' own range boundary;
- each species' colour-state frequency.

No colour label may move between species.

A global non-randomness statistic must therefore be recomputed for each species-conditioned permutation rather than compared with an unconditional shuffle of all photographs.

## Species-level transition boundaries

Only after the species-conditioned null is rejected do we estimate colour-transition zones separately for each species.

Boundary estimation must distinguish:

- a detected colour boundary;
- a location where a boundary could have been detected but was not;
- a location where the sampling geometry makes boundary detection impossible.

The last state is `not evaluable`, not zero.

## Detectability denominator A(x)

For species i and spatial location/cell x, define `D_i(x)=1` only when the observation geometry provides enough support for a transition boundary to have been detectable there.

Crucially, `D_i(x)` must be defined without using observed flower-colour labels. It may use only quantities such as:

- species-specific observation locations;
- local sample support;
- distances among observations;
- geographic support on both sides of a candidate cell;
- a predeclared interpolation/detection radius.

Then:

`A(x) = sum_i D_i(x)`

is the number of species for which a boundary was detectable at x.

If `B_i(x)=1` denotes a detected boundary for species i, shared-boundary strength is:

`S(x) = sum_i D_i(x) B_i(x) / A(x)`.

Cells with `A(x)=0` are `not evaluable`. The minimum value of A(x) required for the primary shared-boundary surface must be fixed before the confirmatory analysis; low-support cells must be masked rather than interpreted as zero shared-boundary strength.

## Null distribution for shared boundaries

Each permutation replicate must rerun the entire downstream boundary procedure:

1. permute flower-colour labels within each species;
2. reconstruct the species-level colour surface using the frozen procedure;
3. redetect species-level transition boundaries;
4. recompute shared-boundary strength using the unchanged, label-independent detectability mask;
5. recompute the preregistered global summary statistic.

The null must therefore capture boundaries produced accidentally by the observed sampling geometry.

Cell-wise uncorrected p-values are not the primary inference. The confirmatory global statistic should be frozen before evaluation, for example a maximum/envelope statistic or a preregistered summary of high-strength boundary area.

## Geographic correspondence is downstream

Environmental or historical variables are not selected first to predict flower colour.

Shared flower-colour boundaries are discovered first. Only then are they compared against a predeclared geographic reference library. Candidate reference classes may include:

- mountain systems and major topographic breaks;
- dry–wet geographic boundaries;
- coast–interior transitions;
- glacial or late-Quaternary historical boundaries;
- biogeographic regions;
- anthropogenic fragmentation layers.

This second stage is an explanatory correspondence analysis, not the primary discovery test. A boundary with no supported correspondence remains an observed boundary and is reported as unresolved rather than discarded.

## Photograph calibration gate

Current development sample: 6 species × 200 photographs = 1,200 photographs.

The confirmatory colour analysis is blocked until the following calibration gate is completed.

### Calibration set

For each species, 80 photographs are assigned to the blinded calibration set (480 total). The remaining 120 photographs per species (720 total) remain untouched during rule construction.

The 480-image calibration freezes:

1. flower visibility / evaluability;
2. flower-region segmentation rules;
3. species-specific colour-state coding;
4. ambiguity and `unresolved` rules;
5. exclusion rules and failure codes.

The calibration split and all rule versions must be hash-manifested before the 720-image evaluation set is opened.

### Evaluation set

After the rules are frozen, they are applied unchanged to the 720 held-out photographs. A post-hoc change prompted by the evaluation set invalidates that set as confirmatory and must be versioned as a new development cycle.

## Current empirical state

As of 2026-08-28:

- development photographs acquired: 1,200 across 6 species;
- spatial, observer, and seasonal overlap controls: passed;
- flower-colour classifications completed: 0;
- P1 literature records: 26;
- full texts recovered: 11;
- manually extracted records: 4;
- analytically usable continuous boundaries from those four: 0;
- existing test/manifest/citation-contract layer: passing before this branch was cut.

Therefore random versus non-random colour placement and shared boundaries are both currently `not evaluated`.

## Interpretation rule

The governing rule for this chapter is:

> Species may disappear from the map display, but species must never disappear from the null model.

An unconditional pooling of all angiosperm photographs would estimate the global map of species composition, not the spatial organization of flower colour.
