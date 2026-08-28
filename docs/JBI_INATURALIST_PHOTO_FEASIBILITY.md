# JBI iNaturalist photo feasibility and calibration gate

## Purpose

This document fixes the operational gate between photograph acquisition and confirmatory global flower-colour spatial inference.

The current development sample contains 6 species with 200 photographs each (1,200 photographs total). Acquisition feasibility and overlap controls have passed, but flower-colour classification has not yet begun. The next valid result is therefore not a spatial boundary result; it is blinded measurement feasibility.

## Split

For each of the six development species:

- calibration: 80 photographs;
- held-out evaluation: 120 photographs.

Totals:

- calibration: 480 photographs;
- held-out evaluation: 720 photographs.

The split must be deterministic, frozen, and hash-manifested before any flower-colour rule is tuned.

## Calibration sequence

Each calibration photograph passes through the following ordered decisions.

### 1. Flower visibility

Assign one of:

- `evaluable`;
- `not_evaluable`.

Failure reasons should be coded explicitly (for example occlusion, distance, blur, overexposure, non-target organ, or insufficient flower area) rather than silently discarded.

### 2. Flower-region segmentation

For evaluable photographs, define the flower region used for colour assessment. A segmentation failure remains a measurement failure and must not be forced into a colour state.

### 3. Species-specific colour coding

Colour-state definitions are calibrated within species. No global binary assumption is imposed. Species may require two or more discrete states, or may remain unsuitable for discrete coding.

### 4. Ambiguity

If the photograph cannot be assigned reliably under the frozen species-specific rules, classify it as `unresolved`. `Unresolved` is not merged with any biological colour state.

## Blindness and leakage control

During calibration, the operator or algorithm used to establish visibility, segmentation, and colour rules must not use downstream boundary outcomes or environmental/geographic explanatory layers to choose thresholds.

The 720-image evaluation set must remain unopened for rule tuning. Any rule changed in response to evaluation-set performance starts a new development version and the affected evaluation set is no longer confirmatory.

## Freeze artifact

Before evaluating the 720 photographs, freeze a manifest containing at least:

- species list;
- photograph identifiers assigned to calibration/evaluation;
- split-generation rule and seed if stochastic;
- visibility codebook version;
- segmentation procedure/version;
- species-specific colour-state codebook version;
- unresolved/not-evaluable rules;
- software commit SHA;
- content hashes of the calibration assignments and codebooks.

## Gate decision

The calibration gate passes only when the six species have a documented measurement rule that is stable enough to apply unchanged to the held-out photographs. A species that cannot support reliable flower visibility, segmentation, or colour-state coding is retained in the audit and marked `not_evaluable` for the relevant downstream analysis rather than being silently removed.

## Downstream constraint

Passing photograph measurement feasibility does not itself imply non-random colour geography. Only after held-out colour states are generated can the species-conditioned random-labelling analysis be run.
