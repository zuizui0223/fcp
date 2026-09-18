# disttrait random-effects comparator diagnostic — 2026-09-19

Status: diagnostic only; the tested DL/Q implementation was **not promoted**.

## Question

Could a simple species-slope random-effects meta-analysis close the v0.7
random-slope/hierarchical comparator gap?

Prototype:

- one OLS signed slope per species;
- per-species sampling variance from the OLS residual mean square;
- DerSimonian-Laird tau^2;
- normal test for the random-effects mean slope;
- Cochran Q test for slope heterogeneity;
- conservative Bonferroni omnibus of mean-slope and heterogeneity p-values.

The prototype was evaluated on the same 24-cell direction-heterogeneity design used
by v0.7, with 40 worlds per cell.

## Diagnostic result

The prototype had the desired power pattern under 50% direction reversal:

- effect 0.4, reversal 0.5:
  - Q heterogeneity detection = **0.95–1.00**;
  - Bonferroni omnibus detection = **0.95–1.00**;
  - common-slope detection = **0.025–0.05**;
  - direction-invariant matched-null detection = **0.40–0.90**.

However, null calibration was unacceptable:

- maximum Q-test null rejection = **0.30**;
- maximum mean-slope null rejection = **0.175**;
- maximum omnibus null rejection = **0.30**;
- the inflation was concentrated in smaller post-missingness samples.

The issue is consistent with using plug-in small-sample slope variances and
large-sample Q/normal reference distributions as though they were exact.

## Decision

Do **not** promote the DL/Q prototype as a disttrait method or comparator.

The v0.8 model-based comparator should instead use an observation-level
random-intercept/random-slope hierarchical model, while retaining explicit
calibration checks on the same frozen synthetic surface.

This diagnostic is a negative methods result and should not be rewritten as
evidence that random-effects models are generally invalid. It only rejects this
specific small-sample DL/Q implementation for the current benchmark.
