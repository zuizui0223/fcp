# H1 observer-disjoint D reliability — frozen result

Date frozen: 2026-09-13 JST

Protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_RELIABILITY_PROTOCOL_20260913.md`

Workflow run: `34737156300` (run #2; success)

Artifact: `10311386392` (`polymorphism-h1-observer-disjoint-D-reliability-20260913`)

Artifact digest: `sha256:ded4d8ae15546321a03a54389abd52cf746464f8cf30037332e9a78ef0eb7b81`

Canonical machine-readable receipt: `results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json`

## Status

**Frozen verdict: `H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`.**

The current four-state polymorphism score is

`D = 1 - sum_k p_k^2`,

using only `white`, `yellow_orange`, `red_pink`, and `blue_purple`; `mixed_uncertain` and unresolved technical states are excluded.

The historical approximate observer-reproducibility value (~0.971) is not used as evidence for this endpoint. The present deterministic observer-disjoint recomputation is the canonical direct test.

## Frozen decision rule

The species-disjoint reserve cohort had to satisfy all three preregistered engineering-style floors:

1. Spearman `rho(D_A, D_B) >= 0.80`;
2. 95% species-bootstrap lower bound `>= 0.70`;
3. Lin CCC `>= 0.75`.

Discovery cannot rescue reserve failure.

## Result

### Discovery calibration

- reliability-eligible species: **369**
- observer leakage: **0**
- Spearman rho: **0.8205485**
- bootstrap 95% CI: **[0.7774661, 0.8548739]**
- Lin CCC: **0.8542231**
- Pearson r: **0.8547401**
- median `|D_A-D_B|`: **0.0603567**
- 90th percentile `|D_A-D_B|`: **0.1832545**
- mean signed difference `D_A-D_B`: **0.0006442**

### Species-disjoint reserve decision cohort

- reliability-eligible species: **363**
- observer leakage: **0**
- Spearman rho: **0.7927277**
- bootstrap 95% CI: **[0.7418679, 0.8323659]**
- Lin CCC: **0.8474290**
- Pearson r: **0.8478833**
- median `|D_A-D_B|`: **0.0637755**
- 90th percentile `|D_A-D_B|`: **0.1881662**
- mean signed difference `D_A-D_B`: **0.0065950**

The bootstrap-lower-bound and CCC floors passed. The reserve Spearman floor failed narrowly: `0.7927277 < 0.80`, a shortfall of approximately `0.00727`. Under the frozen all-conditions rule this is a formal failure and is not rounded upward or rescued by the other diagnostics.

## H1 gate-stability diagnostics

These are diagnostics only; they do not override the D-reliability verdict.

Reserve half-vs-half agreement:

- second-state fraction `>=0.10`: agreement **0.848485**, Cohen kappa **0.69336**
- second-state fraction `>=0.20`: agreement **0.878788**, Cohen kappa **0.68201**

Discovery agreement was 0.81030 / 0.87534 for the same primary / strict gates.

## Interpretation

The result does **not** show that D is noise. Observer-disjoint estimates have substantial rank and absolute agreement in both cohorts, and no observer leakage occurred. What failed is the stronger frozen claim that D clears all prespecified engineering floors for a highly stable species-level measurement.

Therefore the manuscript must not state that H1 measurement reliability was formally validated. It may report the observed reliability estimates and the narrowly missed prespecified gate transparently.

This H1 failure does not retroactively change the already frozen H2 construction or its structured-null tests. H2 is based on full-cohort admission plus label-free continuous palette geometry, not on selecting species by observer-split D agreement. However, H2 must not be rhetorically promoted by claiming that continuous D itself has passed the H1 reliability gate.

## Hard nonclaims

- no global prevalence estimate for the 42,111-species frame;
- no claim that excluded species are monomorphic or intrinsically unreliable;
- no validation of `mixed_uncertain` as a biological morph;
- no resurrection of the historical ~0.971 figure as canonical evidence;
- no post-outcome lowering of the `rho >= 0.80` floor.
