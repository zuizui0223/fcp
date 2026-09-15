# H1 observer-disjoint D — chronology audit

Date: 2026-09-15 JST

This audit runs **no new statistical test**. Its purpose is to reconstruct the immutable order of the H1 observer-disjoint analyses from Git commits and GitHub Actions records, because the 2026-09-14 reconciliation receipt incorrectly identified a later repeated-partition analysis as the first frozen H1 test.

## 1. Immutable chronology

### A. Outcome-blind deterministic split: first frozen H1 test

1. `637d204fa09d21eafc480543953a87f0b4d754ce` — `implement outcome-blind H1 observer split preflight`, 2026-09-12 16:53:10 UTC.
2. Workflow run `34706579235` — `polymorphism H1 observer-split preflight`, 16:53:13–16:53:41 UTC, **success**. The preflight reads only `species`, `observer_id`, and `global_classifiable`; it does not open `morph`, palette variables, D, or D_unbiased.
3. `743e13d3597c88ea7a725c596ea733b30d978024` — `freeze H1 observer-disjoint D protocol before outcome opening`, 16:53:54 UTC.
4. `a77e23bd7a098307c4961d77271fb2e5a294ce88` — `freeze H1 outcome-blind observer split receipt`, 16:54:57 UTC.
5. `daabd3b63d540238a2cd27185330a59ca45e9d2b` — `implement frozen H1 observer-disjoint D test`, 16:55:56 UTC.
6. Workflow run `34706708406` — `polymorphism H1 observer-disjoint D`, 16:55:59–16:56:43 UTC, **success**. The job log records the primary reserve result before any of the later H1 protocols were frozen.

Primary immutable artifact from run `34706708406`:

- artifact ID: `10302770009`
- artifact name: `polymorphism-h1-observer-disjoint-d-20260913`
- digest: `sha256:7a87432d91c257a1c010b9cf075eeb8a5f15e2ea016d55ed6d8cfaddc6f0b64e`

Primary reserve result printed in that run:

- n = 363 species
- Spearman rho(D_A, D_B) = **0.8109164414604769**
- bootstrap 95% percentile interval = **[0.7648383249147064, 0.8474976906537689]**
- permutation p = **4.999750012499375e-05**
- Lin CCC = **0.8582944313924915**
- median |D_A-D_B| = **0.06243496357960443**

Frozen decision gates were reserve rho >= 0.80, bootstrap lower bound > 0.70, and permutation p < 0.001. All pass.

The result was committed into the repository later, but the outcome itself is timestamped and immutable in the successful 2026-09-12 GitHub Actions artifact/log. Repository commit time therefore must not be substituted for outcome-opening time.

### B. Repeated-partition analysis: later robustness analysis

The repeated 200-partition protocol was frozen **after** the deterministic primary outcome was already known:

- protocol commit: `92488ada5f4f535f6731bba6519101d8ed8a7e1a`, 2026-09-12 17:03:49 UTC
- result commit: `927940e5a3b0541b96515b53107cab24cb138950`, 17:15:12 UTC

Reserve across 200 observer-disjoint partitions:

- median paired n = 329
- median rho = **0.7891029925726101**
- q05 rho = **0.7651649940131238**
- q95 rho = **0.8109402582919217**
- median CCC = **0.8548306738463054**

This analysis supports substantial reproducibility across many alternative observer partitions, but it is a later robustness analysis and cannot replace the first-frozen deterministic primary by prospective chronology.

### C. Later strict deterministic stress test

A still later strict test was frozen after the primary H1 result was known:

- protocol commit: `8b2b416265b2ddcad53ecaafc214c80cc16a2648`
- protocol time: 2026-09-13 04:03:25 UTC
- result commit: `9a507a1a0cfbb5316657760c1283865a672ad346`
- result time: 2026-09-13 12:36:06 UTC

Reserve:

- n = 363
- rho = **0.7927276933213221**
- bootstrap 95% interval = **[0.7418679143323839, 0.832365877315184]**
- CCC = **0.8474290237389148**
- its deliberately strict rho >= 0.80 criterion is not met.

This failure is retained as a stress-test caveat. Its threshold is not relaxed, but because the protocol is later than the first-frozen primary, it does not retroactively replace the primary decision.

## 2. Correction to the 2026-09-14 reconciliation

`results/polymorphism_h1_evidence_reconciliation_20260914/result.json` previously stated that the repeated 200-partition analysis was the first frozen H1 primary. That chronology is false: the deterministic protocol and successful outcome run precede the repeated-partition protocol by about ten minutes.

The reconciliation is therefore corrected on chronology only. No H1 statistic, threshold, split assignment, or result is recomputed or changed.

## 3. Canonical interpretation

H1 is positive as a **measurement-validity** result under the first-frozen deterministic observer-disjoint design. The later robustness analyses show that the exact correlation is split-design-sensitive at roughly rho ~0.79–0.81, so the paper must not claim near-perfect or split-invariant reliability.

Allowed wording:

> Under the frozen high-depth photo design, continuous four-state flower-colour polymorphism D is reproducible at the species level across disjoint observer sets in the species-disjoint reserve cohort (first-frozen reserve rho = 0.811, 95% bootstrap interval 0.765–0.847). Later alternative observer-partition analyses yielded rho around 0.79, so near-perfect or split-invariant reliability is not claimed.

The historical project-summary value near 0.971 is not current four-state-D evidence and remains excluded.
