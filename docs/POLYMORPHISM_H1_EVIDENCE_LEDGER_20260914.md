# Polymorphism H1 evidence ledger — observer-disjoint D reliability

Date reconciled: 2026-09-14 JST

## Canonical status

**H1 measurement-validity status: supported under the first-frozen repeated-partition protocol, with a later stricter single-split stress-test caveat.**

Canonical label for manuscript routing:

`H1_MEASUREMENT_VALIDITY_SUPPORTED_WITH_STRICT_STRESS_TEST_CAVEAT`

This ledger resolves two apparently conflicting observer-disjoint analyses by their prospective chronology. It does not add a new statistical test and does not alter either frozen result.

## Outcome definition

Both analyses concern the same continuous four-state flower-colour polymorphism score

`D = 1 - sum_k p_k^2`

using the biological states:

- white
- yellow_orange
- red_pink
- blue_purple

`mixed_uncertain` is not promoted to a biological state.

H1 is a measurement/sampling-stability gate. It is not a global polymorphism-prevalence estimator and is not itself a biological mechanism test.

## 1. First-frozen primary evidence: repeated observer-disjoint partitions

Protocol:
`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`

Protocol freeze commit:
`92488ada5f4f535f6731bba6519101d8ed8a7e1a`

Protocol freeze time:
2026-09-12 17:03:49 UTC

Result:
`results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`

Result commit:
`927940e5a3b0541b96515b53107cab24cb138950`

Result-open time:
2026-09-12 17:15:12 UTC

The protocol used 200 observer-disjoint partitions. Its frozen primary rule required, on reserve, at least 100 median paired species, median split Spearman rho >= 2/3, and 5th-percentile split rho >= 0.5. Sensitivity analyses could not rescue the primary rule.

### Fresh reserve result

- primary half minimum classifiable count: 20
- partitions: 200
- paired species: median 329, range 315–339
- median Spearman rho: **0.789102993**
- 5th percentile rho: **0.765164994**
- 95th percentile rho: 0.810940258
- median Lin CCC: **0.854830674**
- 5th percentile CCC: 0.833728332
- median absolute D difference: 0.078980607
- 95th percentile absolute bias: 0.010581315
- median Spearman-Brown reliability: 0.882121371

Frozen verdict:
`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`

Discovery was consistent with the reserve result but is not needed to rescue reserve.

## 2. Later stricter deterministic single-split stress test

Protocol:
`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_RELIABILITY_PROTOCOL_20260913.md`

Protocol freeze commit:
`8b2b416265b2ddcad53ecaafc214c80cc16a2648`

Protocol freeze time:
2026-09-13 04:03:25 UTC

Result:
`results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json`

Result commit:
`9a507a1a0cfbb5316657760c1283865a672ad346`

Result-open time:
2026-09-13 12:36:06 UTC

This protocol was frozen before its own deterministic split outcomes, but **after the first-frozen repeated-partition H1 result had already been opened**. It therefore cannot retroactively replace the earlier test as an untouched primary H1 test. It is retained as a valid, deliberately stricter post-outcome stress test.

### Fresh reserve result

- reliability-eligible species: 363
- observer leakage count: 0
- Spearman D_A vs D_B: **0.792727693**
- bootstrap 95% interval: **0.741867914–0.832365877**
- Lin CCC: **0.847429024**
- Pearson r: 0.847883340
- median absolute D difference: 0.063775510
- 90th-percentile absolute D difference: 0.188166183
- H1 0.10 gate agreement: **0.848484848**
- H1 0.20 gate agreement: **0.878787879**

The deterministic stress test required reserve Spearman rho >= 0.80. The observed value, 0.792727693, missed that floor. CCC >= 0.75 and bootstrap lower bound >= 0.70 both passed.

Frozen verdict:
`H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`

No threshold is relaxed and this failure is preserved.

## Reconciliation rule

The two results answer related but not identical reliability questions:

1. The first-frozen test asks whether D is stably reproducible across the distribution of observer-disjoint partitions under a prespecified minimum reliability rule.
2. The later stress test asks whether one deterministic split satisfies a stronger rho >= 0.80 requirement.

Because the first protocol and result were frozen/opened before the second protocol existed, the first is the canonical primary evidence. The later failure constrains the strength of the claim rather than reversing the prospective ordering.

Therefore:

- **Supported:** D has observer-disjoint reproducibility sufficient to pass the first-frozen H1 measurement-validity gate.
- **Also true:** a later, stricter deterministic split missed a rho=0.80 criterion by 0.00727.
- **Not supported:** near-perfect reliability, split-invariant rho >= 0.80, or any claim that every reasonable observer split passes an 0.80 floor.

## Historical ~0.971 value

The historical approximately 0.971 observer-reproducibility value is **not used as evidence for current four-state D reliability**. The dedicated deterministic result explicitly records `historical_approx_0_971_used_as_evidence = false`.

The manuscript must not state that current four-state D has observer-disjoint reliability of ~0.971.

## Allowed manuscript claim

> The continuous four-state polymorphism score D was reproducible across repeated observer-disjoint partitions under the first-frozen validation rule in the species-disjoint reserve cohort (median split rho 0.789; 5th percentile 0.765; median CCC 0.855). A later deliberately stricter deterministic split yielded rho 0.793 and failed a prespecified 0.80 floor, so near-perfect or split-invariant reliability is not claimed.

## Hard nonclaims

- no global prevalence estimate from the 369/363 high-depth cohorts
- no globally representative-species claim
- no current-D observer reliability of ~0.971
- no near-perfect or split-invariant reliability claim
- no post-hoc relaxation of the rho=0.80 stress-test floor
- no promotion of `mixed_uncertain` to a biological morph
- no causal ecological, phylogenetic, or H2 geometry inference from H1

## Consequence for paper architecture

H1 remains a **measurement-validity/admission gate**, not the paper's main biological result. H2 white-versus-nonwhite colour-space geometry remains the main positive biological result. H3a and H3b remain bounded fresh-reserve negatives and are not used to rescue H1 or H2.