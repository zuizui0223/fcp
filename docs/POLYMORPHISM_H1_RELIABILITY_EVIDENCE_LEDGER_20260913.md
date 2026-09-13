# H1 reliability evidence ledger

Date frozen: 2026-09-13 JST

## Purpose

Resolve the coexistence of two observer-disjoint reliability analyses for the current four-state species-level polymorphism score

`D = 1 - sum_k p_k^2`

without post-outcome choice of the more favorable rule.

The governing principle is **freeze chronology**: the earliest outcome-blind protocol for this endpoint remains the primary inferential contract. Any later protocol is a sensitivity/audit and cannot replace or rescue that primary result.

## 1. Canonical primary analysis — 200 observer partitions

### Protocol freeze

- protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`
- freeze commit: `92488ada5f4f535f6731bba6519101d8ed8a7e1a`
- commit time: `2026-09-12T17:03:49Z`
- commit message: `freeze H1 observer-disjoint D reliability protocol`

The protocol fixes 200 deterministic outcome-blind observer partitions, a >=20-classifiable-per-half primary gate, reserve as the decision cohort, and the following support rule:

1. median paired species >= 100;
2. median split-half Spearman rho >= 2/3;
3. 5th percentile split-half Spearman rho >= 0.50.

Discovery cannot rescue reserve failure. The >=15-per-half analysis is sensitivity only.

### Implementation-only optimization after freeze

- workflow head commit: `5142f7951af0dde5364bb047a566d67e8c479e51`
- commit message: `optimize H1 observer split runtime without changing frozen estimand`

The change preaggregates observer-level row and four-state counts but retains the frozen observer assignment, state definition, half-size gates, statistics, seeds, and decision rule. It is computational rather than inferential.

### Canonical execution

- workflow run: `34707537360`
- run conclusion: success
- source SHA verification: success
- syntax check: success
- frozen analysis: success
- decision receipt verification: success
- artifact: `10302466831`
- artifact name: `polymorphism-h1-observer-disjoint-reliability-20260913`
- artifact digest: `sha256:19185ad9160b24ba47ef8d99a33d952ab5f5b6b62243ff756d41917bfc587c53`

### Frozen result

Canonical machine-readable result:
`results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`

Canonical result freeze:
`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_RESULT_FREEZE_20260913.md`

Result-freeze commit:
`902768972bdc4855220346c3250d5b94703617ab`

Reserve primary, across 200 partitions:

- full eligible species = 363;
- partitions with defined rho = 200/200;
- median paired species = **329**;
- median rho = **0.789102993**;
- rho 5th percentile = **0.765164994**;
- rho 95th percentile = **0.810940258**;
- median Spearman-Brown projected full-estimate reliability = **0.882121371**;
- median CCC = **0.854830674**;
- median absolute split difference in D = **0.078980607**;
- median signed bias A-B = **-0.000979831**.

All three frozen primary support conditions pass.

**Canonical verdict: `H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED+DISCOVERY_CONSISTENT`.**

Discovery is consistent but non-decisive:

- full eligible species = 369;
- median paired species = 340;
- median rho = 0.811544534;
- rho 5th percentile = 0.788380274;
- median CCC = 0.851426192.

## 2. Later deterministic single-split audit — sensitivity only

### Later protocol freeze

- protocol: `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_RELIABILITY_PROTOCOL_20260913.md`
- freeze commit: `8b2b416265b2ddcad53ecaafc214c80cc16a2648`
- commit time: `2026-09-13T04:03:25Z`
- commit message: `freeze observer-disjoint D reliability protocol before split outcomes`

This protocol was frozen **after the canonical 200-partition primary protocol and after the canonical primary result had already been produced**. It therefore cannot become the primary inferential rule for H1.

Its stricter engineering-style reserve support rule requires simultaneously:

- one deterministic split Spearman rho >= 0.80;
- bootstrap 95% lower bound >= 0.70;
- Lin CCC >= 0.75.

### Execution

- workflow run: `34737156300`
- artifact: `10311386392`
- artifact digest: `sha256:ded4d8ae15546321a03a54389abd52cf746464f8cf30037332e9a78ef0eb7b81`

A first run failed before split-D outcomes because pandas 3 rejected partial string assignment into an integer observer column. The repair changed only dtype handling; source SHA and syntax gates had passed before the failure.

Successful single-split audit result:

Discovery:

- n = 369;
- rho = 0.8205485;
- bootstrap 95% CI = [0.7774661, 0.8548739];
- CCC = 0.8542231.

Reserve:

- n = 363;
- rho = **0.7927277**;
- bootstrap 95% CI = **[0.7418679, 0.8323659]**;
- CCC = **0.8474290**;
- median |D_A-D_B| = 0.0637755.

Under this later protocol, the arbitrary rho >= 0.80 floor is missed narrowly while the CI and CCC floors pass. Its internal verdict is therefore `H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`.

### Governance interpretation

This later result is **not a contradictory primary replication**. It is one deterministic observer partition analyzed under a later and stricter cutoff, whereas the frozen primary estimand is the distribution of reliability across 200 outcome-blind observer partitions.

The point estimate itself is highly consistent with the canonical reserve distribution:

- later single split rho = 0.79273;
- canonical 200-partition median rho = 0.78910;
- canonical 5th–95th percentile range = 0.76516–0.81094.

Therefore the later audit is retained as evidence that the quantitative reliability estimate is stable to an alternative deterministic split, while its post-primary 0.80 pass/fail threshold is **not allowed to overturn the earlier frozen decision rule**.

## 3. Historical ~0.971 value

A previously circulated approximate observer-disjoint reproducibility value of ~0.971 was searched in the repository, but its provenance as this exact current four-state D endpoint could not be established before the direct H1 analyses were opened.

It is therefore **not used as canonical H1 evidence** and should not appear as the principal observer-disjoint reliability estimate in the manuscript.

## 4. Allowed H1 claim

Within the existing high-depth discovery/reserve design, the continuous four-state flower-colour polymorphism score D is reproducibly recoverable from observation sets contributed by disjoint observers. The species-disjoint reserve cohort passes the earliest pre-outcome reliability rule across 200 outcome-blind observer partitions, and a later deterministic split yields a closely matching rho (~0.793).

## 5. Hard nonclaims

H1 reliability does not establish:

- global prevalence of flower-colour polymorphism;
- representativeness of the 369/363 high-depth cohorts;
- correctness of every image-level colour label;
- absence of geographic or temporal morph-frequency structure;
- ecological, phylogenetic, genetic, or adaptive causes of D;
- H2 white-versus-nonwhite geometry;
- a universal biological meaning for any reliability cutoff.

## 6. Manuscript status

H1 is **closed as supported measurement validity** under the canonical chronological evidence hierarchy.

Do not reopen H1 by threshold hunting, alternative split selection, or choosing between the two protocols according to their pass/fail labels. The next inferential work belongs to the prospective expansion of the already frozen H2 white-axis target, not further reliability-rule optimization.
