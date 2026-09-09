# RGFCA sharedness-v2: frozen qualification failure and post-fail identifiability diagnosis

Updated 2026-09-08 (Japan).

## Decision

The high-depth sharedness-v2 lane is **closed as a failed synthetic qualification**. The frozen geometry contains **150 species × 300 photographs = 45,000 photographs**, split species-disjoint into 75 training and 75 evaluation species. The paired continuous estimator compares flower-colour distance with matched background-colour distance and asks whether a geographic boundary learned from training species predicts held-out species.

The prospectively fixed synthetic qualification did **not** pass. Nuisance-world false-positive control was acceptable (maximum approximately **2.4%**), but all three prespecified hard-positive power gates failed. The best hard-positive detection rate was only approximately **3.2%**, against the required **80%** floor. Therefore:

- `qualification_pass = false`;
- empirical image acquisition for this lane remains prohibited;
- observed flower colour, observed background colour and image pixels remain unopened for the v2 empirical test;
- no statistic, threshold, axis grid, species count, photo count, retention arm or power floor may be changed post hoc to rescue the decision.

Canonical result: `docs/supporting/rgfca_sharedness_v2_synthetic_qualification_result_v1.json`.

## Post-fail identifiability audit

A separate descriptive audit was run after the failure without reopening qualification or searching parameters. It reproduced the frozen positive generator and asked whether its true shared boundary was actually observable on the retained real metadata geometry.

Across the six audited hard-positive scenario/threshold-SD combinations (1,500 synthetic worlds), only about **42–45%** of shared species straddled the true injected boundary. Only about **22–26%** had at least 10% of retained observations on the minority side, and about **15–18%** had at least 20% on the minority side. At the species-instance level the median minority-side fraction was **0**.

The fixed geographic axis grid was not the dominant limitation. The nearest frozen axis was typically only about **6–7 degrees** from the true common normal, and the median nearest-axis plus frozen-threshold partition agreement was **1.0**.

Canonical diagnostic: `docs/supporting/rgfca_sharedness_v2_postfail_identifiability_v1.json`.

## Interpretation

The failed power gate should not be described as evidence that shared flower-colour geography is absent in nature. The diagnostic instead shows a **support/identifiability mismatch for the global common-hyperplane estimand on the realized species-range geometry**: in more than half of shared-species instances, the sampled species distribution does not cross the injected common boundary at all, so a within-species contrast cannot identify that boundary regardless of nominal sample depth.

This diagnosis occurs upstream of biological measurement. It does not convert the failed qualification to a pass, authorize empirical colour acquisition, or justify a new post-hoc tuning cycle. The scientifically admissible conclusion is narrower: **this frozen global-common-boundary design is not sufficiently qualified on the available species-range geometry to support the intended empirical sharedness claim.**

## Relation to existing RGFCA evidence

This result strengthens, rather than replaces, the existing separation of claims:

1. a weak within-species photo-derived distance-colour association exists in discovery and replicates directionally in the reserve;
2. the reserve matched flower-minus-background robustness gate is unsupported (`p = 0.087`), so flower-specific replication is not established;
3. the original repeated-field G1 primary test (`p = 0.070`) and species-disjoint commonness test (`p = 0.856`) do not establish a shared global boundary;
4. sharedness-v2 cannot presently adjudicate that biological question because the intended common-hyperplane estimand fails its prospective synthetic qualification on the real metadata geometry.

No biological absence claim follows from item 4.