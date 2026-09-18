# Flower-colour polymorphism Supporting Information map — 2026-09-18

This document is the compact Supporting Information index for the active flower-colour polymorphism paper.

It does not recompute any result. It organizes frozen protocols, machine-readable outputs, figures and provenance into a manuscript-facing evidence map.

## S1. Global frame, cohort roles, acquisition lineage and inferential separation

Global programme/frame provenance:

- `docs/RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md` — conceptual relationship between the upstream Repeated Global Flower-Colour Atlas and the current species-level polymorphism estimands;
- `docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md` — metadata-discovery and 42,111-species opportunity-frame lineage
- metadata-discovered V1 + cache-resistant V2 union = **42,111 unique iNaturalist species**;
- discovery grid = 18 × 9 equal-area cells;
- V2 = 20 metadata-only rounds × 162 cells = **3,240** fixed request attempts;
- V2 request errors = 0;
- no image pixels or flower-colour outcomes used to define the species union;
- high-depth capacity after observer cap: `U100 = 4,730 species`.

The 42,111 frame is an outcome-blind metadata-discovery/opportunity universe, not an estimate of all angiosperm species and not a polymorphism-prevalence denominator.

### Cohort roles

| Cohort / execution | Upstream source | Analyses using it | Biological outcome status | Key denominator |
|---|---|---|---|---:|
| Discovery | Original RGFCA high-depth iNaturalist measurement table; 500 species × 100 raw photos | D; H1 diagnostic; legacy H2 discovery/audit; D–spatial discovery; H3b discovery calibration | opened legacy cohort | 369 D-eligible species |
| Reserve | Species-disjoint RGFCA complement under the same acquisition contract; 500 species × 100 raw photos | H1 primary; legacy H2 validation; D–spatial replication; H3a; H3b reserve replication | opened species-disjoint reserve | 363 D-eligible species; 341 H3a tips |
| P500 | Separate prospective high-depth expansion from the 42,111-species opportunity frame | measurement transport only | measurement PASS; no durable H2 verdict | 499 species, 49,900 rows; 373 measurement-evaluable |
| Third cohort | New species and fresh photo IDs after exclusion of legacy 1,000 species and P500 | untouched prospective q_white/W test | H2 confirmed | 499 species, 49,900 rows; 377 measurement-evaluable |

The original RGFCA discovery/reserve acquisition contract required Research Grade species-rank iNaturalist records with photographs and georeferences, flowering annotation term 12/value 13, positional accuracy <=5 km, unobscured/open coordinates and allowed CC licences. Observer contribution was capped at two photographs per species and deterministic geographic maximin sampling fixed 100 raw photographs per species. No native-range restriction or explicit captive/wild filter was imposed.

Discovery and reserve therefore reuse the same RGFCA measurement resource but answer different species-level polymorphism questions. The third cohort is species- and photo-disjoint from the legacy cohorts where specified, but all cohorts remain within the same iNaturalist source/opportunity universe unless explicitly stated otherwise.

## S2. H1 observer-disjoint measurement validity

Canonical protocol:

- `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`

Canonical machine-readable result:

- `results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`

Frozen result note:

- `docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_RESULT_FREEZE_20260913.md`

Workflow provenance:

- run: `34707537360`
- artifact ID: `10302466831`
- artifact digest: `sha256:19185ad9160b24ba47ef8d99a33d952ab5f5b6b62243ff756d41917bfc587c53`

### Primary reserve result

| Quantity | Frozen value |
|---|---:|
| Partitions | 200 |
| Median paired species | 329 |
| Paired-species range | 315–339 |
| Median Spearman rho | 0.789103 |
| 5th percentile rho | 0.765165 |
| 95th percentile rho | 0.810940 |
| Median Spearman-Brown projected reliability | 0.882121 |
| Median Lin CCC | 0.854831 |
| Median absolute split difference in D | 0.078981 |

Frozen primary verdict:

`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`

Later strict deterministic stress test:

- reserve rho = 0.792727693;
- bootstrap 95% CI = 0.741867914–0.832365877;
- strict floor = 0.80;
- strict verdict = not supported.

The later stress test constrains the claim; it does not rewrite the chronologically earlier primary H1 result.

## S3. Legacy H2 target localization

Target freeze:

- `docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`

Canonical result:

- `results/polymorphism_white_axis_targeted_test_20260912/result.json`

Status:

`targeted_post_audit_decomposition_not_untouched_confirmatory_test`

Fixed target:

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`

Statistic:

`W = mean_i (u_i dot q_white)^2`

### Legacy targeted result

| Tier | Cohort | Species | Observed W | Null median | Upper-tail p |
|---|---|---:|---:|---:|---:|
| 0.10 | Discovery | 152 | 0.514625 | 0.430808 | 0.001 |
| 0.10 | Reserve | 129 | 0.514586 | 0.466546 | 0.001 |
| 0.20 | Discovery | 75 | 0.542355 | 0.443626 | 0.001 |
| 0.20 | Reserve | 65 | 0.510517 | 0.469943 | 0.008 |

Frozen legacy verdict:

`WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`

The named q_white target was isolated after the original broad H2 geometry had been opened. These rows are therefore target-localization evidence, not untouched prospective confirmation.

## S4. Third-cohort prospective H2 chain of custody

Selection protocol:

- `docs/POLYMORPHISM_H2_THIRD_COHORT_SELECTION_PROTOCOL_20260916.md`

Prospective measurement protocol:

- `docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md`

Result/claim freeze:

- `docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`

Machine-readable support result:

- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`

Machine-readable H2 result:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

### Outcome-blind selection

| Quantity | Value |
|---|---:|
| Candidate frame | 3,230 species |
| Frozen selected set | 500 species |
| Candidate CSV SHA256 | `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e` |
| Selected CSV SHA256 | `4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a` |
| Selected manifest SHA256 | `16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59` |
| Selection-freeze commit | `2768b2dd0f8baf4ed1185a1b64c1060768b36c00` |

### Pre-opening technical qualification

| Quantity | Value |
|---|---|
| Qualifying head | `7fcabcda61c90aaac3ad84e7f53ad5cd6c8c5663` |
| Qualification run | `35065814180` |
| Qualification job | `104695722267` |
| Tests | 17 / 17 PASS |
| Artifact ID | `10434700102` |
| Artifact digest | `sha256:2727d23a516cbf162b60606ba36932bd4fc1a4942ddc4242a2a03d18a6692e74` |

This stage contained synthetic-only technical evidence and no biological H2 outcome.

### Biological execution

| Quantity | Value |
|---|---|
| Execution head | `f7582767da237aeb344dfc072452c6daf92e662b` |
| Workflow run | `35177668182` |
| Terminal reassembly/test job | `105194354979` |
| Immutable biological result commit | `7e538e5c51c05a7cc47b2fcf53eea92634c8a863` |
| Final artifact ID | `10496492307` |
| Final artifact digest | `sha256:319c040aaffc30d3cd97dcbcc217409e75010ec0f67bcc86c6bb82d8275779c7` |

### Measurement-support gate

| Quantity | Value |
|---|---:|
| Authorized species | 499 |
| Target rows per species | 100 |
| Terminal rows | 49,900 |
| Unique measurement IDs | 49,900 |
| Duplicate measurement IDs | 0 |
| Terminal partitions | 256 / 256 |
| Classifiable rows | 25,788 |
| Nonclassifiable rows | 24,112 |
| Measurement-evaluable species | 377 |
| Required minimum | 250 |
| Replacement species | 0 |
| Replacement rows | 0 |
| Persisted image pixels | false |
| Support decision | PASS |

### Prospective H2 result

| Tier | Species | Observed W | Null median | Null 95% interval | p | Decision |
|---|---:|---:|---:|---|---:|---|
| Primary 0.10 | 158 | 0.5172457461 | 0.4571428150 | 0.4358491120–0.4752987776 | 0.001 | PASS |
| Strict 0.20 | 86 | 0.5329282123 | 0.4593196659 | 0.4328679570–0.4867224043 | 0.001 | PASS |

Frozen terminal verdict:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

Interpretation boundary:

This is prospective species-disjoint confirmation within the same iNaturalist opportunity universe, not independent-source replication.

## S5. Replicated D–spatial organization

Reporting-only machine-readable receipt:

- `results/polymorphism_spatial_organization_clue_20260918/result.json`

The receipt copies already frozen values from PR #32 / `feat/polymorphism-paper-v0-1-post-step9` and records the original Git blob identities. It performs no new biological analysis.

### Core result

| Cohort / response | Observed association | Geometry-preserving p |
|---|---:|---:|
| Discovery raw D–spatial | rho = 0.0892133 | 0.034 |
| Reserve raw D–spatial | rho = 0.1016008 | 0.025 |
| Discovery, span + clear technical-failure adjusted | partial rho = 0.1266367 | 0.007 |
| Reserve, span + clear technical-failure adjusted | partial rho = 0.0992877 | 0.025 |
| Reserve matched flower-minus-background | partial rho = 0.1162411 | 0.010 |

Reserve uniform ambiguity-endpoint stress tests remain supported:

- primary D_min4: rho = 0.0970781, p = 0.029;
- primary D_max4: rho = 0.1252858, p = 0.008;
- flower-minus-background D_min4: rho = 0.1162986, p = 0.009;
- flower-minus-background D_max4: rho = 0.1327852, p = 0.006.

Interpretation boundary:

Greater D is associated with stronger within-species geographic colour organization. This is a replicated structural correlate and mechanistic clue, not evidence that geographic organization causes D or that any particular climate, pollinator, demographic, gene-flow or selection mechanism has been identified.

## S6. P500 terminal status

Canonical postmortem:

- `docs/P500_PROSPECTIVE_TERMINAL_POSTMORTEM_20260916.md`

P500 completed prospective measurement transport:

- 499 species;
- 49,900 rows;
- 256/256 measurement partitions;
- 373 measurement-evaluable species;
- support gate PASS.

However, the H2 calculation failed during post-calculation serialization before a durable terminal H2 result was written.

Therefore P500 supplies neither prospective H2 confirmation nor prospective H2 refutation and is not replayed under the frozen one-shot rule.

## S7. H3a broad phylogenetic-signal boundary

Canonical manifest:

- `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`

Workflow provenance:

- run: `34677042793`
- artifact ID: `10292218669`
- artifact digest: `sha256:7f95699149f3111e00b9a095d78ddfa32e1111727d12ff1bb3ce9dc9c1dbd892`

Reserve has 341 retained tips in each S1-S3 placement scenario.

| Scenario | K | p(K) | lambda | p(lambda=0) | opportunity-adjusted p |
|---|---:|---:|---:|---:|---:|
| S1 | 0.0710190050 | 0.2716 | 0.0485656664 | 0.1656926113 | 0.3500 |
| S2 | 0.0601475598 | 0.4134 | 0.0511287826 | 0.1555141961 | 0.2464 |
| S3 | 0.0707576659 | 0.2674 | 0.0474641312 | 0.1702377597 | 0.3510 |

Frozen verdict:

`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`

This closes only the tested broad tree-wide signal claim.

## S8. H3b sampled-span replication boundary

Canonical result freeze:

- `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Workflow provenance:

- run: `34677468362`
- artifact ID: `10292399238`
- artifact digest: `sha256:34c5e646725f6865e313b17f1b70f2471db8169f443fb0649da83044d654386c`

| Cohort | n | rho(D, log1p sampled span) | permutation p |
|---|---:|---:|---:|
| Discovery calibration | 369 | 0.1798786 | 0.00089996 |
| Reserve replication | 363 | -0.0025855 | 0.9586021 |

Reserve adjusted partial-rank:

- rho = 0.0055187;
- p = 0.9162042.

S1-S3 reserve rank-PGLS:

- beta = -0.0055931;
- p = 0.9136754.

Frozen verdict:

`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`

Sampled photographic span is not true biological range size.

## S9. Canonical main-text figures

Generated reporting-only from frozen results:

- `docs/figures/polymorphism_20260918/polymorphism_figure1_measurement_frame.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure2_h1_reproducibility.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure3_h2_target_localization.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure4_prospective_h2.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure5_explanatory_boundaries.{png,pdf}`

Canonical figure manifest:

- `docs/figures/polymorphism_20260918/polymorphism_figure_manifest_20260918.json`

Figure-generation script:

- `scripts/analysis/make_polymorphism_manuscript_figures.py`

Figure test:

- `tests/test_make_polymorphism_manuscript_figures.py`

The figure manifest records PNG/PDF SHA256 values and `scientific_claims_changed = false`.

## S10. Manuscript claim guard

Claim-guard test:

- `tests/test_polymorphism_manuscript_claims.py`

Workflow:

- `.github/workflows/polymorphism-manuscript-claim-guard.yml`

Most recent verified post-architecture README run:

- workflow run: `35299258800`
- conclusion: success.

The guard checks the current manuscript, claim ledger, figure plan and README against the frozen third-cohort values and required claim boundaries.

## S11. Literature-positioning audit

Bounded manuscript literature audit:

- `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`

This literature layer supports context and interpretation only. It cannot modify the machine-readable empirical verdicts.

## S12. Hard nonclaims carried into all supplementary material

No main-text or supplementary output may claim:

- global prevalence of flower-colour polymorphism from the high-depth cohorts;
- independent-source replication of third-cohort H2;
- pigment chemistry or pigment-loss/gain mechanism;
- evolutionary direction of white/nonwhite transitions;
- pollinator, climate or other adaptive causation;
- a recurrent non-white hue axis;
- absence of all phylogenetic structure;
- irrelevance of true biological range size;
- near-perfect or split-invariant H1 reliability;
- a durable P500 biological H2 verdict.

If a future document conflicts with a committed machine-readable result, the machine-readable result controls.
