# Flower-colour polymorphism mainline evidence ledger

Frozen: 2026-09-13 JST

Purpose: define the manuscript-facing evidence chain after closing H1, H2, H3a, and H3b. This ledger supersedes informal references to earlier exploratory polymorphism summaries. It does not replace the underlying protocol/result freezes.

## Executive status

The current paper has two positive claims and two bounded negative results.

1. **H1 positive — D is reproducibly measurable in the high-depth design.** The continuous four-state polymorphism score `D = 1 - sum p_k^2` is recovered consistently from observation sets contributed by completely disjoint observers.
2. **H2 positive but narrow/targeted — within-species colour variation is preferentially aligned with a white-versus-non-white axis.** The signal survives a coarse-state-preserving structured null, but the specific white axis was isolated after opening the broader geometry and therefore is retrospective/targeted in the present cohorts.
3. **H3a negative — no replicating broad phylogenetic signal in D.** Fresh reserve tests fail across S1-S3 phylogenetic placement scenarios, including opportunity-adjusted sensitivity.
4. **H3b negative — the discovery sampled-span association does not replicate.** Fresh reserve rho is essentially zero and rank-PGLS sensitivities agree.

The manuscript must therefore be framed as **measurement validity + constrained geometry**, not as a paper explaining which lineages or geographic-range classes evolve polymorphism.

---

## Sampling frame and cohort roles

### Global frame

- 42,111 species form the broad flower-colour sampling frame.
- This frame is useful for baseline colour composition and for defining the population from which high-depth opportunity can be audited.
- It is not currently a 42,111-species polymorphism prevalence estimate because repeated, classifiable, observer-supported information is uneven across species.

### High-depth validation cohorts

- discovery full-D cohort: 369 species;
- reserve full-D cohort: 363 species;
- discovery and reserve are species-disjoint;
- both are high-information validation cohorts, not an unbiased random sample of the 42,111-species frame.

The old 369-species discovery set must not be described as "the global sample". The 42,111-species frame is the population-level sampling universe; hypothesis-specific eligibility determines which species enter each validation analysis.

---

## H1 — direct observer-disjoint reliability of D

Canonical protocol:
`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`

Canonical result:
`results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`

Frozen result note:
`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_RESULT_FREEZE_20260913.md`

Canonical optimized workflow run: `34707537360`

Artifact: `10302466831`

Artifact digest: `sha256:19185ad9160b24ba47ef8d99a33d952ab5f5b6b62243ff756d41917bfc587c53`

### Design

- four biological states only: white, yellow_orange, red_pink, blue_purple;
- `mixed_uncertain` excluded from biological morph counts;
- observer-known classifiable n >= 40 and at least two observers;
- 200 deterministic outcome-blind observer partitions;
- observers never split between halves;
- primary requires >=20 classifiable observations per half;
- reserve is the decision cohort;
- frozen support thresholds: median paired N >=100, median split-half Spearman rho >=2/3, 5th percentile rho >=0.50.

### Reserve primary result

- observer-known eligible species: **363/363**;
- median paired N: **329** (315-339);
- median split-half Spearman rho: **0.789103**;
- rho 5th percentile: **0.765165**;
- rho 95th percentile: **0.810940**;
- median Spearman-Brown projected full-estimate reliability: **0.882121**;
- 5th percentile projected reliability: **0.866961**;
- median CCC: **0.854831**;
- median absolute split difference in D: **0.078981**;
- median signed A-B bias: **-0.000980**.

Frozen verdict:
**`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`**

Discovery is independently consistent:
- median paired N = 340;
- median rho = 0.811545;
- rho 5th percentile = 0.788380;
- projected full-estimate reliability = 0.895970.

Final label:
**`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED+DISCOVERY_CONSISTENT`**

### Implementation-equivalence audit

The original row-level implementation run `34707151264` subsequently completed successfully.

Artifact: `10302481007`

Artifact digest: `sha256:9255740cfe9218f433ebcbe52901a80cbc286d6c26c7e6d1db7a2cd0cb66f36d`

The original row-level `result.json` and optimized observer-preaggregation `result.json` are identical for every statistical value and decision field. The only JSON difference is the optimized implementation's explanatory `implementation_note`. Thus the runtime optimization did not change the frozen estimand or result.

### Allowed H1 claim

> Within the existing high-depth design, the four-state species-level polymorphism score D is reproducibly recoverable from completely observer-disjoint observation sets, supporting D as a measurable species attribute rather than an artifact of particular observers.

### H1 non-claims

H1 does not establish global polymorphism prevalence, global representativeness, perfect image-level classification, absence of real geographic morph-frequency structure, or a biological cause of D.

---

## H2 — recurrent geometry of within-species colour variation

Canonical evidence ledger:
`docs/POLYMORPHISM_H2_EVIDENCE_LEDGER_20260912.md`

### Broad label-free geometry

Primary 10% continuous-mode gate:
- discovery N = 152, lambda1 = **0.54143**, isotropic p approximately 1e-4;
- reserve N = 129, lambda1 = **0.53496**, isotropic p approximately 1e-4;
- reserve projection on frozen discovery axis = **0.52418**, p approximately 1e-4;
- independently fitted discovery/reserve axis absolute dot = **0.98571**.

Strict 20% gate:
- discovery N = 75, lambda1 = **0.57750**;
- reserve N = 65, lambda1 = **0.56185**;
- reserve projection = **0.54826**;
- independent-axis absolute dot = **0.98458**.

### Coarse-state-preserving structured null

Run: `34674420093`
Artifact: `10291618984`

The null preserves H2-selected species, coarse-state composition, and the global coarse-state-to-nine-colour palette mapping while destroying species-specific continuous palette structure.

Primary:
- discovery observed lambda1 0.54143 vs null median 0.47622, **p = 0.001**;
- reserve transport 0.52418 vs null median 0.48269, **p = 0.001**.

Strict:
- discovery 0.57750 vs null median 0.48359, **p = 0.001**;
- reserve transport 0.54826 vs null median 0.47660, **p = 0.001**.

Verdict:
**`H2_SURVIVES_COARSE_STATE_PRESERVING_STRUCTURED_NULL`**

### White-axis localization

Most H2 species have a top-two coarse pair involving white:
- primary discovery 124/152 = 81.58%;
- primary reserve 110/129 = 85.27%;
- strict discovery 65/75 = 86.67%;
- strict reserve 54/65 = 83.08%.

White-axis audit run: `34675375488`, artifact `10292067445`.

After projecting the canonical white-versus-equal-nonwhite contrast out:
- primary discovery residual structured-null p = **1.0**;
- primary reserve residual transport structured-null p = **1.0**;
- strict discovery residual structured-null p = **1.0**;
- strict reserve residual transport structured-null p = **0.962**.

The non-white-only discovery subsets likewise fail to establish a recurrent hue axis under the structured null.

Verdict:
**`H2_DOMINATED_BY_WHITE_VERSUS_NONWHITE_AXIS`**

### Targeted fixed white-axis statistic

Freeze:
`docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`

Canonical result:
`results/polymorphism_white_axis_targeted_test_20260912/result.json`

Run: `34675697583`
Artifact: `10292077657`

The fixed statistic is `W = mean_i (u_i^T q_white)^2`, with the same coarse-state-preserving structured null.

Primary:
- discovery W = **0.514625**, null median 0.430808, **p = 0.001**;
- reserve W = **0.514586**, null median 0.466546, **p = 0.001**.

Strict:
- discovery W = **0.542355**, null median 0.443626, **p = 0.001**;
- reserve W = **0.510517**, null median 0.469943, **p = 0.008**.

Verdict:
**`WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`**

### Allowed H2 claim

> In the existing high-depth validation cohorts, within-species flower-colour polymorphism is disproportionately aligned with an achromatic-chromatic, white-versus-non-white axis, beyond expectation from the frozen coarse-state composition and global coarse-state-to-palette mapping.

### Mandatory qualification

The specific white-axis target was isolated after the broader H2 geometry was opened. Therefore the current evidence is **retrospective/targeted**, not untouched confirmatory evidence for that specific axis. The axis, statistic, and structured null are now frozen for prospective expansion.

### H2 non-claims

Do not claim a recurrent non-white hue axis, a particular pigment pathway, white as ancestral or derived, adaptive causation, pollinator causation, climate causation, or global prevalence.

---

## H3a — broad phylogenetic signal

The H3a tree/preflight and signal analysis were frozen before opening reserve D. V.PhyloMaker2 used fixed upstream SHA `7af3fb5152f691af2e4ec9d5e2e467d1b50505e9` and S1-S3 placement scenarios.

Coverage gate passed:
- discovery: 368/369 tips = 99.73% for each S1-S3;
- reserve: 341/363 tips = 93.94% for each S1-S3.

Fresh reserve raw Blomberg K permutation p values:
- S1: **0.272**;
- S2: **0.413**;
- S3: **0.267**.

Pagel lambda is unsupported in all scenarios, and opportunity-adjusted K remains unsupported across S1-S3.

Frozen verdict:
**`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`**

### Allowed H3a statement

> The present high-depth reserve cohort provides no evidence for a broad phylogenetic signal in D under the frozen S1-S3 tree-placement scenarios.

Do not rewrite this as evidence that polymorphism has no phylogenetic history at any scale. Discovery genus clustering is calibration history and cannot rescue the fresh reserve failure.

---

## H3b — sampled geographic span

Frozen result:
`docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Run: `34677468362`
Artifact: `10292399238`
Artifact digest: `sha256:34c5e646725f6865e313b17f1b70f2471db8169f443fb0649da83044d654386c`

Discovery calibration:
- n = 369;
- rho = **0.179879**;
- 20,000-permutation p = **0.000900**;
- observer/classifiable-adjusted rho = 0.158509, p = 0.002700.

Fresh reserve:
- n = 363;
- rho = **-0.002586**;
- p = **0.958602**;
- unbiased-D rho = -0.002374, p = 0.962902;
- adjusted rho = 0.005519, p = 0.916204.

Frozen S1-S3 rank-PGLS sensitivity:
- beta_span_rank approximately -0.005593;
- p approximately **0.913675** in every placement scenario;
- lambda approximately 1e-7.

Verdict:
**`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`**

### Allowed H3b statement

> The discovery association between sampled geographic span and D does not generalize to the species-disjoint reserve cohort.

Do not equate sampled photographic span with true biological range size, and do not rescue this result with alternative span definitions, thresholds, or post-hoc predictor hunting.

---

## Integrated claim hierarchy

### Claim 1 — measurement foundation

**Supported and reserve-confirmed.** Flower-colour polymorphism can be represented as a continuous species-level four-state diversity score that is reproducible across completely disjoint observers in the high-depth design.

### Claim 2 — geometry

**Supported in both cohorts but specific-axis status is targeted/retrospective.** Variation is strongly concentrated along white versus non-white, and that concentration is not explained by the frozen coarse-state construction null. No residual recurrent non-white hue direction is supported.

### Claim 3 — broad predictors

**Not supported in fresh reserve tests.** Neither broad phylogenetic signal nor sampled photographic span provides a replicating explanation of D in the current validation design.

This asymmetry is scientifically useful: **the phenotype is measurable and geometrically constrained even though the tested broad lineage/geographic predictors do not generalize.**

---

## Manuscript hard boundaries

The manuscript must not claim:

- an unbiased global prevalence of flower-colour polymorphism;
- that the 369/363 species are globally representative;
- a globally recurrent non-white hue axis;
- that white is ancestrally lost or derived repeatedly;
- a specific pigment mechanism;
- adaptation, pollinator mediation, climate causation, or urban effects;
- absence of phylogenetic effects at finer taxonomic scales;
- absence of biological range-size effects;
- that H2 white-axis targeting was prospective in the current cohorts.

The failed global shared-boundary/RGFCA route remains an identifiability result and must not be blended into the positive FCP evidence chain.

---

## What is already closed versus what remains

### Closed

- D definition and high-depth eligibility;
- direct observer-disjoint H1 reliability;
- broad H2 label-free geometry;
- coarse-state-preserving H2 structured null;
- localization of H2 to white versus non-white;
- no-white and non-white-only diagnostics;
- current-cohort targeted white-axis test;
- H3a reserve phylogenetic signal test;
- H3b sampled-span reserve replication.

### Highest-value remaining evidence

The largest remaining inferential upgrade is **not another post-hoc predictor**. It is a prospective expansion using the now-frozen white-axis H2 target in a larger high-depth subset of the 42,111-species frame (e.g. U100/eligibility funnel) without changing `q_white`, W, the mode construction, or the structured null.

A successful prospective replication would convert the current strongest biological pattern from "targeted after audit, replicated across existing cohorts" to a genuinely prospective confirmatory result.