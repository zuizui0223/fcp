# Flower-colour polymorphism — integrated manuscript evidence ledger

Date frozen: 2026-09-13 JST

## Scope

This ledger defines the current manuscript-level evidence hierarchy after closing the H1 measurement-validity, H2 colour-space geometry, H3a phylogenetic-signal, and H3b sampled-span branches.

It is a synthesis of already frozen analyses. It does not create a new biological test and cannot be used to rescue a failed branch by changing thresholds, endpoints, or cohort roles.

## 1. Sampling frame and cohort roles

### Global frame

- global species frame: **42,111 species**
- this frame is the universe for opportunity/accounting and future expansion
- it is **not** valid to infer global polymorphism prevalence from the existing high-depth validation cohorts

### High-depth validation cohorts

- discovery: **369 species** with frozen four-state `n_classifiable >= 40`
- reserve: **363 species** with frozen four-state `n_classifiable >= 40`
- discovery and reserve are species-disjoint

Their role is measurement/geometry validation under deep photographic sampling. They are not a random or representative 732-species sample from the 42,111-species frame.

### U100 expansion frame

Previously audited high-depth opportunity at >=100 raw photographs after the frozen observer cap:

- **4,730 species**

This is the natural prospective expansion frame. It is an opportunity set, not an already observed H1/H2 analysis cohort.

## 2. H1 — can continuous species-level polymorphism be measured reproducibly?

### Frozen quantity

For the four biological coarse states

- white
- yellow_orange
- red_pink
- blue_purple

using only `global_classifiable` rows,

`D = 1 - sum_k p_k^2`.

`mixed_uncertain` is never promoted to a biological state.

### Canonical observer-disjoint reliability evidence

Evidence hierarchy is fixed in:

`docs/POLYMORPHISM_H1_RELIABILITY_EVIDENCE_LEDGER_20260913.md`

The earliest outcome-blind observer-disjoint protocol is the canonical primary analysis:

`docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md`

Protocol-freeze commit:
`92488ada5f4f535f6731bba6519101d8ed8a7e1a`

Canonical run:
`34707537360`

Canonical reserve result across 200 outcome-blind observer partitions:

- defined rho partitions: **200/200**
- median paired species: **329**
- median Spearman rho(D_A,D_B): **0.789103**
- rho 5th percentile: **0.765165**
- rho 95th percentile: **0.810940**
- median diagnostic Lin CCC: **0.854831**
- median Spearman-Brown projected full-estimate reliability: **0.882121**
- median |D_A-D_B|: **0.078981**
- median signed bias A-B: **-0.000980**

All three frozen reserve decision conditions pass.

Canonical verdict:

**`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED+DISCOVERY_CONSISTENT`**

A later single deterministic observer split produced reserve rho **0.792728**, bootstrap 95% CI **[0.741868, 0.832366]**, and CCC **0.847429**. Its later engineering cutoff of rho >=0.80 was narrowly missed; by freeze chronology this is a sensitivity analysis and cannot replace the earlier primary rule. Numerically it falls inside the canonical 200-partition distribution.

The historical approximate reliability value ~0.971 is not canonical because repository audit did not establish that it referred to this exact current four-state D endpoint.

### Allowed H1 claim

Within the existing high-depth validation design, continuous four-state flower-colour polymorphism D is reproducibly recoverable from observation sets contributed by disjoint observers.

### H1 nonclaims

H1 does not establish global prevalence, cohort representativeness, perfect image-level labels, evolutionary cause, adaptation, or environmental drivers.

## 3. H2 — along what colour-space direction does polymorphism recur?

Canonical synthesis:

`docs/POLYMORPHISM_H2_EVIDENCE_LEDGER_20260912.md`

### 3.1 Label-free geometry

Within H1-admitted species, H2 reconstructs two modes without using the four coarse H1 labels, from the normalized nine-colour palette:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze`.

The mode displacement vector is `Delta_i`; the sign-invariant unit axis is `u_i = Delta_i / ||Delta_i||`.

Primary 10% gate:

- discovery H2 N = **152**, lambda1 = **0.54143**
- reserve H2 N = **129**, lambda1 = **0.53496**
- reserve projection on frozen discovery axis = **0.52418**
- independent discovery/reserve leading-axis dot = **0.98571**

Strict 20% gate:

- discovery H2 N = **75**, lambda1 = **0.57750**
- reserve H2 N = **65**, lambda1 = **0.56185**
- reserve projection = **0.54826**
- axis dot = **0.98458**

### 3.2 Coarse-state-preserving structured null

Workflow run:
`34674420093`

Verdict:

**`H2_SURVIVES_COARSE_STATE_PRESERVING_STRUCTURED_NULL`**

Primary:

- discovery observed lambda1 **0.54143** vs null median **0.47622**, p = **0.001**
- reserve transport **0.52418** vs null median **0.48269**, p = **0.001**

Strict:

- discovery **0.57750** vs null median **0.48359**, p = **0.001**
- reserve transport **0.54826** vs null median **0.47660**, p = **0.001**

Thus broad geometry cannot be explained solely by the frozen coarse-state composition and global coarse-state-to-palette mapping.

### 3.3 White-axis audit

Most H2 species have a top-two coarse pair containing white:

- primary discovery: **124/152 = 81.58%**
- primary reserve: **110/129 = 85.27%**
- strict discovery: **65/75 = 86.67%**
- strict reserve: **54/65 = 83.08%**

After projecting out the canonical white-versus-nonwhite contrast, residual geometry does **not** survive the same structured control:

- primary discovery residual structured p = **1.0**
- primary reserve transport structured p = **1.0**
- strict discovery structured p = **1.0**
- strict reserve transport structured p = **0.962**

Verdict:

**`H2_DOMINATED_BY_WHITE_VERSUS_NONWHITE_AXIS`**

There is no supported recurrent non-white hue direction after construction control.

### 3.4 Frozen targeted white-axis decomposition

Canonical fixed contrast:

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`

Target statistic:

`W = mean_i (u_i^T q_white)^2`.

The structured null is the same coarse-state-preserving null.

Freeze:
`docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`

Run:
`34675697583`

Verdict:

**`WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`**

Primary:

- discovery N=152, W=**0.514625**, null median=**0.430808**, p=**0.001**
- reserve N=129, W=**0.514586**, null median=**0.466546**, p=**0.001**

Strict:

- discovery N=75, W=**0.542355**, null median=**0.443626**, p=**0.001**
- reserve N=65, W=**0.510517**, null median=**0.469943**, p=**0.008**

### H2 evidence status

The white axis was isolated after the initial broad H2 geometry was opened. Therefore the current high-depth evidence is a strong **retrospective targeted decomposition with species-disjoint transport**, not an untouched prospective confirmation of the specific white axis.

### Allowed H2 claim now

Within the existing high-depth cohorts, flower-colour polymorphism is disproportionately aligned with an achromatic-chromatic **white-versus-nonwhite** axis beyond expectation from the frozen coarse-state composition and global state-to-palette mapping.

### H2 nonclaims

Do not claim:

- a general recurrent non-white hue axis;
- a specific pigment or genetic pathway;
- that white is ancestral or derived;
- adaptive loss/gain of pigmentation;
- pollinator, climate, or other ecological mechanism;
- global prevalence from the 369/363 cohorts.

## 4. H3a — broad phylogenetic structure of D

Frozen protocol:
`docs/POLYMORPHISM_H3A_PHYLOGENETIC_SIGNAL_PROTOCOL_20260912.md`

The V.PhyloMaker2 preflight passed before D was opened:

- discovery tree coverage: **368/369 = 99.73%** under each S1-S3 scenario
- reserve tree coverage: **341/363 = 93.94%** under each S1-S3 scenario

The fresh reserve primary Blomberg-K permutation tests are unsupported under all three placement scenarios:

- S1 p = **0.272**
- S2 p = **0.413**
- S3 p = **0.267**

Pagel-lambda and opportunity-adjusted sensitivities are also unsupported.

Verdict:

**`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`**

The previously observed discovery genus clustering cannot rescue reserve failure.

### Allowed H3a interpretation

The existing high-depth reserve cohort does not support broad phylogenetic conservatism of continuous D at the scale resolved by the frozen S1-S3 trees.

This is not evidence that all evolutionary history is irrelevant; it is a bounded negative result for this preregistered phylogenetic-signal endpoint.

## 5. H3b — sampled geographic span as a predictor of D

Frozen result:
`docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Discovery calibration:

- rho = **0.179879**
- permutation p = **0.000900**

Fresh species-disjoint reserve:

- rho = **-0.002586**
- permutation p = **0.958602**
- unbiased-D sensitivity p = **0.962902**
- observer/classifiable-adjusted rho = **0.005519**, p = **0.916204**
- S1-S3 rank-PGLS p approximately **0.914** in every scenario

Verdict:

**`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`**

The discovery association does not generalize and is not a manuscript-level predictor claim.

Trait-source coverage gates for pollination/life form/self-incompatibility were not relaxed post hoc; failed source coverage is not interpreted as a biological negative.

## 6. Current manuscript claim hierarchy

### Main empirical claim 1 — supported measurement property

**FCP can be quantified as a reproducible continuous species-level attribute under high-depth image sampling.**

Evidence: H1 observer-disjoint D reliability.

### Main empirical claim 2 — supported geometry, prospective target identified

**Observed FCP is preferentially organized along a white-versus-nonwhite colour-space axis in the existing high-depth cohorts.**

Evidence: label-free geometry, structured null, white-axis audit, targeted white-axis decomposition, reserve transport.

Status qualifier: the specific axis is retrospective because it was isolated after opening broad H2 geometry.

### Boundary result 1 — not supported

**Broad phylogenetic conservatism of D is not supported in the fresh reserve test.**

### Boundary result 2 — not supported

**The discovery sampled-span association is not replicated.**

## 7. What the paper is NOT

The present paper is not evidence for:

1. a global estimate of the fraction of flowering-plant species that are colour polymorphic;
2. a universal global flower-colour boundary;
3. a general non-white hue axis;
4. phylogenetic conservatism of D;
5. sampled geographic span as a general predictor of D;
6. pollinator, climate, pigment-pathway, genetic, or adaptive causation;
7. ancestral/derived direction of white versus pigmented morphs.

## 8. Single next promotion gate

The only analysis that can materially upgrade the positive biological conclusion is a **prospective test of the already frozen white-axis target in a fresh U100 expansion**.

The prospective test must:

- start from the previously defined U100 opportunity frame rather than hand-picked species;
- keep discovery and reserve species outside the fresh confirmatory cohort;
- retain the frozen four-state H1 admission system and nine-colour H2 palette;
- retain the fixed `q_white` axis and `W` statistic;
- retain the same coarse-state-preserving structured null;
- make the 10% gate the single primary endpoint and keep 20% as robustness/sensitivity;
- prohibit residual hue-axis searching if white-axis confirmation fails;
- report the full opportunity -> measurement -> H1 -> H2 eligibility funnel with reason-coded attrition;
- treat insufficient coverage as `UNRESOLVED_BY_COVERAGE`, not biological failure;
- preserve the hard nonclaims above.

Until such a prospective fresh-species test is completed, the specific white-axis claim remains retrospective/targeted despite its strong structured-null and reserve-transport evidence.

## 9. No-rescue rules

Do not:

- change H1 reliability thresholds after observing current outcomes;
- choose the later single-split H1 cutoff over the earlier frozen 200-partition primary based on which label is favorable;
- reinstate historical ~0.971 as canonical without exact endpoint provenance;
- hunt new H3 predictors after H3a/H3b closure;
- reinterpret trait source failures as biological zeros;
- replace the structured H2 null with an easier isotropic null;
- redefine white/nonwhite after the current targeted result;
- use discovery/reserve species as the prospective confirmation cohort;
- report U100 opportunity counts as observed polymorphism prevalence.
