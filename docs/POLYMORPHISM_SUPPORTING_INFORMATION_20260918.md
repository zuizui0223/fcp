# Flower-colour polymorphism Supporting Information map — 2026-09-18

This document is the compact Supporting Information index for the active flower-colour polymorphism paper.

A direct claim-to-input routing table is maintained at `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`.

A self-contained verification bundle is maintained under release tag `fcp-np-provenance-20260926` (https://github.com/zuizui0223/fcp/releases/tag/fcp-np-provenance-20260926). The canonical Git-tracked receipt, `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`, records the source commit, asset byte count, SHA256 and packaged-file count for the currently published asset.

It does not recompute any result. It organizes frozen protocols, machine-readable outputs, figures and provenance into a manuscript-facing evidence map.

### Secondary PAL/WAL maintenance evidence

- Within-*Silene littorea* persistence receipt: `results/polymorphism_silene_decoupling_persistence_20260925/result.json`.
- Cross-species interval/censoring-aware receipt: `results/polymorphism_crossspecies_pal_wal_frequency_20260925/result.json`.
- Source: Del Valle et al. (2019), Supplementary Tables S2 and S1, respectively.

| Contrast | PAL | WAL |
|---|---:|---:|
| *S. littorea* positive-frequency median | 15.5% | 0.21% |
| *S. littorea* positive-frequency range | 8–21% | 0.05–0.86% |
| Cross-system instances | 13 | 13 |
| Cross-system median bound used for comparison | lower bound 5% | numeric upper bound 0.1% |
| Largest quantified WAL upper bound | — | 1.4% |

In *S. littorea*, every positive PAL frequency exceeded the maximum positive WAL frequency. Across the 13+13 literature systems, 7/13 PAL lower bounds exceeded the 1.4% maximum quantified WAL upper bound, and every PAL upper endpoint exceeded it. These analyses are descriptive reanalyses; literature ascertainment, heterogeneous survey effort and censored/qualitative frequency reporting preclude an unbiased meta-analytic or causal interpretation.

## S1. Data sources, selection and inferential necessity

Primary provenance documents:

- `docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md` — outcome-blind global opportunity-frame lineage;
- `docs/POLYMORPHISM_METHODS_CLASSIFICATION_20260918.md` — separation of established statistical components from study-specific sampling, null-model and prospective-confirmation design;
- `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md` — claim-to-input routing, immutable commits, hashes and archive routes.

### Global opportunity frame

Metadata-only iNaturalist discovery across an 18 × 9 equal-area grid produced **42,111 unique species** before candidate image pixels or flower-colour outcomes were opened. A capacity scan identified **4,730 species** with >=100 retainable photographs after observer capping.

**Necessity:** the candidate universe must be defined independently of flower colour so visually conspicuous or previously known polymorphic species are not preferentially admitted. The frame is a sampling opportunity universe, not a flower-colour dataset and not a polymorphism-prevalence denominator.

### Image resources and annotations

| Data resource | Source / selection | Scale | Why required | Role |
|---|---|---:|---|---|
| Discovery cohort | Fixed high-depth iNaturalist acquisition; colour-blind filters, observer cap, geographic maximin sampling | 500 species × 100 photos; 369 D-eligible | Estimate within-species distributions deeply enough for hypothesis discovery | Discovery/calibration |
| Species-disjoint validation cohort | Same acquisition contract, different 500 species | 500 species × 100 photos; 363 D-eligible | Test reliability and spatial associations on species not used for discovery | Validation within same source system |
| Fresh-image D transport | New iNaturalist photo IDs, same frozen measurement system | 136 species evaluable in both executions | Test transport beyond the original observer partitions to a fresh image sample | Same-source transport |
| Prospective confirmation cohort | Previously used high-depth species excluded; 500 species frozen outcome-blind from 3,230 candidates; no replacement | 499 species × 100 photos = 49,900; 377 measurement-evaluable | The white axis was data-derived in the first resource, so a new species/photo cohort was needed for prospective confirmation | Primary H2 confirmation |
| Highlight-validity reacquisition | Same 49,900 prospective rows reacquired only after H2 terminalization | 44,098 rows with highlight metrics | Quantify image-exposure coupling of the measured white state | Post-confirmatory validity bound |
| Climate annotation | WorldClim 2.1 BIO5, BIO14 and solar radiation joined to prospective-cohort coordinates; BIO5 transported to discovery/validation | 281 species, 12,583 rows in primary panel | Test a simple abiotic explanation and whether it generalizes | Secondary explanatory test |
| Phylogenetic annotation | V.PhyloMaker2 S1–S3 placements joined to validation D | 341 validation species | Test broad shared ancestry as an alternative explanation | Alternative-explanation filter |

The high-depth image acquisition contract required Research Grade species-rank iNaturalist records with photographs and georeferences, flowering annotation, positional accuracy <=5 km, unobscured/open coordinates and allowed CC licences. Observer contribution was capped at two photographs per species and deterministic geographic maximin sampling fixed 100 raw photographs per species. No native-range restriction or explicit captive/wild filter was imposed.

The central inferential sequence is therefore **outcome-blind frame → discovery/validation → target freeze → prospective confirmation**. Highlight, climate and phylogenetic data are annotations of those image resources, not additional independent flower-colour cohorts.


### Finite-sample sensitivity of Gini–Simpson D

Machine-readable reporting receipt:

- `results/polymorphism_D_finite_sample_sensitivity_20260928/result.json`

Raw D is the primary estimand. The sensitivity uses

[
D_{\mathrm{corr}} = 1-\sum_k\frac{n_k(n_k-1)}{n(n-1)}
= D\frac{n}{n-1}.
]

Under independent multinomial sampling this is the finite-sample-corrected Gini–Simpson estimator; it targets plug-in bias only and does not correct observational representativeness or image/classification error. With n >= 40, the largest possible multiplier is 1.0256; median n was 56 in both high-depth cohorts.

| Quantity | Discovery | Validation |
|---|---:|---:|
| Species | 369 | 363 |
| n_classifiable, median (min–max) | 56 (40–86) | 56 (40–92) |
| Spearman(raw D, corrected D) | 0.999977 | 0.999978 |
| Mean absolute D change | 0.00486 | 0.00469 |
| Maximum absolute D change | 0.01646 | 0.01740 |
| D–spatial adjusted partial rho, raw | 0.126637 | 0.099288 |
| D–spatial adjusted partial rho, corrected | 0.126256 | 0.099561 |
| Geometry-preserving p, raw / corrected | 0.007 / 0.007 | 0.025 / 0.025 |

In validation, the matched flower-minus-background result was likewise unchanged (raw partial rho = 0.116241, corrected = 0.116155; p = 0.010 for both). The sampled-span replication remained absent after correction (raw rho = -0.0026, p = 0.9586; corrected rho = -0.0024, p = 0.9629). Validation Blomberg-K effect sizes were also nearly unchanged: S1 0.071019 → 0.071224, S2 0.060148 → 0.060262 and S3 0.070758 → 0.070962.

Thus finite-sample correction does not alter the substantive D-based conclusions.

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

### Primary validation-cohort result

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

- validation rho = 0.792727693;
- bootstrap 95% CI = 0.741867914–0.832365877;
- strict floor = 0.80;
- strict verdict = not supported.

The later stress test constrains the claim; it does not rewrite the chronologically earlier primary H1 result.

## S3. Discovery/validation H2 target localization

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

### Original-cohort targeted result

| Tier | Cohort | Species | Observed W | Null median | Upper-tail p |
|---|---|---:|---:|---:|---:|
| 0.10 | Discovery | 152 | 0.514625 | 0.430808 | 0.001 |
| 0.10 | Validation | 129 | 0.514586 | 0.466546 | 0.001 |
| 0.20 | Discovery | 75 | 0.542355 | 0.443626 | 0.001 |
| 0.20 | Validation | 65 | 0.510517 | 0.469943 | 0.008 |

Frozen targeted-test verdict:

`WHITE_AXIS_TARGETED_SUPPORT_PRIMARY_AND_STRICT`

The named q_white target was isolated after the original broad H2 geometry had been opened. These rows are therefore target-localization evidence, not untouched prospective confirmation.

## S4. Prospective-confirmation H2 chain of custody

The exact protocol and machine-readable result identities are routed through the current data-lineage map and self-contained provenance snapshot. The reader-facing quantities needed to verify selection, technical qualification and the terminal H2 result are reported below.

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

This is prospective species-disjoint confirmation within the same iNaturalist opportunity universe, not independent-source replication. The confirmatory contrast is observed W versus the coarse-state-preserving structured null; it does not establish that the coarse white state itself is artifact-free.

## S4b. Post-confirmatory H2 validity and direct highlight audit

Exact protocol/result identifiers are preserved in the data-lineage map and provenance snapshot. The analyses below were performed after the frozen prospective-confirmation H2 verdict and do not replace its prospective decision rule.

### Inferential decomposition

| Quantity | Value |
|---|---:|
| Isotropic 8-D expectation | 0.125 |
| Frozen structured-null median | 0.457143 |
| Observed W | 0.517246 |
| Observed / null median | 1.1315 |
| Primary vector species | 158 |
| White among primary/secondary coarse morphs | 137 / 158 (86.7%) |
| Mean species contribution, white included | 0.5622 |
| Mean species contribution, white absent | 0.2239 |

Thus the frozen prospective test establishes excess alignment beyond a coarse-state-preserving baseline rather than the full observed alignment relative to isotropy.

### Gate-reapplied null

A 299-replicate post-confirmatory sensitivity re-applied the continuous minor-cluster gate in every null world, starting from the 185 coarse-gate species.

- median W = **0.447495**;
- mean W = **0.447583**;
- 95% interval = **0.428986–0.469580**;
- plus-one upper-tail p = **0.00333**;
- retained vector species = 166–183, median 176.

The median is lower than the frozen conditional-null median (0.457143), so the frozen null is slightly more conservative with respect to this gate asymmetry.

### Direct digital-highlight control

The one-shot direct control reacquired all frozen prospective-cohort image rows before joining biological outcomes.

| Quantity | Value |
|---|---:|
| Frozen reacquisition rows | 49,900 |
| Acquisition failures | 0 |
| Source-byte drift | 0 |
| Highlight metrics available | 44,098 |
| High-clip threshold | near_clip_fraction > 0.4933450501 |
| High-clip rows | 2,205 |
| Coupling-model rows | 23,320 |
| Coupling-model species | 461 |
| OR per 1 within-species SD near-clip | 1.444493 |
| 95% CI | 1.389332–1.501845 |

The entire confidence interval is above the predeclared negligible-coupling upper bound of 1.25. This is direct evidence that white classification is coupled to digital highlight exposure within species.

### High-clip exclusion sensitivity

After removing the response-blind high-clip set:

- vector species = **142 / 158 (89.9%)**;
- W = **0.503428**;
- structured-null median = **0.453801**;
- structured-null p = **0.001**;
- H2 support remains **positive**.

The frozen executable required >=90% retention (>=143 vectors) for the sensitivity gate to clear. Because 142 vectors remained, its generated terminal state was `INDETERMINATE`.

The separately frozen prose protocol contains a FLAGGED coupling clause because the complete OR interval exceeds 1.25. Since the prose clause and executable precedence were both frozen before outcome opening, the machine state is not retrospectively recoded. The empirical interpretation is reported directly: exposure coupling is present, while H2 excess alignment survives high-clip exclusion.

### disttrait audit

The later generic `disttrait` implementation is not bitwise identical to the frozen study-specific H2 code on the full prospective-cohort data:

- frozen W = **0.5172457461**;
- `disttrait.two_mode_axis` W = **0.5183899565**;
- `disttrait.structured_alignment_null` observed route = **0.5157982982**.

Near-tied deterministic initializations can switch after defensive row renormalization. The paper's numerical results therefore remain controlled by the frozen study-specific pipeline.

## S4c. Post-confirmatory environmental filter and BIO5 transport

Machine-readable environmental and transport result identities are routed through the current data-lineage map and provenance snapshot.

This analysis was specified after the prospective H2 result had been terminalized. It is therefore a secondary mechanistic/ecological analysis and cannot alter the frozen H2 verdict.

### Prospective-cohort environmental filter

After response-blind high-clip exclusion:

| Variable | Eligible species | Median white-minus-nonwhite contrast | Holm-adjusted species-level p | Conditional OR per within-species SD | Conditional p | Frozen gate |
|---|---:|---:|---:|---:|---:|---|
| BIO5 | 281 | +0.0690 SD | 0.0354 | 1.073 | 0.000919 | PASS |
| BIO14 | 281 | -0.0187 SD | 0.692 | 0.993 | 0.749 | FAIL |
| Mean solar radiation | 281 | -0.00683 SD | 0.692 | 0.987 | 0.544 | FAIL |

BIO5 was therefore the only prespecified environmental variable to pass the prospective-cohort mechanism gate.

### Species-disjoint BIO5 transport

| Cohort | Eligible species | Median white-minus-nonwhite BIO5 contrast | Species-level p | Conditional OR | Conditional p | Primary support |
|---|---:|---:|---:|---:|---:|---|
| Discovery | 271 | -0.0068 SD | 0.743 | 1.033 | 0.131 | false |
| Validation | 260 | +0.0662 SD | 0.0541 | 1.048 | 0.0319 | false |

Frozen transport decision:

**Not supported under the prespecified cross-cohort rule.**

Interpretation boundary:

The prospective confirmation cohort contains a prospectively specified within-cohort association between white states and warmer BIO5 environments, but the effect does not transport as a common rule across the original species-disjoint cohorts. This is consistent with context-dependent environmental sorting and does not establish causal heat selection, a universal temperature effect, or independence from the known exposure coupling of the white classifier.

WorldClim provenance: exact BIO/SRAD input archives are mirrored under release tag `fcp-worldclim-2.1-10m-20260925`; checksum authority is `archive/fcp_submission_20260925/worldclim_checksums.txt`.

## S5. Replicated D–spatial organization

Reporting-only machine-readable receipt:

- `results/polymorphism_spatial_organization_clue_20260918/result.json`

The reporting receipt copies previously frozen values and performs no new biological analysis.

### Core result

| Cohort / response | Observed association | Geometry-preserving p |
|---|---:|---:|
| Discovery raw D–spatial | rho = 0.0892133 | 0.034 |
| Validation raw D–spatial | rho = 0.1016008 | 0.025 |
| Discovery, span + clear technical-failure adjusted | partial rho = 0.1266367 | 0.007 |
| Validation, span + clear technical-failure adjusted | partial rho = 0.0992877 | 0.025 |
| Validation matched flower-minus-background | partial rho = 0.1162411 | 0.010 |

Validation-cohort uniform ambiguity-endpoint stress tests remain supported:

- primary D_min4: rho = 0.0970781, p = 0.029;
- primary D_max4: rho = 0.1252858, p = 0.008;
- flower-minus-background D_min4: rho = 0.1162986, p = 0.009;
- flower-minus-background D_max4: rho = 0.1327852, p = 0.006.

Interpretation boundary:

Greater D is associated with stronger within-species geographic colour organization. This is a replicated structural correlate and mechanistic clue, not evidence that geographic organization causes D or that any particular climate, pollinator, demographic, gene-flow or selection mechanism has been identified.

## S6. H3a broad phylogenetic-signal boundary

Canonical manifest:

- `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`

Workflow provenance:

- run: `34677042793`
- artifact ID: `10292218669`
- artifact digest: `sha256:7f95699149f3111e00b9a095d78ddfa32e1111727d12ff1bb3ce9dc9c1dbd892`

The validation cohort has 341 retained tips in each S1-S3 placement scenario.

| Scenario | K | p(K) | lambda | p(lambda=0) | opportunity-adjusted p |
|---|---:|---:|---:|---:|---:|
| S1 | 0.0710190050 | 0.2716 | 0.0485656664 | 0.1656926113 | 0.3500 |
| S2 | 0.0601475598 | 0.4134 | 0.0511287826 | 0.1555141961 | 0.2464 |
| S3 | 0.0707576659 | 0.2674 | 0.0474641312 | 0.1702377597 | 0.3510 |

Frozen verdict:

`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`

This closes only the tested broad tree-wide signal claim.

## S7. H3b sampled-span replication boundary

The exact result/protocol identity is routed in the data-lineage map. Workflow provenance:

- run: `34677468362`
- artifact ID: `10292399238`
- artifact digest: `sha256:34c5e646725f6865e313b17f1b70f2471db8169f443fb0649da83044d654386c`

| Cohort | n | rho(D, log1p sampled span) | permutation p |
|---|---:|---:|---:|
| Discovery calibration | 369 | 0.1798786 | 0.00089996 |
| Validation replication | 363 | -0.0025855 | 0.9586021 |

Validation adjusted partial-rank:

- rho = 0.0055187;
- p = 0.9162042.

S1-S3 validation rank-PGLS:

- beta = -0.0055931;
- p = 0.9136754.

Frozen verdict:

`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`

Sampled photographic span is not true biological range size.

## S8. Canonical main-text figures

Generated reporting-only from frozen results:

- `docs/figures/polymorphism_20260918/polymorphism_figure1_measurement_frame.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure2_h1_reproducibility.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure3_h2_target_localization.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure4_prospective_h2.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figure5_explanatory_boundaries.{png,pdf}`
- `docs/figures/polymorphism_20260918/polymorphism_figureS9_secondary_mechanism_evidence.{png,pdf}`

Canonical figure manifest:

- `docs/figures/polymorphism_20260918/polymorphism_figure_manifest_20260918.json`

Figure-generation script:

- `scripts/analysis/make_polymorphism_manuscript_figures.py`

Figure test:

- `tests/test_make_polymorphism_manuscript_figures.py`

The figure manifest records PNG/PDF SHA256 values and `scientific_claims_changed = false`.

## S9. Manuscript claim guard

Claim-guard test:

- `tests/test_polymorphism_manuscript_claims.py`

Workflow:

- `.github/workflows/polymorphism-manuscript-claim-guard.yml`

Most recent verified post-architecture README run:

- workflow run: `35299258800`
- conclusion: success.

The guard checks the current manuscript, claim ledger, figure plan and README against the frozen prospective-confirmation values and required claim boundaries.

## S10. Literature-positioning audit

Bounded manuscript literature audit:

- `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`

This literature layer supports context and interpretation only. It cannot modify the machine-readable empirical verdicts.

## S11. Hard nonclaims carried into all supplementary material

No main-text or supplementary output may claim:

- global prevalence of flower-colour polymorphism from the high-depth cohorts;
- independent-source replication of prospective H2;
- pigment chemistry or pigment-loss/gain mechanism;
- evolutionary direction of white/nonwhite transitions;
- pollinator, climate or other adaptive causation;
- a universal or replicated BIO5–white association across cohorts;
- a recurrent non-white hue axis;
- absence of all phylogenetic structure;
- irrelevance of true biological range size;
- near-perfect or split-invariant H1 reliability;
- that the coarse white state is free of digital exposure, background-context or ROI-contamination effects;
- numerical identity between the frozen biological H2 implementation and the later generic disttrait package.

If a future document conflicts with a committed machine-readable result, the machine-readable result controls.
