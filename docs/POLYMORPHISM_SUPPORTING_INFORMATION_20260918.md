# Flower-colour polymorphism Supporting Information map — 2026-09-18

This document is the compact Supporting Information index for the active flower-colour polymorphism paper.

A direct claim-to-input routing table is maintained at `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`.

A self-contained verification bundle is maintained under release tag `fcp-np-provenance-20260926` (https://github.com/zuizui0223/fcp/releases/tag/fcp-np-provenance-20260926). The canonical Git-tracked receipt, `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`, records the source commit, asset byte count, SHA256 and packaged-file count for the currently published asset.

It does not recompute any result. It organizes frozen protocols, machine-readable outputs, figures and provenance into a manuscript-facing evidence map.

## How to read this Supporting Information

The reader-facing order mirrors the paper's inferential chronology:

**sampling frame -> measurement validity -> geometry discovery/prospective confirmation -> geographic allocation of ITV -> spatial-process decomposition -> alternative explanations -> bounded mechanism evidence.**

Workflow IDs, immutable hashes and release receipts are retained where they establish provenance, but they are not additional biological results. The self-contained release bundle remains the deeper reproducibility layer.

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

The central data chronology remains **outcome-blind frame → discovery/validation → target freeze → prospective confirmation**. The later distributed-polymorphism and IBD/IBE-like analyses are explicitly post-outcome analyses applied across these frozen cohorts; highlight, climate and phylogenetic data are annotations rather than additional independent flower-colour cohorts.


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

### Post hoc observer sensitivities

These analyses were performed only after the primary environmental result was opened and therefore cannot alter its frozen gate.

| Sensitivity | Species / strata | Median BIO5 contrast | Wilcoxon p | Conditional OR | Conditional p | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| Observer-paired | 106 species / 144 paired species-observer strata | 0.000 SD | 0.485 | 0.787 (95% CI 0.544–1.138) | 0.202 | null under strongest observer control |
| Observer-balanced | 352 species | +0.0541 SD | 0.0750 | — | — | positive direction, weaker statistical support |

The observer-balanced sign test was p = 0.0484. Together, these sensitivities show that the prospective BIO5 association is not invariant to observer conditioning.

Machine-readable receipt:
- `results/polymorphism_white_environment_observer_sensitivity_20260925/result.json`.

### Species-disjoint BIO5 transport

| Cohort | Eligible species | Median white-minus-nonwhite BIO5 contrast | Species-level p | Conditional OR | Conditional p | Primary support |
|---|---:|---:|---:|---:|---:|---|
| Discovery | 271 | -0.0068 SD | 0.743 | 1.033 | 0.131 | false |
| Validation | 260 | +0.0662 SD | 0.0541 | 1.048 | 0.0319 | false |

Frozen transport decision:

**Not supported under the prespecified cross-cohort rule.**

Interpretation boundary:

The frozen primary prospective analysis contains a within-cohort association between white states and warmer BIO5 environments, but the association weakens under observer controls and does not transport as a common rule across the original species-disjoint cohorts. This is consistent with context-dependent environmental sorting and does not establish causal heat selection, a universal temperature effect, or independence from the known exposure coupling of the white classifier.

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

D summarizes species-wide sampled four-state colour diversity and can combine local coexistence with differentiation among sampled locations. Greater D is associated with stronger within-species geographic colour organization in both species-disjoint cohorts. This is a replicated structural correlate, not a variance partition and not evidence that geographic organization causes D or that any particular climate, pollinator, demographic, gene-flow or selection mechanism has been identified. The frozen vertex-permutation null tests geographic organization against an unstructured assignment; an equal-strength spatial-process measurement-opportunity question is treated separately as a post hoc diagnostic below.

### Post hoc equal-structure measurement-opportunity diagnostic

The frozen D–spatial test does not ask whether species with different D values would yield different measured spatial rho if their **true spatial ordering strength were held constant**. We therefore added a post-outcome diagnostic on the 363-species validation cohort. This diagnostic does not alter any frozen decision.

The diagnostic first reproduced the frozen validation quantities directly from the photo-level measured table:

- mean species spatial rho = **0.0254826061**;
- Spearman(D, species spatial rho) = **0.1016008447**.

For every species, the diagnostic then preserved the observed coordinates and the exact multiset of observed four-state colour vectors. Colour vectors were ordered by the first classical-MDS coordinate of their pairwise Jensen–Shannon dissimilarities and assigned along the first spatial principal axis plus deterministic Gaussian rank noise. The same noise parameter k was applied to every species within a replicate; larger k gives weaker common spatial ordering.

| Equal-structure setting | Replicates | Mean of species-mean rho | Mean D–rho | D–rho range | Replicates >= observed D–rho |
|---|---:|---:|---:|---:|---:|
| k = 4 | 100 | 0.02890 | -0.0359 | -0.1454 to +0.05371 | 0 |
| k = 6 | 25 | 0.01350 | -0.0197 | -0.07924 to +0.04803 | 0 |
| k = 10 | 25 | 0.00423 | -0.0178 | -0.1592 to +0.1081 | 1 |

The **k = 4** setting is the primary diagnostic comparison because its mean spatial rho (0.02890) is closest to the observed validation mean (0.02548). None of its 100 replicates reached the observed cross-species D–rho of 0.1016. Across all sensitivity settings, one of 150 replicates reached or exceeded the observed association.

This result argues against the narrow explanation that the observed positive D–spatial association arises simply because different D values create different measurement opportunity under this linear-gradient equal-structure process. It is **not** a preregistered test, does not justify a universal claim that high-D species are “disproportionately geographically partitioned”, and does not cover patchy, threshold, multimodal or environmentally warped spatial structures. The frozen discovery/validation D–spatial decisions remain unchanged.

Canonical implementation and receipt:
- `scripts/analysis/run_polymorphism_D_spatial_equal_structure_posthoc_20261007.py`;
- `results/polymorphism_D_spatial_equal_structure_posthoc_20261007/result.json`.

The reporting receipt records the exact measured-table SHA256, canonical workflow run **37564518153**, artifact **11457707244**, artifact digest and the fact that confirmatory decisions are unchanged.

## S5a. Distributed polymorphism: local diversity is depleted relative to species-wide composition

Canonical reporting receipt:

- `results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json`

The analysis was defined after the earlier D–spatial result and is therefore post hoc. It tests a different estimand: whether nearby observations are more homogeneous than expected **conditional on each species' exact overall colour composition**.

### Primary 50-km result

| Cohort | Evaluable species | Species-wide pair diversity | Local pair diversity | Mean depletion | Matched p |
|---|---:|---:|---:|---:|---:|
| Discovery | 166 | 0.22139 | 0.20087 | **0.02053** | **0.005** |
| Validation | 181 | 0.22866 | 0.20999 | **0.01867** | **0.005** |
| Third cohort | 204 | 0.24820 | 0.23352 | **0.01468** | **0.005** |

The matched vertex null fixes coordinates, sample size and every species' exact coarse-state counts. Therefore the positive depletion does not arise simply because a high-D species has more possible mismatched pairs.

Mean local depletion was also supported at **25, 100 and 250 km in all three cohorts**.

### Direct falsification tests at 50 km

| Test | Discovery | Validation | Third cohort |
|---|---:|---:|---:|
| Different-observer pairs only | 0.01998, p=.005 | 0.01690, p=.005 | 0.01457, p=.005 |
| Nonwhite-only three-state | 0.01947, p=.005 | 0.00729, p=.025 | 0.00891, p=.040 |
| Continuous nine-colour JSD | 0.01518, p=.005 | 0.01366, p=.005 | 0.01301, p=.005 |
| Quarter-stratified composition null | 0.02287, p=.005 | 0.02061, p=.005 | 0.01560, p=.005 |
| Cross-year local pairs | 0.02095, p=.005 | 0.01437, p=.005 | 0.01377, p=.005 |

Thus geographic local homogeneity is not restricted to same-observer comparisons, the exposure-sensitive white state, the four coarse categories, calendar-quarter composition or same-year sampling.

The stronger cross-species claim that high-D species are always more strongly locally partitioned is **not** supported at the fixed 50-km discovery–validation criterion: discovery rho=0.0390 (p=.135), validation rho=0.1171 (p=.025), third rho=0.1118 (p=.005). That relationship is scale dependent and is not promoted to a primary claim.

Interpretation boundary:

> The replicated result is **distributed polymorphism**: species-wide flower-colour variation is allocated among geographic localities. The test does not identify selection, gene flow, drift, history or plasticity as the cause.

## S5b. Phenotypic IBD versus IBE-like decomposition

Canonical reporting receipts:

- `results/polymorphism_phenotypic_IBD_IBE_posthoc_20261007/result.json`
- `results/polymorphism_BIO5_colour_sorting_robustness_20261007/result.json`
- `results/polymorphism_selected_environmental_sorting_transport_20261007/result.json`

For each species, continuous nine-colour Jensen–Shannon dissimilarity was related symmetrically to geographic distance and absolute BIO5 difference.

| Cohort | IBD-like: colour ~ geography \| BIO5 | IBE-like: colour ~ BIO5 \| geography | IBE - IBD |
|---|---:|---:|---:|
| Discovery | **0.02760**, p=.005 | **0.00925**, p=.005 | -0.01835 |
| Validation | **0.02414**, p=.005 | **0.00826**, p=.010 | -0.01589 |
| Third cohort | **0.02225**, p=.005 | **0.01019**, p=.005 | -0.01206 |

Both components are positive in all three species-disjoint cohorts, but the IBD-like component is larger in every cohort. Descriptive commonality partitioning likewise assigns more unique rank variance to geography than to BIO5.

### Technical ceiling on the IBE-like component

The positive BIO5 residual is reproducible in continuous flower-colour measurements, but stronger technical controls do not establish it as a flower-specific environmental response:

- validation flower-minus-background: mean partial rho = **-0.00328**, p = **0.82**;
- same-observer-pair discovery: -0.1538, p=.995;
- same-observer-pair validation: -0.0737, p=.89;
- same-observer-pair third cohort: -0.0633, p=.83.

Accordingly, the manuscript uses **phenotypic IBE-like** only as a bounded diagnostic. It is not genetic IBE, and it does not demonstrate local adaptation.

### Environmental heterogeneity does not generally explain total D

A multivariate climate-heterogeneity prediction did not replicate from discovery to validation. Solar-radiation heterogeneity was positive in the 500+500 screen, but fixed transport to the third cohort failed (partial rho = **0.0171**, p = **0.382**). The paper therefore does not infer that broad climatic heterogeneity generally generates species-wide flower-colour diversity.

### Relationship to the earlier signed BIO5 result

The earlier post-confirmatory prediction that warmer environments should specifically favour the white state passed within the third cohort but weakened under observer controls and failed the frozen discovery/validation transport rule. The newer unsigned BIO5-distance result asks a different question: whether colour changes **somehow** along temperature differences, irrespective of direction. Even that result remains technically bounded by the controls above.

## S5c. Verified post-freeze fourth-cohort 199,700-photo extension (17 October 2026)

**Historical snapshot only:** the 447,500-photo data were unfinished at the 17 October intermediate freeze described below; the **actual completed five-cohort result is now authenticated in S5f**. Preserve the old description to document chronology, not as current status.

This is a **later fixed-method replication** and must not be merged chronologically with the original 1,499-species prospective-confirmation protocol. Source authority: [real 199,700-photo original-method release](https://github.com/zuizui0223/fcp/releases/tag/fcp-real-original-199700-wave4-50km-H1-W-20261017), its SHA-pinned release receipt, and [the explicit final-U100/pending-status evidence ledger](FCP_U100_REAL_199700_EMPIRICAL_EVIDENCE_AND_447500_PENDING_20261017.md). The [actual wave4 ROI run #38992543115](https://github.com/zuizui0223/fcp/actions/runs/38992543115) finished **256/256 genuine image workers**. The [original-method spatial append replay #39023781282](https://github.com/zuizui0223/fcp/actions/runs/39023781282) matched the frozen original three-cohort effects.

| Cohort | Measured original photographs | True 50-km eligible species | Species-equal global-minus-local pair diversity | Upper p (199 vertex nulls) |
|---|---:|---:|---:|---:|
| Discovery (unchanged) | 50,000 | 166 | 0.0205293 | 0.005 |
| Validation (unchanged) | 50,000 | 181 | 0.0186730 | 0.005 |
| Prospective cohort (unchanged) | 49,900 | 204 | 0.0146849 | 0.005 |
| **Fourth wave (newly measured)** | **49,800** | **190** | **0.0172936** | **0.005** |
| **Four-cohort pooled** | **199,700** | **741** | **0.0176372** | **0.005** |

Eligibility is **at least 40 biologically classifiable photographs and 30 local photograph pairs** per species. Geographic separation is within 50 km, a radius selected *post hoc* before these cohort outcomes were opened. The null preserves each species' full colour-state counts and its original coordinate graph and uses 199 complete-vertex colour-label assignments (smallest possible plus-one p = 0.005). The pooled mean global pair diversity was **0.2342544**; local pair diversity was **0.2166172**. This is within-species *photographic* colour structure rather than population morph frequency or selection.

Other unchanged frozen analyses in the genuine fourth wave:

| Analysis | Fourth-wave result | Inferential boundary |
|---|---|---|
| Observer-disjoint D, 200 splits | Median Spearman rho **0.805596**; 5th–95th percentile **0.775457–0.832733**; median paired species **316**, median CCC **0.844141** | Same-provider measurement reliability, not unbiased wild population morph frequencies |
| Fixed white-axis W, primary second-mode >=0.10 | n **134**; W **0.536816** vs structured-null median **0.436254**; upper p **0.001** | Squared axis alignment only, not directional whitening |
| Fixed white-axis W, strict second-mode >=0.20 | n **67**; W **0.547280** vs structured-null median **0.458086**; upper p **0.001** | Preserves coarse-state structure; 999 nulls |
 
Neither the H1 replication nor the W result is a newly preregistered prospective test. The original three-cohort nonwhite-only, different-observer, calendar-quarter, cross-year, continuous-distance and 25/100/250-km sensitivities were **not all rerun for the fourth wave**; only the exact 50-km fixed-method comparison and listed H1/W components are confirmed here.

**Unmeasured fifth cohort — do not borrow results:** an independently source-verified 4,475-species / **447,500-source-photograph** analysis-eligibility ledger was frozen after whole-taxon exclusion of three genuine source-collision taxa. Source eligibility alone is not actual ROI-v4 colour measurement. At the 17 October evidence freeze, the new-U100 pixel job [#38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400) remained in production, and final five-cohort combined null/H1/W estimates were **not available**. GitHub Actions workflows can display aggregate success with downstream scientific jobs skipped; only authenticated completed pixel terminals plus executed inferential jobs authorize a final biological claim.

## S5d. Post hoc species-level leverage of the replicated 50-km effect (17 October 2026)

**Question:** does a repeated positive *cohort mean* imply that most species have comparably strong local photographic colour clustering? The independently run [real-photo species-leverage audit](https://github.com/zuizui0223/fcp/actions/runs/39042459866) recovered the unchanged, [SHA256-frozen 199,700-photo original-method release](https://github.com/zuizui0223/fcp/releases/tag/fcp-real-original-199700-wave4-50km-H1-W-20261017), verified its species-level table and the previously frozen 50-km mean and permutation p separately for all four cohorts and the pooled result, then described the distribution of the original species-specific depletion values. Machine-readable audit: `results/fcp_real_199700_spatial_species_leverage_20261017.json`. **No new label permutations, radius selection, inferential p-values, or image classification were performed.**

| Cohort | Species with positive depletion / evaluable species | Original mean depletion | Median species depletion | Mean after dropping largest 10% effects | Mean after dropping largest 20% effects |
|---|---:|---:|---:|---:|---:|
| Discovery | 91 / 166 (54.8%) | 0.02053 | 0.00580 | 0.00326 | -0.00406 |
| Validation | 104 / 181 (57.5%) | 0.01867 | 0.00740 | 0.00098 | -0.00739 |
| Prospective | 113 / 204 (55.4%) | 0.01468 | 0.00564 | 0.00078 | -0.00779 |
| Fourth cohort | 109 / 190 (57.4%) | 0.01729 | 0.00556 | 0.00250 | -0.00656 |
| All four pooled | 417 / 741 (56.3%) | 0.01764 | 0.00597 | 0.00191 | -0.00655 |

Of **741** species meeting the unchanged >=40 classifiable-photo and >=30 within-50-km local-pair gates, **417** (56.3%) had positive global-minus-local pair dissimilarity, **61** (8.2%) exactly zero and **263** (35.5%) negative. The pooled original mean was **0.0176372**, whereas the median per-species depletion was only **0.005974**. Removing the top 10% largest positive species-level effects (75 species) leaves a descriptive mean **0.001911** among 666 species. Removing the top 20% (149 species) makes that descriptive mean **−0.006553** among 592 species.

**Inference boundary:** These effect-tail exclusions were devised *after seeing the pooled positive result*. The raw remaining means in this first audit did not themselves have matched-null p-values; Section S5e now supplies an explicitly post hoc test that reapplies each exclusion separately inside every original permutation world. Neither stage is independent confirmation. The original p=0.005 in every cohort and the pooled result continues to describe a departure in the **species-equal-weighted cohort mean under the original null**; it does not establish that the majority of species individually show statistically supported clustering. All four cohorts are from the same iNaturalist photographic source. The observed mixture of positive, zero and negative species effects limits universal claims and motivates species-specific mechanisms and quality-adjusted field validation. The fifth U100 cohort is not included here.

## S5e. Robust original-null replay: removing large-effect species does not erase excess spatial colour sorting (18 October 2026)

**Inferential correction to S5d.** A negative observed mean after dropping the species with the largest positive effects does *not* imply that the remaining species are null-equivalent: deleting top-tail observations also shifts the matched null statistic downward. We therefore re-ran the **exact unchanged 199 original within-species vertex-label permutations** from the [SHA256-frozen genuine 199,700-photo release](https://github.com/zuizui0223/fcp/releases/tag/fcp-real-original-199700-wave4-50km-H1-W-20261017), using the historical photo-ID ordering, cohort-specific seed rules, 40 classifiable-photo and 30 within-50-km pair gates. The original four 50-km species-equal mean depletion statistics and the pooled effect, with their original p = 0.005, were reproduced numerically **before** computing new statistics. [Authenticated execution and result artifact](https://github.com/zuizui0223/fcp/actions/runs/39045286620) (full output SHA256 `543e4a544310086b800d31e410eed176240c88eb7ce2544a25eaef324488547c`).

Four *exploratory* distributional estimands were computed for each cohort and for all four pooled: (i) species-level depletion median, (ii) fraction of species with positive depletion, (iii) mean after removing the highest 10% of species-specific effects, and (iv) mean after removing the highest 20%. **Every one of the 199 null worlds independently re-selects its own highest-effect species for removal**; the observed top-effect species were *not* carried over into the null. One-sided upper p is `(1 + number of null statistics >= observed statistic)/200`, with minimum 0.005. Holm step-down adjustment was applied **within these four newly examined robust tests per cohort**. This does not correct all earlier exploratory analyses, and the source remains one citizen-science image ecosystem.

| Pooled 50-km statistic (741 eligible species) | Observed | Matched null mean | Upper p | Within-four Holm p |
|---|---:|---:|---:|---:|
| Species median depletion | **+0.005974** | −0.000010 | 0.005 | **0.020** |
| Fraction of positive-depletion species | **0.56275** (417/741) | 0.45007 | 0.005 | **0.020** |
| Drop highest-effect 10% (75 species) | **+0.001911** | −0.011139 | 0.005 | **0.020** |
| Drop highest-effect 20% (149 species) | **−0.006553** | −0.018278 | 0.005 | **0.020** |

The four species-disjoint, same-iNaturalist-source cohorts gave the following exploratory within-four Holm p-values:

| Original photo cohort (50-km eligible species) | Median | Positive fraction | Drop top 10% | Drop top 20% |
|---|---:|---:|---:|---:|
| Discovery (166) | 0.030 | 0.030 | 0.020 | 0.020 |
| Validation (181) | 0.020 | 0.020 | 0.020 | 0.020 |
| Prospective (204) | 0.020 | 0.020 | 0.020 | 0.020 |
| Fourth wave (190) | 0.020 | 0.020 | 0.020 | 0.020 |

**Ecological inference ceiling:** The broad rightward shift of the *species-level depletion distribution*, compared with species-conditioned permutations, is not explained exclusively by the largest observed positive effects. Yet a substantial number of species (263 of 741) had negative individual depletion, and this analysis cannot say how many species would individually pass a valid, multiplicity-corrected test. Nor can it establish wild-population morph frequencies, adaptive divergence, pollinator selection, genotype frequencies, or any evolutionary transition direction. All four distributional tests, their chosen trimming levels, and the 50-km radius are post-outcome, and Holm correction was not applied across all hypotheses ever considered in the project.

## S5f. Actual final five-cohort 447,500 original photo measurements: fixed-method H1, geographic colour partitioning and W (18 October 2026)

**This section supersedes only the pending-data status described in the historical 17 October snapshot in S5c; it does not overwrite any earlier frozen test result.** [Actual source-pixel original-method Release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017); [completed genuine image ROI #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400); [successful unchanged original five-cohort science #39049598596](https://github.com/zuizui0223/fcp/actions/runs/39049598596); [independently audited permanent release #39050465077](https://github.com/zuizui0223/fcp/actions/runs/39050465077). Source-verified machine-readable summaries and immutable release asset SHA256 are indexed by \`results/fcp_original_447500_five_cohort_verified_reporting_20261018.json\` and \`docs/FCP_FINAL_447500_REAL_FIVE_COHORT_RESULTS_20261018.md\`.

The original outcome-blind U100 opportunity frame comprised 4,730 taxa. **2,481** new H9 taxa actually completed 100 authentic photo-ID source assignments each and all **256 real blind ROI-v4/EfficientSAM partitions** succeeded. Three newly sampled taxa genuinely reused historical source photo/observation IDs; following frozen conservative rules they were excluded *whole*, not patched with substitute photos. A further 252 taxa lacked 100 genuine original photos. Together with the 1,499 old discovery/validation/prospective taxa and 498 real wave4 taxa, the valid analysis therefore comprised **4,475 distinct taxa × 100 actual original photographic terminals = 447,500**. Of these, **227,107** photos were classifiable into the four fixed coarse colour states, and **220,393** were nonclassifiable. The original D classifiability gate included 3,297 species; the primary spatial 50-km null required an additional >=30 local pairs, leaving **1,761** taxa. **Do not substitute 4,475 or 3,297 for the 50-km test denominator.**

| Source cohort | Genuine 100-photo taxa | Actual original photo terminals | Four-state classifiable rows | 50-km eligible taxa | Mean global-minus-local dissimilarity | Original upper p (199 permutations) |
|---|---:|---:|---:|---:|---:|---:|
| Discovery (old untouched) | 500 | 50,000 | 25,377 | 166 | 0.020529255 | 0.005 |
| Validation (old untouched) | 500 | 50,000 | 24,885 | 181 | 0.018672972 | 0.005 |
| Prospective (old untouched) | 499 | 49,900 | 25,788 | 204 | 0.014684920 | 0.005 |
| Fourth-wave unchanged-method | 498 | 49,800 | 24,585 | 190 | 0.017293596 | 0.005 |
| **New fifth-wave U100** | **2,478** | **247,800** | **126,472** | **1,020** | **0.017598698** | **0.005** |
| **All five source cohorts** | **4,475** | **447,500** | **227,107** | **1,761** | **0.017614902** | **0.005** |

The original 50-km global mean pairwise dissimilarity was **0.2243045205**, versus the within-50-km local mean **0.2066896183**. The frozen species-equal-weighted pooled null mean was **0.0000445943**, central 95% interval **−0.0020032432 to 0.0024735998**, and upper Monte Carlo p **0.005**. The genuine fifth-wave subset alone (n=1,020) had depletion **0.0175986977**, upper p **0.005**, close to the four-cohort pooled depletion **0.0176372079**. The new-wave 25-/100-/250-km sensitivity means were also positive with upper p=0.005 (n=710, 1,361 and 1,772, respectively). **The 50-km radius remains post hoc**, and additional radii are post hoc tests, not alternative optimized primary scales.

The entire study was also replayed through the **unchanged original 200-partition observer-disjoint H1/D test**: median Spearman rho **0.8035357424** (5th–95th percentiles **0.7931469330–0.8094055921**), median paired D species **3,044**, median CCC **0.8553147636**, median Spearman–Brown **0.8910671671**. H1 describes repeatability of photograph-sample species-level D, not unbiased natural morph frequencies.

The frozen **white-versus-eight-other-colour nine-component axis** and original **999 coarse-state-preserving structured null worlds** gave, for the combined actual five-cohort data:

| Frozen second-mode threshold | Vector-eligible species | Observed squared W | Structured-null median | Upper p |
|---|---:|---:|---:|---:|
| Primary 0.10 | **1,223** | **0.5187948657** | **0.4553650050** | **0.001** |
| Strict 0.20 | **658** | **0.5254168292** | **0.4612009436** | **0.001** |

The original **499-species** third-cohort **untouched prospective confirmation** of the fixed white axis remains separate: primary n=158, W=0.5172457461 versus null median 0.4571428150, p=0.001. The later 4,475-species pooling cannot retrospectively upgrade this original confirmatory test, and squared W is **sign-invariant**: it does *not* indicate flower-colour evolutionary direction, anthocyanin reactivation or a synthesis-cost fitness trade-off.

**Limits:** these are photographed, location-blind ROI outputs from the *same iNaturalist source and fixed classifier*, not independent photographic-provider evidence or genotype/population morph frequency estimates. Exactly which selection/history/plasticity/pollinator mechanisms generated local sorting remains unresolved. The original three-cohort battery included continuous nine-colour and background-specific controls that were **not** identically rerun for every new fifth-wave taxon. However, the simpler **different-observer, cross-year and nonwhite-only restrictions were explicitly extended across the completed five measured cohorts in S5g**; the nonwhite-only contrast is not supported in every separate cohort. The original-four-cohort robust-tail tests in S5d/e were not independently repeated on all 2,478 new taxa.

## S5g. FCP commonality across the full five genuinely measured cohorts: morph dominance and geographic allocation (18 October 2026)

This explicitly **post hoc** analysis asks about general *within-species FCP architecture* rather than global proportions of named flower hues. Authentic source is the [447,500 real original photographs](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). Independent execution [#39091197295](https://github.com/zuizui0223/fcp/actions/runs/39091197295) reproduced the original 50-km species-equal mean/null numerically, verified original source hashes, passed regression tests, and generated the complete species-level table and source-backed commonality result (artifact \`fcp-447500-real-polymorphism-commonality-source-backed\`). Full output JSON SHA256: **\`3a51d537fd99766e64846ad318c53d3e70f1d1fe0a1222ecfcae2c6854e551e2\`**. A compact machine-readable source-identity summary is \`results/fcp_original_447500_polymorphism_commonality_20261018.json\`; the numerical/biological interpretation ledger is \`docs/FCP_WITHIN_SPECIES_POLYMORPHISM_COMMONALITY_447500_20261018.md\`.

**Species-level composition measurement:** Among 3,297 species with >=40 four-state classifiable original photographic rows (227,107 classifiable rows out of 447,500 photo terminals), Gini–Simpson \`D=1-sum(p_k²)\` had median **0.172363** (IQR **0.063424–0.392938**). The median photographic share of the most common coarse colour state was **0.905660**, with 51.4% of species having >=90% of their classified photos in one coarse state. The counts of species with 1, 2, 3, 4 observed coarse states were **268 / 1,411 / 1,080 / 538**. This is **photograph-sampled skew**, NOT genetically validated FCP prevalence, true wild population morph frequency or evidence for pigment-cost selection.

**Original geographic null and identity-free extension:** The original 1,761 taxon 50-km local-pair test reproduced pooled mean global-minus-local four-state diversity **0.017614902**, with upper p=**0.005** under the exact original **199 complete-vertex within-species label permutations**. Replacing named states *within each species* by the **most frequent photographed state versus all other states** preserved its sample composition/coordinates/null worlds and gave **0.016485145**, p=**0.005**. Genus-equal weighting of 832 inferred genera gave four-state **0.015251210** (p=0.005) and majority-versus-alternatives **0.014632192** (p=0.005). Genus equality is not a phylogenetic covariance model.

| Photographically most frequent state in a species | Eligible species | Original four-state mean 50-km depletion | Upper p |
|---|---:|---:|---:|
| White | 807 | 0.018551954 | 0.005 |
| Yellow/orange | 521 | 0.010326010 | 0.005 |
| Red/pink | 203 | 0.025927013 | 0.005 |
| Blue/purple | 230 | 0.023501657 | 0.005 |

The species-global D rank quartiles (n=441/440/440/440) also showed positive mean original 50-km depletions **0.002904, 0.011920, 0.018917, 0.036752** respectively (all upper p=0.005). The much smaller effect in low-D species follows limited variation available to be segregated and does *not* identify a selective cause. All stratification choices and grouped comparisons are **post hoc**; the 50-km radius itself was originally selected after relevant outcomes. Subgroup tests from the same photographed data and colour label assignments are not independent replication or a project-wide multiplicity-corrected family.

**Same-source falsifications now extended to the fifth cohort (distinct from the more extensive original three-cohort technical battery):**

| Restriction | Eligible species | Species-equal mean depletion | Upper p |
|---|---:|---:|---:|
| Original full four-state 50-km test | 1,761 | 0.017614902 | 0.005 |
| Exclude local pairs with the same observer | 1,745 | 0.017070095 | 0.005 |
| Retain cross-year local photo pairs only | 1,643 | 0.015844575 | 0.005 |
| Retain species with >=60 classifiable photos | 845 | 0.017065745 | 0.005 |
| Remove all white photographed states; three nonwhite states only | 763 | 0.006594523 | 0.005 |

**Critical nonwhite exception:** Although the *pooled* five-cohort nonwhite-only analysis is positive, the original prospective cohort (n=78) did not support that analysis (upper p=**0.14**), and the fourth-cohort n=79 result was negative (**−0.002277**, upper p=**0.72**). The other cohort estimates also vary. Thus the effect cannot honestly be described as a universal chromatic-only geographic pattern. This photo-state/observer/year filter set does not test the full flower-minus-background ROI context, local plasticity, field morph frequencies or historic gene flow. White classification still has a measured clipping/exposure association.

**Conclusion and hard ceiling:** FCP is best summarized as species-conditional **(i) skewed photographed alternative-state composition** combined with **(ii) above-null average local spatial sorting, even when specific hue labels are discarded**. This is a comparative regularity in a nonrandom citizen-science photograph sample, not a known balancing-selection equilibrium, a recurrent genetic loss/re-expression process, directional anthocyanin evolution, a demonstrated pigment fitness trade-off or a population demographic mechanism. The separate originally frozen, untouched prospective H2 W confirmation and earlier climate/phylogenetic boundaries are unchanged by the new synthesis.

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

## S8. Bounded secondary mechanism evidence

This section deliberately houses evidence that is **not part of the main empirical Results spine**. The *Silene littorea* transcriptomic/HPLC evidence is published work used as mechanistic context, whereas the PAL/WAL frequency comparison is the present study's descriptive reanalysis of published tables. Main text cites their bounded interpretation in Discussion only. Neither layer upgrades the prospective H2 test or provides a direct molecular bridge to the cross-species H2 geometry.

### Secondary *Silene littorea* molecular anchor

Machine-readable structured extraction:
- `results/polymorphism_silene_molecular_anchor_20260928/result.json`.

Source:
- Casimiro-Soriguer et al. (2016), doi:10.3389/fpls.2016.00204.

The current study did not re-run the raw RNA-seq. It deterministically extracted the source-reported quantitative evidence most directly relevant to petal-specific anthocyanin loss: all p<0.05 bud-stage pigmented-versus-white ABP expression contrasts, significant dark-versus-light regulatory contrasts, the expanded sequence follow-up and petal HPLC results.

Key extracted results:
- F3h1 dark/white = **49.0×**, p = **0.039**;
- F3h1 light/white = **42.2×**, p = **0.049**;
- Myb1a dark/white = **5.1×**, p = **0.009**;
- F3h1 is the only locus significant in both pigmented-versus-white bud contrasts;
- the reported F3h1 sequence table contains zero SNPs;
- expanded sequencing of 38 individuals did not identify a SNP consistently differentiating colour morphs;
- petal HPLC differences in rutin, quercetin and isovitexin are concordant with a blockage near F3h1.

This provides a source-derived single-species molecular anchor for the PAL interpretation. It does not establish a causal Myb1a mutation or the molecular basis of the cross-species H2 white axis.

### Secondary PAL/WAL maintenance evidence

Evidence-status ledger: `docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md`.

- Within-*Silene littorea* persistence receipt: `results/polymorphism_silene_decoupling_persistence_20260925/result.json`.
- Cross-species interval/censoring-aware receipt: `results/polymorphism_crossspecies_pal_wal_frequency_20260925/result.json`.
- Source: Del Valle et al. (2019), Supplementary Tables S2 and S1, respectively.
- Biochemical anchor: in *S. littorea*, HPLC-DAD-MS^n profiling of petals, calyces, leaves and stems distinguished petal anthocyanin-loss (PAL; anthocyanins absent from petals but retained in photosynthetic tissues) from whole-plant anthocyanin-loss (WAL; anthocyanins absent from petals and photosynthetic tissues). The *S. littorea* frequency contrast is therefore biochemically anchored rather than based on visual white-flower labels alone. The cross-system Table S1 entries retain the source paper's PAL/WAL classifications; those systems were not independently re-profiled by the current study.

| Contrast | PAL | WAL |
|---|---:|---:|
| *S. littorea* positive-frequency median | 15.5% | 0.21% |
| *S. littorea* positive-frequency range | 8–21% | 0.05–0.86% |
| Cross-system instances | 13 | 13 |
| Cross-system median bound used for comparison | lower bound 5% | numeric upper bound 0.1% |
| Largest quantified WAL upper bound | — | 1.4% |

In *S. littorea*, every positive PAL frequency exceeded the maximum positive WAL frequency. Across the 13+13 literature systems, 7/13 PAL lower bounds exceeded the 1.4% maximum quantified WAL upper bound, and every PAL upper endpoint exceeded it. These analyses are descriptive reanalyses; literature ascertainment, heterogeneous survey effort and censored/qualitative frequency reporting preclude an unbiased meta-analytic or causal interpretation.

A direct PAL/WAL-to-FCP bridge was also audited rather than assumed. Taxonomy-aware resolution preserved the source strings while matching accepted names. The natural-frequency subset overlaps only four D-eligible reserve species: PAL *Gymnadenia rhellicani* (D = 0.6515) and *Silene gallica* (D = 0.2344), plus WAL source *Delphinium nelsonii* resolved to *D. nuttallianum* (D = 0.1659) and source *Silene dioca* resolved as an orthographic candidate to accepted *S. dioica* (D = 0.2055). Discovery overlap and prospective selected overlap are both zero. *Mimulus guttatus* resolves to *Erythranthe guttata* (D = 0) but its source frequency is greenhouse-only and is excluded from the natural-frequency bridge. With only two PAL and two natural-context WAL overlaps, all in one reserve cohort, no inferential PAL-versus-WAL comparison of D or H2 geometry is estimable. Machine-readable receipt: `results/polymorphism_pal_wal_h2_overlap_20260928/result.json`. Taxonomic crosswalk: `docs/POLYMORPHISM_PAL_WAL_TAXONOMIC_ALIAS_AUDIT_20260928.md`.

## S9. Canonical main-text figures

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

## S10. Manuscript claim guard

Claim-guard test:

- `tests/test_polymorphism_manuscript_claims.py`

Workflow:

- `.github/workflows/polymorphism-manuscript-claim-guard.yml`

Most recent verified post-architecture README run:

- workflow run: `35299258800`
- conclusion: success.

The guard checks the current manuscript, claim ledger, figure plan and README against the frozen prospective-confirmation values and required claim boundaries.

## S11. Literature-positioning audit

Bounded manuscript literature audit:

- `docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`

This literature layer supports context and interpretation only. It cannot modify the machine-readable empirical verdicts.

## S12. Hard nonclaims carried into all supplementary material

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
