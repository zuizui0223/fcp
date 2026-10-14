# Flower-colour variation — one unified New Phytologist Supporting Information map (14 October 2026)

This is the single Supporting Information index for the integrated New Phytologist Full Paper. The original high-depth evidence and its frozen release receipt are preserved; the newly added sparse worldwide colour–environment source has its own verified numerical ledger and different statistical grain.

A direct claim-to-input routing table is maintained at `docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`.

A self-contained verification bundle is maintained under release tag `fcp-np-provenance-20260926` (https://github.com/zuizui0223/fcp/releases/tag/fcp-np-provenance-20260926). The canonical Git-tracked receipt, `archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`, records the source commit, asset byte count, SHA256 and packaged-file count for the currently published asset.

It does not recompute any result. It organizes frozen protocols, machine-readable outputs, figures and provenance into a manuscript-facing evidence map.

## How to read this Supporting Information

The reader-facing order mirrors the unified paper's inferential chronology, with original S1–S12 retained and new S13–S16 added:

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

## S13. Sparse worldwide photo-environment evidence: source and response denominators

### Exact repeat-photo availability and why estimands have different species counts

The independent outcome-blind photograph-capacity scan in `main/docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md` establishes these observer-capped source ceilings: at least 1 retained photo for 41,322 taxa; at least 10 for **18,301**; at least 20 for **12,985**; at least 40 for **8,753**; at least 60 for **6,770**; at least 80 for **5,558**; and at least 100 for **4,730**. These counts were recorded before phenotypic classification. They explain why all approximately 20,000 cross-species flower-colour/environment units cannot automatically supply 100 independent individual photos and the same W/minor-cluster and within-locality photo-pair estimands.

Retrospectively, original flower-image pixels had **already been measured at sparse depth for 42,111 nominal species (18,457 classifiable single-photo outcomes)** and 85,337 species×geographic-cell states (39,075 classifiable). These coarse colour classifications are real measurements, not merely unprocessed photo-ID metadata, but do not furnish 100 images per species or a reliable species-specific two-colour axis. Original-image synthesis: FCP PR #134/#135, verified archived source commit `2b5390ea74f8d196012499d2d35dd477f9795938`. A failed prior P500 allocated 500 names, 499 measurement attempts, but no admissible prospective H2 conclusion. The full historical 4,730-species ≥100-photo opportunity frame had 2,000 identities allocated (1,000 discovery/validation + 500 failed P500 + 500 selected third), leaving **2,730 truly unused**. Another **2,000 selected species and 730 reserve taxa** have since been frozen on PR #133 as *metadata only* with no new biological image measurement. Therefore the paper's 1,499 is **admitted high-depth sampled species**, not all measured flower-colour species worldwide and not the physical or historical maximum.

The global genus-conditioned BIO12 and vapour-pressure white-photo slopes use **20,694 named source species in 1,934 source genera**. This is smaller than the 25,028 species with climate-complete four-colour records, and differs from the genus/geography prediction heldout species count (17,183 climate-only or 13,961 all36). At the other end, the 1,499 species in high-depth resources have 100 sampled photographs each; 377 of 499 new cohort species satisfy the frozen >=40 classifiable-photo gate, 158 satisfy the prospective achromatic W-vector construction and 166/181/204 across the three cohorts qualify for the specified 50-km local-pair spatial estimand.

A minimum of ten *original observed photos* is not ten *classifiable flower organs* or a confirmed second colour cluster, and even 100 photographs do not guarantee enough distinct local plants, multiple sites or a valid colour minor-mode vector. Pooling low-count taxa toward a genus could estimate broad class probabilities with model uncertainty, but unobserved biological morphs must be marked **UNRESOLVED**. Full capacity, detection probabilities, and a non-confirmatory tiered modelling extension are separated in `docs/POLYMORPHISM_20K_SCALING_SUPPORT_20261014.md`. No new W or local geographic inference at 18,301 species has been carried out.


The additional worldwide macroecological comparison is **not a fourth high-depth hypothesis-test cohort**. It uses the original 100,543 photographed opportunities of 42,111 named species (45,953 four-colour classifiable photos; 54,590 classified as unknown). Of these, 45,836 photos of 25,028 taxa have all 26 nonsoil predictors; 31,925 photos of 20,786 taxa additionally have all ten soil predictors. These are labelled photo, species and covariate denominators rather than true global flower-colour population frequencies.

The 36 measured variables are organized as 11 temperature, eight precipitation, two solar-radiation, two wind, two atmospheric water-vapour, one elevation and ten soil variables. Genus-heldout photo-colour composition prediction gained Brier +0.002943 with 26 predictors (17,183 test species) and +0.001432 with all 36 (13,961 test species) using genus-equal weighting and different complete-case sets. Those are **prediction differences, not directional fitness or selection coefficients**.

Source photo archive SHA256: `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`. Original 36-abiotic archive SHA256: `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`. Archived parent analyses and their seven-block score interpretation: FCP PR #165–169, including Shapley model #166; the unified manuscript uses prediction only as background, not a new select-one-variable causal result.

**Supporting Information Fig. S12** will visualize source representation, outcome-classification missingness, seven correlated environmental blocks, genus/region holdout and the different sets of eligible taxa; Figure 6a contains the two headline Brier score improvements.

## S14. Signed colour directions, photographic quality and spatial context

FCP PR #193 evaluated six semantically fixed ecological axes (BIO5, BIO12, BIO14, BIO15, annual solar radiation, mean water-vapour pressure) × four mutually exclusive photographed colour categories, adding a seventh soil available-water proxy in the soil-complete population. Global and four 90° longitude sectors, two source-quality regimes and both equal-species/equal-genus weightings yielded **1,280 distinct directional rows**. Every source-specific family used 999 matched genus bootstrap draws and approximate simultaneous max-t inference over that scope's axis×four hue set. The axes were defined after earlier exploratory FCP outcomes and do not constitute a historically prospective preregistration.

The global genus-equal photographic-white slope for annual precipitation was +2.255 percentage points / one original global BIO12 SD (95% scope-max-t interval +0.319 to +4.190); annual mean water-vapour was +3.952 (+1.199 to +6.704). Strict photo quality changed the BIO12 slope to +1.229 (−2.046 to +4.504), and genus-by-coarse-geographic-cell comparison to +0.808 (−3.556 to +5.172). East 90–180°E broad support showed +4.716 (+1.565 to +7.866), but its strict-photo and local spatial comparison intervals included zero. This is **source-observational geography/lineage-dependent covariation**, not natural selection.

Complete numerical source `results/fcp_macro_direction_20261013/direction_ledger.csv` (1,280 rows) and `results/fcp_macro_direction_20261013/decision_receipt.json`, reproduced from original data in hosted [Actions #38502660130](https://github.com/zuizui0223/fcp/actions/runs/38502660130), including 352 independently refitted raw-source regression coefficients.

**Supporting Information Fig. S13** will display all colour-by-environment directions and their source-specific bootstrap/simultaneous intervals, rather than selecting favourable colours.

## S15. Paired same-species comparison of spatial conditioning

FCP PR #194 fixes identical original colour and environmental row values across spatially conditioned models, separating eligibility from the changed geographic adjustment: A uses all eligible species in matched genera; B restricts to species in qualified genus-by-coarse-cell groups under smooth coordinate adjustment; C uses precisely those B species with genus-by-cell centring. The world B/C comparison retained **9,370 species in 845 genera** and East B/C **2,886 species in 280 genera**.

For genus-weighted photographic white–BIO12, B/C slopes were **+1.326/+0.808** (world) and **+1.102/−2.076** (East), in pp per original-global precipitation SD. East B−C difference **+3.178** had nominal genus-bootstrap 95% [+1.077,+5.351] but six-axis×four-colour max-t interval **[−0.258,+6.613]**. No conclusion of a universal climate effect or definitive causal weakening follows. The complete matched-stage and paired-difference ledger has **576 rows**, stored at `results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv`, with independent source recomputation of 288 source-row WLS coefficients (Actions #38502660130).

**Supporting Information Fig. S14** will show stage A/B/C support and identical-taxon conditional slopes, with their predeclared uncertainty and sensitivity HOLDs.

## S16. External descriptive phenotype and inference boundaries

A second colour measurement derives from Natural History Museum botanical prose, DOI `10.5519/p3dm31kc`, strict and consistent single-flower/petal/corolla white or blue/purple descriptions and exact binomials. After matching to the **same FCP photo-site geography**, 2,476 white/blue-purple taxa had 8,120 source photo-site environmental opportunities; the fixed East 90–180° comparison retained **58 shared genera / 670 named species**. The East-minus-other white–blue rainfall contrast was species-equal +2.07 (95% cluster CI −13.87 to +19.07) and genus-equal +3.56 (−15.34 to +22.60). Four-sector familywise inference did not confirm a positive direction.

The original NHM descriptive source had 59,369 records and some mixed or ambiguous descriptions; all excluded labels remain explicitly unresolved, not biological monomorphism. FCP PR #191 and [Actions #38405818980](https://github.com/zuizui0223/fcp/actions/runs/38405818980) record original-source and computational verification. The original photo-site geographic sampling remains shared with the primary FCP source, hence **phenotype measurement** is external but **distribution sampling** is not.

**Supporting Information Fig. S15** will document exact colour-label matching, mixed/unknown attrition and the nondirectional wide-interval external result.

### Supplemental table availability and non-pooling rule

**Supporting Information Table S1** distinguishes the four levels: 149,900 high-depth images and 1,499 named taxa; 100,543 sparse photographed taxon opportunities and 45,953 classified photo colours; 26/36 environmental complete cases and differing test-species supports; 2,476 external-colour-matched source taxa with original FCP photo geographic sites. The groups are not added to one independent biological sample size.

**Supporting Information Table S2** links each published numeric claim to the corresponding frozen JSON/CSV, original photo/environment SHA hashes, Actions replay and photographic or spatial sensitivity. The canonical reporting-only claim-to-source map is `docs/POLYMORPHISM_UNIFIED_EVIDENCE_AND_CLAIM_LEDGER_20261014.md`.

## S17. One-paper scientific stopping rules and submission archiving

The unified Full Paper can claim replicated geographic allocation of photographed flower colour and prospectively controlled phenotype-space alignment in the high-depth cohort. The environmental macro source may claim weak between-species predictive information and particular source-specific directions, with their simultaneous confidence limits, support dependence and alternative phenotype nonconfirmation. A global universal flower-colour boundary, climatic natural selection and genetic frequency-dependent morph maintenance are not identified and must not be presented as empirical conclusions.

The historical FCP high-depth release does **not** contain all post-October macro results. A single archive containing the unified manuscript, fully rendered Figure 6 and remaining figures, extended SI tables and the new source-hash verified analysis outputs must be frozen for submission. Author metadata and journal-facing file export remain human gates; no publication or submission is claimed.

## S18. Global photographically near-half white and the unverified ecological balance hypothesis

This is a **retrospective reporting-only integration of already-executed source tests**, not a new fitness experiment, independently confirmed dataset or prospective hypothesis. Historical PRs [#158](https://github.com/zuizui0223/fcp/pull/158), [#160](https://github.com/zuizui0223/fcp/pull/160), and [#168](https://github.com/zuizui0223/fcp/pull/168) control their original source inferences. Their exact original JSON receipts are copied unchanged to `results/np_unified_20261014/source_pr*.json`; a standalone arithmetic-only reconstruction and machine guard are `scripts/analysis/build_np_ecological_balance_bridge_20261014.py`, `results/np_unified_20261014/ecological_balance_reporting_reconciliation.json` and `tests/test_np_ecological_balance_20261014.py`. This report does NOT reprocess photographs or rerun historic 999-null simulations.

### Classified photo-white is near one half, but phenotypic missingness makes latent frequency unidentified

After excluding **26,905 duplicate photo IDs** across the archived 42,111-photo breadth and 85,337-photo species×cell measurements, the global source union has **100,543 unique original photos** from **42,111 nominal species**, with 45,953 colour-classifiable and 54,590 unclassified. Of classifiable original photos, **21,574/45,953 = 46.94797% photo-white**; this is photo-weighted, not species-equal, natural morph-frequency or individual-genotype frequency. With NO restrictions on the unknown classifications, white among all photos can range from **21,574/100,543 = 21.45749%** to **(21,574+54,590)/100,543 = 75.75266%**. The published independent one-photo-per-species breadth lineage has **8,597/18,457** white among classifiable species; because source imagery overlaps, agreement is not an independent biological replication.

| Original observed field | Source observation | Source-conditional inferential limit |
|---|---:|---|
| World classifiable photo-white | 21,574 / 45,953 = 46.95% | Photos, not equal-weight species |
| All original-photo unclassifiable hues | 54,590 / 100,543 = 54.30% | True flower hue unknown (MNAR possible) |
| Absolute-latitude 0–30° white | 9,247 / 19,642 = 47.08% | Only source classifiable geotagged photos |
| Absolute-latitude 30–60° white | 11,713 / 25,033 = 46.79% | Distinct photographic flora/effort |
| Absolute-latitude 60–90° white | 590 / 1,204 = 49.00% | Sparse high-latitude source |
| Latitude spread vs original within-species colour-preserving null | Observed 2.213 pp; null mean 2.687 pp; P = .214 | No unexpectedly constant global white share |
| Species-composition × local within-species deviation | Covariance +0.00018144; negative-compensation P = .997 | **No supported geographic cancellation maintaining a half-white proportion** |
| Original wettest–driest opportunity rain deciles | 53.253% minus 41.799% = +11.454 pp photo-white | Observational, selected on classifier success |
| MNAR tipping (scenario, not estimate) | If missing dry photos are 50% white, missing wet photos need only 41.264% to eliminate all-photo wet–dry contrast | Unobserved photo hue may explain it |
| Mid-rainfall local photo-white and nonwhite mixing | 5 evaluable tests, minimum unadjusted one-sided P = .387; all Holm P = 1.0 | No support; **zero species** with local-pair support in all three dry/mid/wet bins |

The original within-species 999 shuffle worlds preserve each source species' exact photographic white count and the true photo-site opportunity positions when assessing latitude. The geographical-decomposition test operates on 8,620 repeated-cell original photo species and 96 qualifying cells; the covariance is **positive**, opposite the compensatory negative expectation. A post hoc positive co-clustering permutation tail P=.004 weakens under source minimum-five-photo species and cell bootstrap, so it cannot be substituted as a new universal colour rule.

### Rainfall covariation is between species in the broad design, not an identified within-population selection gradient

The existing global genus-adjusted signed white/BIO12 association used 20,694 qualifying species and had slope +2.255 pp per source SD with family-maxT interval [+0.319,+4.190]. Strict image quality and genus×coarse-cell conditioning gave wider intervals including zero. Another original-source analysis of raw white photo labels found broad genus-adjusted wet–dry contrast +11.454 pp, whereas equal-size photo-source missing-colour sensitivity admits nonrandom hue failure sufficient to erase a *latent* contrast; these are different estimands.

The original unseen-genus matched soil-cohort comparison (PR #168) already evaluates **the same 20,786 soil-available species** under 26 predictors versus 36 including SoilGrids and finds an incremental held-out Brier gain **+0.000075729** (95% interval **−0.000516 to +0.000653**). It is unnecessary to launch a redundant generic 26-versus-36 comparison, though this is not a universal soil exclusion test.

A source-matched within-versus-between-sample study (PR #156) found slightly different rainfall Brier improvements on the same 14,776 original test photos; paired gain difference CI included zero. Original PR #159 had a positive broad rainfall relation for continuous photographic white-pixel fraction that did not survive strict ROI quality. Original PR #161 found a nominal rainfall environmental niche-separation difference among 422 species but nothing survived the 12-hypothesis Bonferroni family. PR #160's locally mixed photographs did **not** show a consistent intermediate-rainfall peak after colour-count-preserving tests. Different response grains, photo/source supports, inferential models and measurement error stop these from being pooled into one universal coefficient.

### Required biological experiment for any pigment benefit–cost conclusion

A true balance mechanism would need independently verified floral morph identity, actual sympatric reproducing populations, genotype or plant-individual tracking, pigment chemistry (including UV-absorbing colourless compounds), water-stress manipulations, pollinator/florivore interactions and **morph-specific relative male/female reproductive fitness**. The leading falsifiable prediction is a reproducible water-treatment × colour-morph interaction in relative fitness, with confounding by taxonomy, flowering season and colour-photo quality controlled. Neither near-50% photographic composition nor an aligned nine-colour axis alone requires that interaction; **a universal 50:50 world frequency is not predicted by the hypothesis**.

The prior FCP NEXT-2000 project (PR #133) has only frozen **2,000 new source-species identities and 730 unopened reserve identities**. It has no new biological colour measurement, and its preallocated taxa must not be reselected in response to this retrospectively constructed moisture-hypothesis narrative.

**Claim gate:** broad moisture–photo-white association = observational with sensitive transport; prospective recurrent white–chromatic axis = supported phenotype geometry conditional on photographic construction and exposure; geographic negative compensation / extraordinary global half = unsupported; measured pigment benefit–cost equilibrium = NOT TESTED.
