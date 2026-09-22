# Flower-colour polymorphism current claim ledger — 2026-09-18

Date: 2026-09-18 JST

This ledger supersedes `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260916.md` for manuscript claim-status purposes.

The authoritative third-cohort biological result is:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

The authoritative third-cohort measurement/support receipt is:

- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`

If this prose conflicts with those machine-readable files, the machine-readable files control.

## 0. Programme lineage — RGFCA versus the current paper

RGFCA means **Repeated Global Flower-Colour Atlas**. It is the upstream global sampling, measurement and species-specific spatial-analysis framework from which the current discovery/reserve resource was inherited; it is not the present paper's final biological claim.

RGFCA originally asked whether independent species repeatedly place strong within-species flower-colour discontinuities in the same broad geographic regions. That shared-geography estimand did not become the positive biological spine of the current paper. The upstream frozen record includes primary recurrent-field G1 p = **0.070** and species-disjoint commonness p = **0.856**, while later sharedness qualification identified a support/identifiability limitation.

The current paper therefore changes the level at which cross-species generality is sought:

- **amount:** reproducible species-level D;
- **phenotype geometry:** prospectively confirmed excess alignment with a frozen white-versus-nonwhite axis relative to a coarse-state-preserving null;
- **spatial realization:** stronger species-specific geographic organization at higher D;
- **broad explanation:** not reducible to sampled geographic span or broad tree-wide phylogenetic conservation.

The supported conceptual interpretation is that generality is stronger in **phenotype space than in shared geographic location** under the tested designs. This does not imply absence of flower-colour biogeography or identify the ecological mechanism maintaining polymorphism.

Canonical programme interpretation:

- `docs/RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md`

Canonical 42,111-frame provenance:

- `docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md`

## 1. Paper mainline

The paper now has three positive/structural contributions and two bounded alternative-explanation results:

1. **H1 — measurement validity.** The continuous four-state species-level flower-colour polymorphism score
   `D = 1 - sum_k p_k^2` is reproducible across observer-disjoint photo sets under the first-frozen repeated-partition validation rule. A later deliberately stricter deterministic split missed a prespecified rho = 0.80 floor, so near-perfect or split-invariant reliability is not claimed.
2. **H2 — geometry.** Existing discovery/reserve cohorts localized recurrent within-species colour variation to a white-versus-nonwhite axis after a construction-preserving audit. In the pre-frozen species-disjoint third cohort, observed alignment exceeded the already white-aligned coarse-state-preserving structured null. The confirmatory quantity is therefore the increment above that construction baseline, not the entire white-axis signal relative to isotropy.
3. **Spatial organization — replicated structural correlate.** Species with greater D also tend to show stronger within-species geographic colour organization. This association replicated in the species-disjoint reserve and survived sampled-span plus clear technical-failure adjustment, a matched flower-minus-background contrast, and uniform ambiguity-endpoint stress tests. It is a structural correlate, not a causal mechanism.
4. **H3a — phylogeny.** Broad tree-wide phylogenetic signal in D is not supported in the species-disjoint reserve cohort across any of the three frozen tree-placement scenarios.
5. **H3b — sampled geographic span.** The discovery association between sampled photographic span and D does not replicate in reserve.

The paper is therefore about **measurement + recurrent geometry + spatial organization**, with H3 used to reject two simple broad explanations. It is not a predictor-hunting paper, a global prevalence paper, or a shared-boundary paper.

## 2. Global frame and sampling boundary

The global flower-colour sampling frame contains 42,111 species.

The high-depth cohorts are hypothesis-specific validation/measurement cohorts, not a probability sample from those 42,111 species. Their species counts must not be interpreted as a global estimate of flower-colour polymorphism prevalence.

The 42,111-species frame and later opportunity subsets define where high-depth sampling could be attempted. Measurement-gate failure or insufficient classifiable observations are missing/underidentified states, not biological monomorphism.

## 3. H1 — positive measurement-validity claim with stress-test caveat

Canonical source:

- `docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md`

Primary first-frozen repeated observer-disjoint validation in reserve:

- 200 partitions;
- median paired species = 329;
- median split Spearman rho = **0.789102993**;
- 5th percentile rho = **0.765164994**;
- median Lin CCC = **0.854830674**;
- median Spearman-Brown reliability = **0.882121371**.

Frozen primary verdict:

`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`

Later deterministic stress test:

- reliability-eligible reserve species = 363;
- observer leakage = 0;
- Spearman rho = **0.792727693**;
- bootstrap 95% interval = **0.741867914–0.832365877**;
- Lin CCC = **0.847429024**.

The stricter test required rho >= 0.80 and therefore returned:

`H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`

### Allowed H1 claim

> The continuous four-state polymorphism score D was reproducible across repeated observer-disjoint partitions under the first-frozen reserve validation rule (median split rho 0.789; 5th percentile 0.765; median CCC 0.855). A later deliberately stricter deterministic split yielded rho 0.793 and missed a prespecified 0.80 floor, so near-perfect or split-invariant reliability is not claimed.

## 4. H2 — legacy discovery/reserve localization

Canonical pre-prospective sources:

- `docs/POLYMORPHISM_H2_EVIDENCE_LEDGER_20260912.md`
- `docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`
- `results/polymorphism_white_axis_targeted_test_20260912/result.json`

The existing cohorts first established recurrent label-free directional concentration and then localized the construction-controlled signal to the fixed contrast

`q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`.

For

`W = mean_i (u_i dot q_white)^2`

under the frozen coarse-state-preserving structured null:

Primary 0.10 tier:

- discovery: N = 152, W = **0.514625**, null median = 0.430808, p = **0.001**;
- reserve: N = 129, W = **0.514586**, null median = 0.466546, p = **0.001**.

Strict 0.20 tier:

- discovery: N = 75, W = **0.542355**, null median = 0.443626, p = **0.001**;
- reserve: N = 65, W = **0.510517**, null median = 0.469943, p = **0.008**.

After projection of q_white, residual directional concentration no longer exceeds the construction-preserving null. The non-white-only diagnostic likewise does not support a recurrent hue direction.

### Chronology boundary

The specific white/nonwhite axis was isolated after the broad H2 geometry had already been opened in the original discovery/reserve cohorts. Those cohorts therefore supply discovery/audit/target-localization evidence, not untouched prospective confirmation of the named axis.

## 5. P500 — closed prospective measurement transport, no durable H2 verdict

Canonical status:

- `docs/P500_PROSPECTIVE_TERMINAL_POSTMORTEM_20260916.md`
- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260916.md`

P500 carried a separately frozen 499-species / 49,900-photo expansion through all 256 location-blind measurement partitions and passed its predeclared support gate with 373 evaluable species.

Its H2 calculation reached an in-memory validated terminal object but failed during post-calculation JSON serialization before a durable H2 result was written.

Therefore:

- P500 supports **measurement-pipeline transport**;
- P500 supplies **no durable confirmatory evidence for or against H2**;
- P500 must not be replayed or retrospectively reclassified as an untouched prospective H2 test.

## 6. H2 — untouched prospective third-cohort confirmation

Canonical sources:

- `docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`
- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`
- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

Frozen verdict:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

Status:

`untouched_prospective_test_of_previously_frozen_axis`

### Selection and measurement

- outcome-blind candidate frame = 3,230 species;
- frozen selected cohort = 500 species;
- fresh metadata yielded 499 authorized species with 100 rows each;
- terminal rows = **49,900**;
- unique measurement IDs = **49,900**;
- duplicate measurement IDs = 0;
- terminal partitions = **256 / 256**;
- classifiable rows = **25,788**;
- nonclassifiable rows = 24,112;
- measurement-evaluable species (>=40 classifiable) = **377**;
- minimum required evaluable species = 250;
- support decision = **PASS**;
- replacement species = 0;
- replacement rows = 0;
- persisted image pixels = false.

The support gate passed before H2 was opened.

### Primary 0.10 tier

- vector species = **158**;
- observed W = **0.5172457461**;
- structured-null repetitions = 999;
- null median = **0.4571428150**;
- null 95% interval = **0.4358491120–0.4752987776**;
- observed / null median = **1.131475174**;
- upper-tail p = **0.001**;
- decision = **PASS**.

### Strict 0.20 tier

- vector species = **86**;
- observed W = **0.5329282123**;
- structured-null repetitions = 999;
- null median = **0.4593196659**;
- null 95% interval = **0.4328679570–0.4867224043**;
- observed / null median = **1.160255595**;
- upper-tail p = **0.001**;
- decision = **PASS**.

### Allowed H2 main claim

> In a pre-frozen species-disjoint third cohort drawn from the same iNaturalist opportunity universe, continuous within-species colour displacement showed excess alignment with the previously specified white-versus-nonwhite axis relative to the unchanged coarse-state-preserving structured null. At the primary tier, W = 0.517 versus null median 0.457 (158 species; p = 0.001); the strict tier also passed.

### Scope boundary

This is a prospective species-disjoint confirmation/transport test **within the same iNaturalist source and opportunity universe**. It is not an independent-source replication.

## 6b. Post-confirmatory H2 validity and direct highlight control

Canonical sources:

- `docs/POLYMORPHISM_H2_POSTHOC_VALIDITY_DIAGNOSTICS_20260922.md`
- `results/polymorphism_h2_posthoc_validity_diagnostics_20260922/result.json`
- `docs/POLYMORPHISM_H2_THIRD_COHORT_HIGHLIGHT_VALIDITY_PROTOCOL_20260922.md`
- `results/polymorphism_h2_third_cohort_highlight_validity_20260922/result.json`
- `docs/POLYMORPHISM_H2_THIRD_COHORT_HIGHLIGHT_DECISION_ADJUDICATION_20260923.md`

These analyses were performed after the prospective H2 result was opened and **do not alter the frozen prospective H2 verdict or estimand**.

### Construction and gate diagnostics

- isotropic expectation = **0.125**;
- frozen structured-null median = **0.457143**;
- observed W = **0.517246**;
- 137/158 primary vector species have white as one of the two leading coarse morphs;
- mean per-species white-axis contribution = **0.5622** for those species versus **0.2239** for the 21 without white among their two leading coarse morphs;
- a 299-replicate null that re-applies the continuous minor-cluster gate has median **0.447495**, 95% interval **0.428986–0.469580**, plus-one p = **0.00333**.

### Direct one-shot digital-highlight control

The third-cohort direct validity control reacquired the complete frozen 49,900-row denominator before opening biological outcomes.

- reacquired rows = **49,900 / 49,900**;
- acquisition failures = **0**;
- source-byte drift = **0**;
- highlight metrics available = **44,098** rows;
- response-blind high-clip threshold = near_clip_fraction > **0.4933450501**;
- high-clip rows = **2,205**;
- coupling model = **23,320 rows / 461 species**;
- OR per one within-species SD near-clip increase = **1.444493**;
- 95% CI = **1.389332–1.501845**.

The complete OR interval lies above the predeclared negligible-coupling upper boundary of **1.25**. Therefore the measured coarse white state is directly demonstrated to be exposure-coupled; it must not be described as artifact-free.

### High-clip-exclusion H2 sensitivity

After removing the response-blind high-clip set:

- primary vector species = **142**;
- retention = **142 / 158 = 0.898734 (89.9%)**;
- W = **0.503428**;
- structured-null median = **0.453801**;
- structured-null p = **0.001**;
- H2 support = **retained**.

The sensitivity retained one fewer vector than the >=143 vectors required for the frozen 0.90 retention threshold.

### Decision-precedence boundary

The frozen prose protocol contains a FLAGGED clause when the complete coupling-model CI lies outside OR [0.80, 1.25]. The frozen executable, however, checks <0.90 H2-vector retention before reaching that coupling-CI branch. Both were frozen before reacquisition.

The generated terminal machine state therefore remains:

`INDETERMINATE`

and is **not retrospectively recoded** after seeing the outcome.

The scientific interpretation is nevertheless fixed:

> Near-clipping is substantially coupled to frozen white classification within species. Removing the response-blind high-clip set does not remove H2 support, but the measured coarse white state cannot be treated as cleared of image-exposure effects. The H2 claim is therefore retained only as excess alignment conditional on the measured colour-state construction, not as proof of a purely biological white-versus-nonwhite axis.

## 7. Replicated spatial organization of D — positive structural clue

Canonical reporting receipt:

- `results/polymorphism_spatial_organization_clue_20260918/result.json`

This receipt is reporting-only and reproduces previously frozen results from PR #32 / `feat/polymorphism-paper-v0-1-post-step9`; it performs no new biological analysis.

Raw discovery and reserve associations:

- discovery: rho(D, within-species spatial organization) = **0.0892133**, p = **0.034**;
- reserve: rho = **0.1016008**, p = **0.025**.

After controlling for sampled geographic span and the rate of clear ROI/flip technical failures:

- discovery primary partial rho = **0.1266367**, geometry-preserving null p = **0.007**;
- reserve primary partial rho = **0.0992877**, p = **0.025**;
- reserve matched flower-minus-background partial rho = **0.1162411**, p = **0.010**.

The reserve result also survives exact uniform ambiguity-endpoint stress tests under the four-state completion model:

- primary D_min4: rho = **0.0970781**, p = **0.029**;
- primary D_max4: rho = **0.1252858**, p = **0.008**;
- flower-minus-background D_min4: rho = **0.1162986**, p = **0.009**;
- flower-minus-background D_max4: rho = **0.1327852**, p = **0.006**.

### Allowed spatial-organization claim

> Species with greater four-state flower-colour diversity tend to show stronger within-species geographic colour organization across the discovery and species-disjoint reserve high-depth cohorts. The reserve association persists after sampled-span and clear technical-failure adjustment, a matched flower-minus-background contrast and uniform ambiguity-endpoint stress tests.

### Scope boundary

This does **not** show that geographic organization causes high D, nor does it identify climate, pollinators, habitat, gene flow, drift, mating system or any other maintenance mechanism. It establishes a replicated spatial correlate that narrows the mechanistic interpretation of the between-species differences in D.

## 8. H3a — broad phylogenetic signal not supported

Canonical source:

- `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`

Reserve retains 341 tips on each frozen S1-S3 placement scenario.

- S1: K = **0.0710190**, p = **0.2716**; lambda = 0.04857, p(lambda=0) = 0.1657; opportunity-adjusted p = 0.3500.
- S2: K = **0.0601476**, p = **0.4134**; lambda = 0.05113, p(lambda=0) = 0.1555; opportunity-adjusted p = 0.2464.
- S3: K = **0.0707577**, p = **0.2674**; lambda = 0.04746, p(lambda=0) = 0.1702; opportunity-adjusted p = 0.3510.

Frozen verdict:

`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`

This does not imply that phylogeny is biologically irrelevant. It only closes the tested broad tree-wide signal claim under the frozen design.

## 9. H3b — sampled-span association does not replicate

Canonical source:

- `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Discovery calibration:

- n = 369;
- Spearman rho(D, log1p sampled span) = **0.1798786**;
- p = **0.00089996**.

Fresh reserve replication:

- n = 363;
- rho = **-0.0025855**;
- p = **0.9586021**;
- adjusted partial-rank rho = 0.0055187, p = 0.9162042;
- S1-S3 rank-PGLS beta = -0.0055931, p = 0.9136754.

Frozen verdict:

`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`

The frozen predictor is sampled photographic span, not true biological range size.

## 10. Current manuscript claim ceiling

The strongest defensible paper-level statement is now:

> Species-level flower-colour diversity can be measured reproducibly from high-depth citizen-science photographs under observer-disjoint validation, although reliability is not split-invariant or near perfect. In a pre-frozen species-disjoint third cohort, continuous colour displacement showed excess alignment with the frozen white-versus-nonwhite axis relative to a coarse-state-preserving structured null. The increment is robust to reapplying the continuous-cluster gate within null worlds, but digital exposure/background-context confounding of the measured coarse white state remains unresolved. Across the original high-depth cohorts, species with greater D also show stronger within-species geographic colour organization. Broad tree-wide phylogenetic conservation is not detected, and the discovery sampled-span association collapses in reserve.

This upgrades the former post-audit H2 claim to a prospective confirmation **of excess alignment with the fixed axis relative to the frozen structured null**, while preserving both the discovery chronology and the unresolved white-state measurement-validity boundary.

## 11. Working title authorization

The prospective third-cohort result permits a title centered on achromatic–chromatic alignment only if it makes the construction-controlled nature of the inference clear.

Preferred working title:

**Within-species flower-colour variation shows achromatic–chromatic alignment beyond coarse colour-state composition**

A safer alternative emphasizing measurement:

**Within-species flower-colour variation is reproducible and shows excess achromatic–chromatic alignment**

## 12. Hard nonclaims

The current evidence does not establish:

- global prevalence of flower-colour polymorphism among the 42,111-species frame;
- an independent-source replication of H2;
- pigment chemistry or pigment-loss/gain mechanism;
- direction of evolutionary transitions between white and non-white states;
- pollinator, climate, habitat, or other adaptive causation;
- a recurrent non-white hue axis;
- absence of all phylogenetic or taxonomic structure;
- irrelevance of true geographic range size;
- near-perfect or split-invariant H1 reliability;
- that the coarse white state is free of digital exposure, scene-background or ROI-contamination effects.

No post-confirmatory analysis may change q_white, W, construction/admissibility gates, the structured null, cohort definition, or decision rule and still be described as the untouched prospective third-cohort test.
