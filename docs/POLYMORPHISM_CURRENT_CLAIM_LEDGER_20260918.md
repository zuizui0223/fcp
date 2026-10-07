# Flower-colour polymorphism current claim ledger — 2026-09-18

Date: 2026-09-18 JST

This ledger supersedes `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260916.md` for manuscript claim-status purposes.

The authoritative prospective-confirmation biological result is:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

The authoritative prospective-confirmation measurement/support receipt is:

- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`

If this prose conflicts with those machine-readable files, the machine-readable files control.

## 0. Reader-facing data architecture and inferential necessity

The current paper uses one explicit inferential sequence.

1. **Outcome-blind opportunity frame.** Metadata-only iNaturalist discovery defines a 42,111-species candidate universe before flower colour is examined; 4,730 species can support >=100 retained photographs after observer capping. This stage is required to keep species entry independent of the focal phenotype and is not a prevalence denominator.
2. **Discovery + species-disjoint validation.** A fixed high-depth resource contains 500 discovery and 500 validation species with 100 photographs per species. This stage is required both to estimate within-species distributions rather than modal colours and to separate pattern discovery from validation on different species.
3. **Prospective confirmation.** Because the white-versus-nonwhite axis was identified after the original geometry was opened, that axis, W statistic, admissibility thresholds and structured null were frozen before newly sampled species and photographs were opened. A 3,230-species unused candidate frame yielded 499 species × 100 new photographs after outcome-blind selection. This fresh cohort is required for the named axis to have prospective confirmatory status.
4. **Post-confirmatory annotations.** Highlight metrics, WorldClim variables, phylogenetic placements and sampled span test measurement coupling and simple alternative explanations. These analyses constrain interpretation but do not create or upgrade the prospective H2 confirmation.

The paper's public story is therefore **measurement validity → geographic allocation of ITV → bounded spatial-process decomposition → complementary prospective phenotype-space confirmation**.

## 1. Paper mainline

The paper now has four positive/structural contributions and several explicit limits.

1. **Measurement validity.** Species-wide sampled four-state diversity D is reproducible across repeated observer-disjoint partitions and shows strong fresh-image same-system transport. D is a range-wide sample property, not local polymorphism magnitude.

2. **Distributed polymorphism — central ecological result.** At the fixed post hoc 50-km scale, local colour-state diversity is lower than expected from each species' exact overall colour composition in discovery, validation and the third cohort (mean depletion 0.02053, 0.01867 and 0.01468; matched p=0.005 in each). The same direction is supported at 25, 100 and 250 km. The result survives exclusion of same-observer pairs, removal of white records and continuous nine-colour reanalysis. It establishes geographic partitioning of ITV, not adaptation.

3. **Spatial-process decomposition.** Continuous nine-colour turnover contains both positive IBD-like and BIO5-associated IBE-like components in all three cohorts, with IBD-like effects consistently larger. The BIO5 component is technically bounded because validation flower-minus-background and strict same-observer sensitivities are unsupported. Geography therefore remains the dominant robust spatial correlate; environmental sorting is a candidate, not a causal conclusion.

4. **Prospectively confirmed phenotype geometry.** The original 1,000 species localize recurrent displacement to an achromatic–chromatic axis. In the pre-frozen species-disjoint prospective cohort, W=0.517 exceeds the already white-aligned structured-null median 0.457 (p=0.001). This confirms excess alignment beyond the coarse-state construction baseline, not a universal pigment or adaptive mechanism.

Important negative/bounded results:
- higher D does not show a scale-invariant relationship with local depletion;
- broad multivariate environmental heterogeneity does not replicate as a general explanation for D;
- solar heterogeneity selected in the 500+500 screen fails third-cohort transport;
- a universal signed warm-to-white BIO5 rule fails cross-cohort transport and weakens under observer controls;
- broad tree-wide phylogenetic signal in D is unsupported under the tested validation design;
- the discovery sampled-span association collapses in validation;
- the coarse white classifier is exposure-coupled.

The paper is therefore about **where flower-colour ITV resides, how its spatial structure decomposes, and which phenotype-space directions recur**. It is not a local-adaptation proof, a predictor-hunting paper, a global prevalence paper or a universal temperature-rule paper.

## 2. Global frame and sampling boundary

The global flower-colour sampling frame contains 42,111 species.

The high-depth cohorts are hypothesis-specific validation/measurement cohorts, not a probability sample from those 42,111 species. Their species counts must not be interpreted as a global estimate of flower-colour polymorphism prevalence.

The 42,111-species frame and later opportunity subsets define where high-depth sampling could be attempted. Measurement-gate failure or insufficient classifiable observations are missing/underidentified states, not biological monomorphism.

## 3. H1 — positive measurement-validity claim with stress-test caveat

Canonical source:

- `docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md`
- `results/polymorphism_fresh_D_transport_20260925/result.json`

Primary first-frozen repeated observer-disjoint validation in the validation cohort:

- 200 partitions;
- median paired species = 329;
- median split Spearman rho = **0.789102993**;
- 5th percentile rho = **0.765164994**;
- median Lin CCC = **0.854830674**;
- median Spearman-Brown reliability = **0.882121371**.

Frozen primary verdict:

`H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED`

Later deterministic stress test:

- reliability-eligible validation species = 363;
- observer leakage = 0;
- Spearman rho = **0.792727693**;
- bootstrap 95% interval = **0.741867914–0.832365877**;
- Lin CCC = **0.847429024**.

The stricter test required rho >= 0.80 and therefore returned:

`H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED`

### Allowed H1 claim

> The continuous four-state polymorphism score D was reproducible across repeated observer-disjoint partitions under the first-frozen validation-cohort rule (median split rho 0.789; 5th percentile 0.765; median CCC 0.855). In a later fresh-image execution within the same iNaturalist/FCP measurement system, 136 overlapping species showed Spearman rho 0.968, Lin CCC 0.972, calibration slope 0.969 and median absolute D change 0.0148. This strengthens same-system transport of D but is not independent-source replication; near-perfect or split-invariant reliability is not claimed.

## 4. H2 — discovery/validation target localization

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
- validation: N = 129, W = **0.514586**, null median = 0.466546, p = **0.001**.

Strict 0.20 tier:

- discovery: N = 75, W = **0.542355**, null median = 0.443626, p = **0.001**;
- validation: N = 65, W = **0.510517**, null median = 0.469943, p = **0.008**.

After projection of q_white, residual directional concentration no longer exceeds the construction-preserving null. The non-white-only diagnostic likewise does not support a recurrent hue direction.

### Chronology boundary

The specific white/nonwhite axis was isolated after the broad H2 geometry had already been opened in the original discovery/validation cohorts. Those cohorts therefore supply discovery/audit/target-localization evidence, not untouched prospective confirmation of the named axis.

## 5. H2 — prospective confirmation in newly sampled species

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

> In a pre-frozen species-disjoint prospective confirmation cohort drawn from the same iNaturalist opportunity universe, continuous within-species colour displacement showed excess alignment with the previously specified white-versus-nonwhite axis relative to the unchanged coarse-state-preserving structured null. At the primary tier, W = 0.517 versus null median 0.457 (158 species; p = 0.001); the strict tier also passed.

### Scope boundary

This is a prospective species-disjoint confirmation/transport test **within the same iNaturalist source and opportunity universe**. It is not an independent-source replication.

## 5b. Post-confirmatory H2 validity and direct highlight control

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

The prospective-cohort direct validity control reacquired the complete frozen 49,900-row denominator before opening biological outcomes.

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

## 5c. Post-confirmatory environmental clue — BIO5 passes in the prospective cohort but does not transport

This analysis was specified after the prospective H2 result had been opened and therefore cannot alter H2 or be described as part of its untouched confirmation.

In the prospective confirmation cohort, after the response-blind high-clip exclusion:
- eligible species = **281**;
- median within-species white-minus-nonwhite BIO5 contrast = **+0.0690113 SD**;
- species-level Wilcoxon p = **0.0118162**;
- Holm-adjusted p across BIO5, BIO14 and mean solar radiation = **0.0354487**;
- species-stratified conditional OR per one within-species SD BIO5 = **1.073475**;
- 95% CI = **1.029395–1.119441**;
- p = **0.0009188**;
- frozen BIO5 mechanism gate = **PASS**.

BIO14 and mean solar radiation did not pass their frozen gates.

Post hoc observer sensitivities weaken the prospective BIO5 result:
- observer-paired design: 106 species / 144 paired species-observer strata, median BIO5 contrast = **0.000 SD**, Wilcoxon p = **0.484962**, conditional OR = **0.786636** (95% CI **0.543940–1.137618**), p = **0.202310**;
- observer-balanced design: 352 species, median BIO5 contrast = **+0.054098 SD**, Wilcoxon p = **0.074964**, sign-test p = **0.048443**.

These were opened after the primary environmental result. They do not overwrite the frozen within-cohort gate, but they show observer-conditioning sensitivity.

A later fixed BIO5-only transport test in the original species-disjoint cohorts failed the joint replication rule.

**Discovery**
- 271 eligible species;
- median contrast = **-0.006777 SD**;
- species-level p = **0.743252**;
- conditional OR = **1.033114**, p = **0.130829**.

**Validation**
- 260 eligible species;
- median contrast = **+0.066165 SD**;
- species-level p = **0.054067**;
- conditional OR = **1.048240**, p = **0.031950**.

Frozen transport decision:

**Not supported under the prespecified cross-cohort rule.**

Allowed interpretation:

> The frozen primary prospective analysis supports a within-cohort association between white states and warmer BIO5 environments, but post hoc observer controls weaken that association and the result does not support a common cross-cohort BIO5 rule.

This is an observer-sensitive environmental sorting association, not causal heat selection, and it remains within the same iNaturalist/FCP source and measurement system.

### Evidence-placement boundary

- **BIO5 is an empirical result of the present study** and remains in the main Results, with its observer sensitivity and failed cross-cohort transport reported alongside it.
- **F3h1/Myb1a expression and HPLC evidence are published results from Casimiro-Soriguer et al. (2016)** and belong in Discussion as single-species mechanistic context, not as a present-study molecular result.
- **The PAL/WAL frequency comparison is the present study's reanalysis of published tables** and is reported in Supporting Information Section S8; main text may cite its bounded interpretation but not elevate it to a primary result.

## 6. Replicated spatial organization of D — positive structural clue

Canonical reporting receipt:

- `results/polymorphism_spatial_organization_clue_20260918/result.json`

This receipt is reporting-only and reproduces previously frozen results; it performs no new biological analysis.

Raw discovery and validation associations:

- discovery: rho(D, within-species spatial organization) = **0.0892133**, p = **0.034**;
- validation: rho = **0.1016008**, p = **0.025**.

After controlling for sampled geographic span and the rate of clear ROI/flip technical failures:

- discovery primary partial rho = **0.1266367**, geometry-preserving null p = **0.007**;
- validation primary partial rho = **0.0992877**, p = **0.025**;
- validation matched flower-minus-background partial rho = **0.1162411**, p = **0.010**.

The validation result also survives exact uniform ambiguity-endpoint stress tests under the four-state completion model:

- primary D_min4: rho = **0.0970781**, p = **0.029**;
- primary D_max4: rho = **0.1252858**, p = **0.008**;
- flower-minus-background D_min4: rho = **0.1162986**, p = **0.009**;
- flower-minus-background D_max4: rho = **0.1327852**, p = **0.006**.

### Allowed spatial-organization claim

> Species with greater **species-wide sampled four-state colour diversity** tend to show stronger within-species geographic colour organization across the discovery and species-disjoint validation high-depth cohorts. D summarizes the retained range-wide sample and can combine local coexistence with differentiation among sampled locations. The validation association persists after sampled-span and clear technical-failure adjustment, a matched flower-minus-background contrast and uniform ambiguity-endpoint stress tests.

### Scope boundary

This does **not** show that geographic organization causes high D, nor does it identify climate, pollinators, habitat, gene flow, drift, mating system or any other maintenance mechanism. The frozen spatial null tests organization relative to an unstructured assignment; it does not test equal-strength true spatial processes across different D values. It establishes a replicated spatial correlate that narrows the mechanistic interpretation of the between-species differences in D.

## 7. H3a — broad phylogenetic signal not supported

Canonical source:

- `results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`

Validation retains 341 tips on each frozen S1-S3 placement scenario.

- S1: K = **0.0710190**, p = **0.2716**; lambda = 0.04857, p(lambda=0) = 0.1657; opportunity-adjusted p = 0.3500.
- S2: K = **0.0601476**, p = **0.4134**; lambda = 0.05113, p(lambda=0) = 0.1555; opportunity-adjusted p = 0.2464.
- S3: K = **0.0707577**, p = **0.2674**; lambda = 0.04746, p(lambda=0) = 0.1702; opportunity-adjusted p = 0.3510.

Frozen verdict:

`H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`

This does not imply that phylogeny is biologically irrelevant. It only closes the tested broad tree-wide signal claim under the frozen design.

## 8. H3b — sampled-span association does not replicate

Canonical source:

- `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`

Discovery calibration:

- n = 369;
- Spearman rho(D, log1p sampled span) = **0.1798786**;
- p = **0.00089996**.

Fresh validation replication:

- n = 363;
- rho = **-0.0025855**;
- p = **0.9586021**;
- adjusted partial-rank rho = 0.0055187, p = 0.9162042;
- S1-S3 rank-PGLS beta = -0.0055931, p = 0.9136754.

Frozen verdict:

`H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`

The frozen predictor is sampled photographic span, not true biological range size.

## 9. Current manuscript claim ceiling

The strongest defensible paper-level statement is now:

> Across three species-disjoint high-depth cohorts, species-wide flower-colour variation is geographically partitioned: nearby conspecific observations contain less colour diversity than expected from each species' fixed overall colour composition. This local depletion persists after excluding same-observer pairs, removing the exposure-sensitive white class and using continuous nine-colour distances. Continuous colour turnover contains both a larger IBD-like component and a smaller BIO5-associated residual component, but the latter fails stricter background/observer controls and is not evidence of local adaptation. Independently, an untouched prospective species-disjoint test confirms excess achromatic–chromatic displacement alignment relative to a coarse-state-preserving structured null.

The paper therefore separates **evidential status** from **biological centrality**: H2 remains the strongest prospectively protected named-axis result, while distributed polymorphism is the broadest replicated ecological pattern.

## 9a. Allowed ecological synthesis — intraspecific variation as a distributional comparative trait

This is an inte## 9a. Allowed ecological synthesis — ITV as a geographically distributed species property

This is an interpretation of the reported results, not an additional statistical test.

The paper distinguishes four properties of the same within-species distribution:

- **overall sampled diversity** — species-wide four-state D;
- **spatial allocation** — local depletion relative to the fixed species-wide composition;
- **spatial process** — IBD-like and bounded IBE-like turnover;
- **phenotype-space direction** — continuous displacement geometry.

The allowed general synthesis is:

> A species-wide ITV value can conceal a distributed polymorphism: flower-colour variants are not merely pooled within species but are repeatedly partitioned among geographic localities.

The allowed evolutionary interpretation is:

> Spatially varying selection is one plausible generator of this distributed polymorphism, but drift, restricted dispersal, colonisation history and phenotypic plasticity remain viable alternatives. The current photographs identify a general spatial pattern that requires an evolutionary explanation; they do not prove local adaptation.

The BIO5 residual may be described only as **phenotypic IBE-like** and technically bounded. It must not be called genetic IBE, causal temperature selection or local adaptation.

The achromatic–chromatic result remains a complementary statement about phenotype-space accessibility. Published pigment-pathway evidence can motivate developmental accessibility but cannot establish a shared molecular cause across the sampled species.

This synthesis does **not** claim:
- the first distributional treatment of ITV;
- that all range-wide flower-colour diversity is among populations;
- a scale-invariant D-to-partitioning law;
- a universal ecological driver;
- a universal achromatic–chromatic law across angiosperms;
- adaptive differentiation without fitness/genetic evidence.

## 10. Current title

**Species-wide flower-colour variation is geographically partitioned across plant species**

This wording centers the most replicated ecological pattern without upgrading it to local adaptation. The prospective achromatic–chromatic confirmation remains explicit in the Summary, Results and Discussion as a complementary phenotype-space result.

## 11. Hard nonclaims

The current evidence does not establish:

- global prevalence of flower-colour polymorphism among the 42,111-species frame;
- an independent-source replication of H2;
- pigment chemistry or pigment-loss/gain mechanism;
- direction of evolutionary transitions between white and non-white states;
- pollinator, climate, habitat, or other adaptive causation;
- a universal or replicated BIO5–white effect across cohorts;
- whether geographic realization is shared or species-specific across species;
- a recurrent non-white hue axis;
- absence of all phylogenetic or taxonomic structure;
- irrelevance of true geographic range size;
- near-perfect or split-invariant H1 reliability;
- that the coarse white state is free of digital exposure, scene-background or ROI-contamination effects.
- that geographic partitioning is caused by local adaptation rather than drift, dispersal limitation, history or plasticity;
- genetic isolation by environment;
- that the BIO5 residual is flower-specific or observer-independent under current controls;
- that higher-D species are universally more strongly partitioned at every spatial scale;

No post-confirmatory analysis may change q_white, W, construction/admissibility gates, the structured null, cohort definition, or decision rule and still be described as the untouched prospective-confirmation test.
