# FCP — Global flower-colour polymorphism

This repository's active mainline asks two linked questions: **can within-species flower-colour polymorphism be measured reproducibly as a species phenotype, and is that variation geometrically constrained in colour space?**

The former shared-boundary / Repeated Global Flower-Colour Atlas (RGFCA) programme is retained for provenance, but it is no longer the positive biological mainline. Its non-support and identifiability limits remain valid historical results.

## Current paper

Working title:

**A recurrent achromatic–chromatic axis structures within-species flower-colour polymorphism across plant species**

The current paper separates three inferential layers:

1. **measurement validity** — whether a continuous species-level polymorphism score survives observer separation;
2. **geometry** — whether within-species colour variation repeatedly follows a fixed direction in colour space;
3. **bounded explanation** — whether broad phylogenetic signal or sampled photographic span explain that phenotype.

The decisive upgrade is now complete: a previously frozen white-versus-nonwhite axis was tested prospectively in a pre-frozen species-disjoint third cohort drawn from the same iNaturalist opportunity universe and was supported at both the primary and strict tiers.

## Start here

- **Active manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT.md`](docs/POLYMORPHISM_MANUSCRIPT.md)
- **Current claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
- **Current figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
- **Paper architecture:** [`docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260918.md`](docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260918.md)
- **Supporting Information:** [`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`](docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md)
- **Literature audit:** [`docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md`](docs/POLYMORPHISM_LITERATURE_AUDIT_20260918.md)
- **Canonical figures + SHA manifest:** [`docs/figures/polymorphism_20260918/`](docs/figures/polymorphism_20260918/)
- **H1 reconciliation:** [`docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md`](docs/POLYMORPHISM_H1_EVIDENCE_LEDGER_20260914.md)
- **Frozen white-axis target:** [`docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md`](docs/POLYMORPHISM_H2_WHITE_AXIS_TARGET_FREEZE_20260912.md)
- **Third-cohort protocol:** [`docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md`](docs/POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md)
- **Third-cohort result/claim freeze:** [`docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`](docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md)
- **Machine-readable prospective H2 result:** [`results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`](results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json)
- **Machine-readable measurement/support receipt:** [`results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`](results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json)

## Frozen result snapshot

### H1 — observer-disjoint measurement validity

For the first-frozen repeated observer-disjoint reserve test:

- 200 partitions;
- median paired species = 329;
- median split Spearman rho = **0.7891**;
- 5th percentile rho = **0.7652**;
- median Lin CCC = **0.8548**;
- median Spearman-Brown reliability = **0.8821**.

A later deliberately stricter deterministic split yielded rho = **0.7927** and missed its prespecified 0.80 floor. The supported claim is reproducibility under the first-frozen validation rule, not near-perfect or split-invariant reliability.

### H2 — recurrent achromatic–chromatic geometry

The original discovery/reserve cohorts localized the recurrent construction-controlled signal to a fixed white-versus-equal-nonwhite contrast. Because that target was isolated after the original broad geometry had been opened, those cohorts are discovery/audit evidence for the named axis.

The prospective third cohort then tested the already frozen axis without retuning:

- outcome-blind candidate frame = 3,230 species;
- selected before biological opening = 500 species;
- terminal cohort = **499 species × 100 rows = 49,900 rows**;
- terminal partitions = **256 / 256**;
- classifiable rows = **25,788**;
- measurement-evaluable species = **377** (required >=250);
- replacements = **0**;
- support gate = **PASS**.

Primary 0.10 tier:

- vector species = **158**;
- observed W = **0.5172457461**;
- structured-null median = **0.4571428150**;
- p = **0.001**.

Strict 0.20 tier:

- vector species = **86**;
- observed W = **0.5329282123**;
- structured-null median = **0.4593196659**;
- p = **0.001**.

Frozen verdict:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

This is a **prospective species-disjoint confirmation within the same iNaturalist source/opportunity universe**, not an independent-source replication.

### H3a — broad phylogenetic signal

Reserve Blomberg-K tests were unsupported on all three frozen placement scenarios:

- S1 p = 0.2716;
- S2 p = 0.4134;
- S3 p = 0.2674.

Frozen verdict: `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

### H3b — sampled photographic span

Discovery: rho = 0.1798786, p = 0.00089996.

Species-disjoint reserve: rho = -0.0025855, p = 0.9586021.

Frozen verdict: `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

## Claim boundaries

The 42,111-species global frame is a sampling/opportunity frame, **not a prevalence denominator**.

The present evidence does not establish:

- global prevalence of flower-colour polymorphism;
- an independent-source H2 replication;
- pigment chemistry or pigment-loss/gain mechanism;
- evolutionary direction of white/nonwhite transitions;
- pollinator, climate or other adaptive causation;
- a recurrent non-white hue axis;
- absence of all phylogenetic structure;
- irrelevance of true biological range size.

P500 remains a successful prospective measurement-transport exercise with **no durable H2 biological verdict** because its post-calculation serialization failed before a terminal H2 result was written. It is not replayed or retroactively classified.

## Archived RGFCA shared-boundary programme

The earlier Repeated Global Flower-Colour Atlas programme remains available for provenance and negative/identifiability results:

- [research goal](docs/RGFCA_RESEARCH_GOAL.md)
- [research status](docs/RGFCA_RESEARCH_STATUS.md)
- [RGFCA manuscript](docs/RGFCA_MANUSCRIPT.md)
- [reserve replication results](docs/RGFCA_RESERVE_REPLICATION_RESULTS.md)
- [measurement/qualification limits](docs/RGFCA_ROI_QUALIFICATION_AUDIT.md)
- [global barrier-atlas protocol](docs/GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md)

It must not be mixed into the polymorphism paper as positive evidence for a shared global colour boundary.

### Legacy: six-species Chapter 1 held-out spatial analysis

- **Manuscript draft:** [`docs/JBI_CHAPTER1_MANUSCRIPT.md`](docs/JBI_CHAPTER1_MANUSCRIPT.md)
- **Current decision and execution status:** [`docs/JBI_CHAPTER1_SPATIAL_STATUS.md`](docs/JBI_CHAPTER1_SPATIAL_STATUS.md)
- **Frozen protocol:** [`docs/JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md`](docs/JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md)
- **Numerical results:** [`docs/JBI_CHAPTER1_RESULTS.md`](docs/JBI_CHAPTER1_RESULTS.md)
- **Figure plan and legends:** [`docs/JBI_CHAPTER1_FIGURE_PLAN.md`](docs/JBI_CHAPTER1_FIGURE_PLAN.md)
- **Canonical figure products:** [`docs/figures/jbi_ch1_figure_c1_stage_a_global.png`](docs/figures/jbi_ch1_figure_c1_stage_a_global.png) through the C4 and C-S2 products
- **Figure manifest:** [`docs/supporting/jbi_ch1_figure_manifest_v1.json`](docs/supporting/jbi_ch1_figure_manifest_v1.json)
- **Boundary CI:** [`.github/workflows/jbi-global-colour-boundaries.yml`](.github/workflows/jbi-global-colour-boundaries.yml)

### Legacy: frozen 34-species comparative paper

- **Manuscript:** [`docs/jbi_manuscript.md`](docs/jbi_manuscript.md)
- **Pipeline and evidence reduction:** [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md)
- **Figure plan:** [`docs/FIGURE_PLAN.md`](docs/FIGURE_PLAN.md)
- **Canonical figures:** [`docs/figures/`](docs/figures/)
- **Supporting Information map:** [`docs/jbi_supporting_information_index.md`](docs/jbi_supporting_information_index.md)
- **Remaining submission gates:** [`docs/jbi_submission_completion_checklist.md`](docs/jbi_submission_completion_checklist.md)
- **Canonical frozen input:** [`data/frozen/frozen_34species_five_metric_dataset.csv`](data/frozen/frozen_34species_five_metric_dataset.csv)
- **Reproduction workflow:** [`.github/workflows/34species-paper.yml`](.github/workflows/34species-paper.yml)

## Legacy Chapter 1: frozen photograph analysis

### Design

The frozen development set contains six species:

- *Antirrhinum majus*;
- *Dactylorhiza sambucina*;
- *Gentiana lutea*;
- *Ipomoea purpurea*;
- *Lysimachia arvensis*;
- *Raphanus sativus*.

For every species, 200 georeferenced photographs were acquired and assigned outcome-blind to 80 calibration and 120 evaluation observations. The complete design is therefore:

```text
1,200 photographs
   ├─ 480 calibration = 80/species
   └─ 720 held-out evaluation = 120/species
```

Calibration geometry did not justify forcing one universal discrete-morph scheme. The primary representation was frozen before evaluation as a **species-specific continuous colour vector**, standardized with calibration-only parameters. All 720 held-out photographs were processed successfully, with no localization failure, no final discrete biological label and no post-opening rule change.

### Ordered inference

```text
frozen continuous representation
        ↓
Stage A: within-species spatial organization
        ↓ prospective gate passed
Stage B: cross-species shared-transition concentration
        ↓ primary gate not passed
no confirmatory geographic-cause overlay
```

Stage A constructs a colour-blind spherical nearest-neighbour graph within each species and permutes complete vectors strictly within species. The primary equal-species result is:

- `k = 5`;
- observed Q = 1.39114;
- null mean = 1.42943;
- standardized clustering deficit = 2.3113;
- lower-tail `p = 0.0113`.

The direction is retained at `k = 3` (`p = 0.0066`) and `k = 8` (`p = 0.0065`).

Stage B uses label-blind geometry to define where transitions were detectable, ranks transition intensity within species and tests concentration of the shared surface under a complete species-conditioned permutation pipeline. The primary 500-km/36×18 analysis gives:

- observed concentration = 0.0082315;
- null mean = 0.0056757;
- standardized excess = 1.4389;
- upper-tail `p = 0.0906`.

The supported conclusion is therefore:

> Continuous flower colour is spatially organized within species, but the present six-species sample does not establish one universal global geography of the strongest transitions.

One coarser sensitivity configuration is nominally below 0.05, but the remaining supports are not; it is retained as exploratory scale sensitivity rather than substituted for the prospectively selected primary analysis.

### Chapter-1 production entry points

- `scripts/data/extract_jbi_ch1_florence_evaluation_features.py`
- `scripts/analysis/run_jbi_ch1_stage_a_continuous_graph.py`
- `scripts/analysis/run_jbi_ch1_stage_b_shared_transition.py`
- `scripts/analysis/make_jbi_ch1_spatial_figures.py`
- `scripts/analysis/make_jbi_ch1_spatial_figures_qa.py`

The governing rule is:

> Species may disappear from the map display, but species must never disappear from the null model.

## Legacy: frozen 34-species comparative paper

### Final paper dataset

The canonical statistical input is committed at:

`data/frozen/frozen_34species_five_metric_dataset.csv`

It is checksum-locked and contains:

- **34 species**;
- **25 plant families**;
- **20** within-population flower-colour polymorphism cases;
- **14** geographically structured flower-colour variation cases;
- minimum **20 occupied climate cells** per species;
- five symmetric climatic-niche metrics.

The labels are currently **source-traceable, rule-derived classifications**. Completed independent blinded human review is not claimed unless completed reviewer sheets are supplied.

### Evidence chain

```text
1,075 retained literature works
        ↓ species mapping + high-recall screening
664 candidate species (140 families)
        ↓ direct evidence screening
72-species initial review queue
        ↓ targeted follow-up + evidence aggregation
111-species resolved review queue
        ↓ unambiguous binary classification + climate eligibility
34 frozen species
   ├─ 20 within-population
   └─ 14 among-population
```

The remembered historical “~180” stage is not used as a formal manuscript count because the repository does not preserve a unique screening unit corresponding to that number.

A later systematic-map search used 15 query blocks and 52 shards and recovered 79,242 deduplicated bibliographic records. It is retained as broader search-completeness infrastructure; it is not presented as a deterministic parent of the original 34-species freeze, and its unreviewed expanded sets are not primary manuscript data.

Full provenance: [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md).

### Comparative pipeline

```text
literature discovery / provenance
        ↓
evidence screening + spatial classification
        ↓
GBIF occurrences + WorldClim occupied climates
        ↓
durable 34-species five-metric freeze
        ↓
five climatic-niche GLMs
        ↓
9,999 permutations + leave-one-family-out + collinearity
        ↓
OpenTree + dated phylogenetic sensitivity
        ↓
CR2/Satterthwaite + power/precision diagnostics
        ↓
canonical figures + manuscript + Supporting Information
```

Every primary metric uses:

```text
among ~ metric_z + effort_z
```

with family-clustered sandwich uncertainty, 9,999 label permutations and leave-one-family-out refits. Holm-adjusted probabilities across the five metrics are multiplicity context. VIF/condition-number diagnostics, OpenTree and time-scaled phylogenetic models, CR2/Satterthwaite inference and design-based power/precision simulation are sensitivity analyses.

All five climatic-niche point estimates are below one. Moisture breadth shows the largest contrast, but multiplicity-adjusted and phylogenetic intervals do not support a unique moisture mechanism. The paper emphasizes effect sizes and directional consistency: geographically structured colour variation tends to occur toward the narrower end of sampled occupied climatic niche breadth than within-population coexistence.

### Shared package and production entry points

- `fcp_pipeline/constants.py` — frozen metrics, counts and model specification
- `fcp_pipeline/evidence.py` — source-traceable spatial-evidence rules
- `fcp_pipeline/models.py` — standardized GLM, permutation and family-deletion helpers
- `fcp_pipeline/validation.py` — dataset, checksum and output invariants
- `scripts/run_34species_models.py`
- `scripts/run_34species_phylogenetic.R`
- `scripts/run_34species_power_precision.py`
- `scripts/run_34species_cr2.R`
- `scripts/make_paper_figures.py`

Install locally with:

```bash
python -m pip install -e .
```

## Repository boundary

The active repository mainline is the species-level flower-colour polymorphism paper defined by the frozen H1-H3 ledgers and the prospective third-cohort H2 result.

A file belongs to the active polymorphism mainline when it supports:

1. the 42,111-species sampling/opportunity frame or high-depth cohort construction;
2. the continuous four-state D phenotype and its observer-disjoint validation;
3. the label-free H2 geometry, white-axis audit or prospectively frozen q_white/W test;
4. the prospective third-cohort measurement, support and H2 chain of custody;
5. the bounded H3a/H3b replication tests;
6. manuscript, figure, audit or reproducibility products for those analyses.

RGFCA shared-boundary work, the six-species Chapter 1 lane and the frozen 34-species comparative paper remain recoverable legacy programmes. Their samples, estimands, protocols and claims are distinct and are not pooled with the active polymorphism inference merely because they coexist in this repository.
