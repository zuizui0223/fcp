# FCP — Repeated Global Flower-Colour Atlas (RGFCA)

This repository is the **geographic-space arm** of a broader programme on the spatiotemporal organization of flower-colour variation. `fcp` asks how intraspecific colour diversity is maintained or sorted across space; the complementary [`chun`](https://github.com/zuizui0223/chun) project asks how similar flower-colour states are repeatedly generated through evolutionary time. See [`docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md`](docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md).

The active research mainline is **RGFCA**: image-first measurement of global flower-colour variation, balanced repeated world maps and species-conditioned spatial and ecological inference. The six-species Chapter 1 analysis and 34-species literature comparison are **legacy studies**, retained for provenance and reproduction rather than used as the biological foundation of the new atlas.

RGFCA has measured **50,000 photographs from 500 species**. Its eligible inferential frame contains **21,424 classifiable photographs from 369 species**. Balanced realizations each use 250 species and 20 photos per species, repeated 200 times with the full schedule preserved in the species-conditioned null. Repetitions measure sampling stability, not additional biological replication.

The research objective is an ecologically interpretable observed signal that survives independent validation. **A new exploratory 369-species test detects a small within-species distance–colour association: mean Spearman rho = 0.0270, permutation p = 0.001.** This is a candidate photo-derived spatial signal, not a shared-boundary or causal result. Completed non-support decisions remain unchanged: primary G1 concentration `p = 0.070`, species-disjoint commonness `p = 0.856`. The next stage is a prospectively fixed replication in the other **500 species / 50,000 previously unmeasured candidate photographs**, with disjoint observation/photo IDs and seasonal/observer sensitivity tests. See the [current status](docs/RGFCA_RESEARCH_STATUS.md), [replication protocol](docs/RGFCA_RESERVE_REPLICATION.md) and [research goal](docs/RGFCA_RESEARCH_GOAL.md).

## Start here

The two frozen inferential lanes retained below are now legacy studies. Their samples, response variables, null models and claims are distinct.

- **Active goal and validation route:** [`docs/RGFCA_RESEARCH_GOAL.md`](docs/RGFCA_RESEARCH_GOAL.md)
- **Current results, execution and next questions:** [`docs/RGFCA_RESEARCH_STATUS.md`](docs/RGFCA_RESEARCH_STATUS.md)
- **Active discovery manuscript (validation pending):** [`docs/RGFCA_MANUSCRIPT.md`](docs/RGFCA_MANUSCRIPT.md)
- **Discovery map, effects and figure provenance:** [`docs/RGFCA_PUBLICATION_FIGURES.md`](docs/RGFCA_PUBLICATION_FIGURES.md)
- **Verified real flower photo bar and source credits:** [`docs/RGFCA_PHOTO_BAR.md`](docs/RGFCA_PHOTO_BAR.md)
- **Supporting numerical evidence and submission gaps:** [`docs/RGFCA_SUPPORTING_EVIDENCE.md`](docs/RGFCA_SUPPORTING_EVIDENCE.md)
- **Audited core literature and claim limits:** [image/ecology precedents](docs/RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md), [statistical interpretation](docs/RGFCA_STATISTICAL_LITERATURE_AUDIT.md)
- **Core scientific-software citation audit:** [`docs/RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md`](docs/RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md)
- **Measurement providers, model identities and reuse limits:** [`docs/RGFCA_MEASUREMENT_PROVIDER_AUDIT.md`](docs/RGFCA_MEASUREMENT_PROVIDER_AUDIT.md)
- **Training data and limits of flower-region validation:** [JRC source and rights](docs/RGFCA_TRAINING_SOURCE_AUDIT.md), [100-image qualification audit](docs/RGFCA_ROI_QUALIFICATION_AUDIT.md)
- **RGFCA protocol:** [`docs/GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md`](docs/GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md)
- **Repeated-atlas method:** [`docs/REPEATED_GLOBAL_FLOWER_COLOUR_ATLAS_METHOD.md`](docs/REPEATED_GLOBAL_FLOWER_COLOUR_ATLAS_METHOD.md)
- **Methodological contribution:** [`docs/RGFCA_METHOD_NOVELTY_POSITIONING.md`](docs/RGFCA_METHOD_NOVELTY_POSITIONING.md)
- **Observation-bias audit:** [`docs/RGFCA_OBSERVATION_BIAS.md`](docs/RGFCA_OBSERVATION_BIAS.md)
- **Programme position:** [`docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md`](docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md)

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

Active RGFCA material belongs to the repeated-atlas research programme described above. Discovery, technical qualification, exploratory inference and independent biological validation have different claims. The two legacy studies below retain their original samples, protocols and results, and are not pooled with RGFCA.

### Legacy Chapter 1 spatial lane

A file belongs here only if it supports:

1. frozen photograph acquisition, calibration or evaluation;
2. species-conditioned local colour organization;
3. label-blind transition detectability and shared concentration;
4. audit, figures or manuscript reporting for those analyses.

### Legacy 34-species comparative lane

A file belongs here only if it supports:

1. discovery and classification of documented flower-colour cases;
2. construction of the frozen 34-species occupied-climate dataset;
3. the five comparative niche models and required robustness analyses;
4. reproduction, audit or submission of that paper.

Historical phase-theory work, unreviewed expanded-set experiments and exploratory geographic-cause overlays are not promoted into either confirmatory main line merely because they remain recoverable in repository history.
