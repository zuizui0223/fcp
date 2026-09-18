# FCP — flower-colour variation across species and space

This repository is the **geographic-space arm** of a broader programme on the spatiotemporal organization of flower-colour variation. `fcp` now contains an active global flower-colour polymorphism mainline together with two older frozen spatial/comparative lanes.

## Active mainline — global flower-colour polymorphism

### What RGFCA was, and what changed

**RGFCA = Repeated Global Flower-Colour Atlas.** It was the upstream image-first sampling and measurement programme that asked whether independent plant species repeatedly place strong within-species flower-colour discontinuities in the same broad geographic regions. Its core design used balanced repeated world-map realizations, equal-species weighting, geographic opportunity correction and species-conditioned colour permutations.

The shared-geography estimand did not become the positive core of the present paper: the primary recurrent-field test was not supported at the prespecified level, cross-species boundary commonness did not transfer, and later sharedness qualification exposed an identifiability limitation. What did remain strongly useful was the infrastructure and the species-specific heterogeneity: globally broad species discovery, high-depth discovery/reserve cohorts, location-blind colour measurement, paired background controls and per-species spatial-organization statistics.

The current paper therefore changes the level of generality from **a shared place on the world map** to **shared structure in phenotype space plus species-specific spatial organization**. A full interpretation is frozen in `docs/RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md`.

The active paper asks three linked questions: **can within-species flower-colour polymorphism be measured reproducibly as a species phenotype, is that variation geometrically constrained in colour space, and do more polymorphic species organize that variation more strongly across geography?**

The current frozen result is:

- observer-disjoint high-depth validation supports reproducible species-level polymorphism `D`, while a later stricter deterministic split shows that reliability is not near-perfect or split-invariant;
- the original discovery/reserve geometry localized to a white-versus-nonwhite achromatic–chromatic axis;
- an already frozen version of that axis was then tested in a pre-frozen species-disjoint third cohort from the same iNaturalist opportunity universe;
- third-cohort measurement completed **49,900 rows from 499 species**, with **377** species passing the predeclared measurement-support gate;
- H2 passed at the primary 0.10 tier (**158 species, W = 0.51725, p = 0.001**) and strict 0.20 tier (**86 species, W = 0.53293, p = 0.001**);
- frozen verdict: `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`;
- across the original high-depth discovery/reserve cohorts, greater D is also associated with stronger within-species geographic colour organization; the reserve association survives sampled-span plus clear technical-failure adjustment, a matched flower-minus-background contrast and ambiguity-endpoint stress tests;
- broad reserve phylogenetic signal and the discovery sampled-span association are not supported as simple general explanations.

The prospective H2 result is **species-disjoint within the same iNaturalist source/opportunity universe**. It is not an independent-source replication, a global prevalence estimate, or a pigment/pollinator/climate mechanism result. The earlier P500 run retains no durable H2 biological verdict and is not rescued retrospectively.

### Reusable methods package

The reusable, flower-colour-independent inference layer now lives in:

- **`packages/disttrait/`** — Python package for validated species-level distributional trait inference;
- **`docs/DISTTRAIT_PACKAGE_V0_7_20260919.md`** — active v0.7 validation/direction-heterogeneity boundary;
- **`docs/POLYMORPHISM_METHODS_CLASSIFICATION_20260918.md`** — standard statistics versus study-specific inference architecture.

`disttrait` currently exposes categorical diversity, observer-disjoint reliability, Hellinger/two-mode geometry, fixed-contrast alignment, construction-preserving nulls, categorical and scalar-continuous species-specific spatial organization, equal-species omnibus inference, matched focal-minus-background structure, and distribution–spatial association. v0.7 retains exact FCP-equivalence fixtures, false-positive/power and model-comparator surfaces, the external 1,600-tree non-flower empirical transport, and adds a direction-heterogeneity benchmark showing that common signed-slope models and direction-invariant distance–dissimilarity analyses answer different cross-species questions. Flower-colour acquisition/segmentation and the frozen white/non-white target remain application-specific.

### Start here

- **RGFCA → polymorphism interpretation:** [`docs/RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md`](docs/RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md)
- **Canonical manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT.md`](docs/POLYMORPHISM_MANUSCRIPT.md)
- **New Phytologist submission draft:** [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md)
- **Current claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
- **Post-confirmation architecture:** [`docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260918.md`](docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260918.md)
- **Figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
- **Supporting evidence map:** [`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`](docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md)
- **Submission-readiness audit:** [`docs/POLYMORPHISM_NEW_PHYTOLOGIST_SUBMISSION_READINESS_20260918.md`](docs/POLYMORPHISM_NEW_PHYTOLOGIST_SUBMISSION_READINESS_20260918.md)
- **Canonical figures:** [`docs/figures/polymorphism_20260918/`](docs/figures/polymorphism_20260918/)
- **Full analysis-history branch:** `analysis/h2-third-cohort-preopening-20260916`

## Retained frozen lanes

The repository also preserves **two frozen inferential lanes** whose samples, estimands and claims are distinct from the active polymorphism paper.

1. **Chapter 1 spatial-photograph lane:** present-day continuous flower-colour observations, within-species spatial organization and cross-species transition geography.
2. **34-species comparative lane:** occupied climatic-niche breadth versus literature-documented local coexistence or geographic colour differentiation.

Their samples, response variables, null models and claims are distinct. These lanes remain valid provenance and separate papers; they are not pooled with the active global polymorphism inference.

### Retained-lane canonical entry points

#### Chapter 1 spatial-photograph lane

- **Manuscript:** [`docs/JBI_CHAPTER1_MANUSCRIPT.md`](docs/JBI_CHAPTER1_MANUSCRIPT.md)
- **Status:** [`docs/JBI_CHAPTER1_SPATIAL_STATUS.md`](docs/JBI_CHAPTER1_SPATIAL_STATUS.md)
- **Frozen protocol:** [`docs/JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md`](docs/JBI_CHAPTER1_SPATIAL_STATE_DISTRIBUTION_PROTOCOL.md)
- **Numerical results:** [`docs/JBI_CHAPTER1_RESULTS.md`](docs/JBI_CHAPTER1_RESULTS.md)
- **Figure plan:** [`docs/JBI_CHAPTER1_FIGURE_PLAN.md`](docs/JBI_CHAPTER1_FIGURE_PLAN.md)

#### Frozen 34-species comparative lane

- **Manuscript:** [`docs/jbi_manuscript.md`](docs/jbi_manuscript.md)
- **Pipeline:** [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md)
- **Supporting Information:** [`docs/jbi_supporting_information_index.md`](docs/jbi_supporting_information_index.md)
- **Submission checklist:** [`docs/jbi_submission_completion_checklist.md`](docs/jbi_submission_completion_checklist.md)
- **Frozen dataset:** [`data/frozen/frozen_34species_five_metric_dataset.csv`](data/frozen/frozen_34species_five_metric_dataset.csv)
- **Reproduction workflow:** [`.github/workflows/34species-paper.yml`](.github/workflows/34species-paper.yml)

## Chapter 1: frozen photograph analysis

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

## Frozen 34-species comparative paper

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

The retained legacy material below remains assigned to its original lane and must not be pooled with the active polymorphism paper.

### Chapter 1 spatial lane

A file belongs here only if it supports:

1. frozen photograph acquisition, calibration or evaluation;
2. species-conditioned local colour organization;
3. label-blind transition detectability and shared concentration;
4. audit, figures or manuscript reporting for those analyses.

### 34-species comparative lane

A file belongs here only if it supports:

1. discovery and classification of documented flower-colour cases;
2. construction of the frozen 34-species occupied-climate dataset;
3. the five comparative niche models and required robustness analyses;
4. reproduction, audit or submission of that paper.

Historical phase-theory work, unreviewed expanded-set experiments and exploratory geographic-cause overlays are not promoted into either confirmatory main line merely because they remain recoverable in repository history.
