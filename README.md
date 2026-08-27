# FCP — spatial organization of intraspecific flower-colour variation

This repository is the **geographic-space arm** of a broader programme on the spatiotemporal organization of flower-colour variation. `fcp` asks how intraspecific colour diversity is maintained or sorted across space; the complementary [`chun`](https://github.com/zuizui0223/chun) project asks how similar flower-colour states are repeatedly generated through evolutionary time. See [`docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md`](docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md).

This repository supports a comparative Journal of Biogeography paper asking:

> **Does occupied climatic niche differ according to whether intraspecific flower-colour variation occurs as local within-population coexistence or as geographic differentiation among populations?**

The active paper is a **frozen 34-species comparative analysis**. Historical mathematical phase-theory work and unreviewed expanded-set experiments are no longer part of the active analysis path.

## Start here

- **Programme position:** [`docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md`](docs/FLOWER_COLOUR_VARIATION_SPATIAL_PROGRAM.md)
- **Manuscript:** [`docs/jbi_manuscript.md`](docs/jbi_manuscript.md)
- **Pipeline and evidence reduction:** [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md)
- **Figure plan:** [`docs/FIGURE_PLAN.md`](docs/FIGURE_PLAN.md)
- **Canonical figures:** [`docs/figures/`](docs/figures/)
- **Supporting Information map:** [`docs/jbi_supporting_information_index.md`](docs/jbi_supporting_information_index.md)
- **Remaining submission gates:** [`docs/jbi_submission_completion_checklist.md`](docs/jbi_submission_completion_checklist.md)
- **Canonical frozen input:** [`data/frozen/frozen_34species_five_metric_dataset.csv`](data/frozen/frozen_34species_five_metric_dataset.csv)
- **Reproduction workflow:** [`.github/workflows/34species-paper.yml`](.github/workflows/34species-paper.yml)

## Final paper dataset

The canonical statistical input is committed at:

`data/frozen/frozen_34species_five_metric_dataset.csv`

It is checksum-locked and contains:

- **34 species**
- **25 plant families**
- **20** within-population flower-colour polymorphism cases
- **14** geographically structured flower-colour variation cases
- minimum **20 occupied climate cells** per species
- five symmetric climatic-niche metrics

The labels are currently **source-traceable, rule-derived classifications**. Completed independent blinded human review is not claimed unless completed reviewer sheets are actually supplied.

## How the 34 species were reached

The exact preserved evidence chain is:

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

A later systematic-map search used 15 query blocks and 52 shards and recovered 79,242 deduplicated bibliographic records. It is retained as broader search-completeness infrastructure; it is not presented as a direct deterministic parent of the original 34-species freeze, and its unreviewed expanded sets are not primary manuscript data.

Full provenance: [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md).

## Active analysis pipeline

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
9,999 permutations + LOFO + collinearity
        ↓
OpenTree + dated phylogenetic sensitivity
        ↓
CR2/Satterthwaite + power/precision diagnostics
        ↓
canonical figures + manuscript + Supporting Information
```

### Shared package

- `fcp_pipeline/constants.py` — frozen metrics/counts/model specification
- `fcp_pipeline/evidence.py` — source-traceable spatial-evidence rules and normalization
- `fcp_pipeline/models.py` — standardized GLM, permutation and LOFO helpers
- `fcp_pipeline/validation.py` — dataset, checksum and output invariants

Install locally with:

```bash
python -m pip install -e .
```

### Production entry points

- `scripts/run_34species_models.py`
- `scripts/run_34species_phylogenetic.R`
- `scripts/run_34species_power_precision.py`
- `scripts/run_34species_cr2.R`
- `scripts/make_paper_figures.py`

The canonical CI/reproduction entry is `.github/workflows/34species-paper.yml`. It reads the committed freeze directly and therefore does not depend on expiring Actions artifacts.

## Statistical design

Every main metric uses the same formula:

```text
among ~ metric_z + effort_z
```

with family-clustered sandwich uncertainty, 9,999 label permutations and leave-one-family-out refits. The permutation implementation first canonicalizes species order so finite Monte Carlo p-values are invariant to input row ordering. Holm-adjusted p-values across five metrics are reported as multiplicity context. VIF/condition-number diagnostics, OpenTree and time-scaled phylogenetic models, CR2/Satterthwaite inference and design-based power/precision simulation are sensitivity analyses.

## Current result

All five climatic-niche point estimates are negative. Moisture breadth shows the largest contrast (non-phylogenetic OR ≈ 0.41), but multiplicity-adjusted and phylogenetic intervals do not support claiming a unique moisture mechanism. The paper therefore emphasizes **effect sizes and directional consistency**: geographically structured colour variation tends to occur toward the narrower end of sampled occupied climatic niche breadth than within-population coexistence.

## Canonical paper figures

The main figures are selected from the biological question backward, not from whichever test gives the smallest p-value:

1. `figure1_geographic_context` — geographic scope of the 34 focal species; broader exact GBIF records are context/QC, not the exact primary climate-metric sample.
2. `figure2_five_metric_forest` — central result: all five production-model odds ratios on one common effect-size display.
3. `figure3_raw_species_metrics` — the 34 species themselves across all five standardized climatic metrics.
4. `figureS1_34_species_distribution_context` — supporting species-by-species geographic occurrence audit.

See [`docs/FIGURE_PLAN.md`](docs/FIGURE_PLAN.md) for the selection logic and interpretation boundary.

## Repository boundary

Active material should answer one of four questions:

1. How were documented flower-colour cases found and classified?
2. How were the frozen 34 species and their occupied climates constructed?
3. What do the five symmetric comparative models and required robustness checks show?
4. What is needed to reproduce, audit or submit that paper?

Exploratory analyses that do not answer one of these questions should remain outside the active path or be recoverable from Git history rather than duplicated in the submission tree.
