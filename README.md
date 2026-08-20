# FCP — spatial organization of intraspecific flower-colour variation

This repository supports a comparative Journal of Biogeography paper asking:

> **Does occupied climatic niche differ according to whether intraspecific flower-colour variation is maintained as local within-population coexistence or expressed as geographic differentiation among populations?**

The active paper is a **34-species frozen comparative analysis**, not a universal meta-analysis of all angiosperms and not the older mathematical phase-theory project that originally occupied this repository.

## Final paper dataset

- **34 species**
- **25 plant families**
- **20** within-population flower-colour polymorphism cases
- **14** geographically structured flower-colour variation cases
- minimum **20 occupied climate cells** per species
- five symmetric climatic-niche metrics

The labels are currently **source-traceable, rule-derived classifications**. Completed independent blinded human review is not claimed unless reviewer sheets are actually added.

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
        ↓ unambiguous spatial classification + climate eligibility
34 frozen binary species
   ├─ 20 within-population
   └─ 14 among-population
```

This is the reproducible chain supported by repository QC. A remembered historical “~180” working stage is **not used as a manuscript count** because the active repository does not preserve it as a uniquely defined final screening unit.

A later systematic-map search (15 query blocks, 52 shards, 79,242 deduplicated records) is retained as broader search-completeness infrastructure. It was developed after the original 34-species evidence path and does not replace the frozen manuscript dataset with its unreviewed expanded sets.

Full details: [`docs/PIPELINE_34SPECIES.md`](docs/PIPELINE_34SPECIES.md).

## Active analysis pipeline

```text
literature discovery / provenance
        ↓
evidence screening + spatial classification
        ↓
GBIF occurrences + WorldClim occupied climates
        ↓
34-species frozen dataset
        ↓
five climatic-niche GLMs
        ↓
permutation + LOFO + collinearity
        ↓
OpenTree + dated phylogenetic sensitivity
        ↓
CR2/Satterthwaite + power/precision diagnostics
        ↓
main manuscript + Supporting Information
```

### Shared code

- `fcp_pipeline/constants.py` — frozen metrics/counts/model specification
- `fcp_pipeline/models.py` — standardized model preparation, GLM, permutation and LOFO helpers
- `fcp_pipeline/validation.py` — hard dataset/model invariants

### Executable entry points

- `scripts/run_34species_models.py`
- `analysis_34species_environmental_niche_phylogenetic.R`
- `scripts/run_34species_power_precision.py`
- `analysis_34species_cluster_small_sample.R`

The canonical CI entry point is `.github/workflows/34species-paper.yml`.

## Statistical design

Every main metric uses the same formula:

```text
among ~ metric_z + effort_z
```

with family-clustered sandwich uncertainty, 9,999 label permutations and leave-one-family-out refits. Holm-adjusted p-values across five metrics are reported as multiplicity context. VIF/condition-number diagnostics, OpenTree and time-scaled phylogenetic models, CR2/Satterthwaite small-cluster inference and design-based power/precision simulation are sensitivity analyses.

## Current result

All five climatic-niche point estimates are negative. Moisture breadth shows the largest contrast (non-phylogenetic OR ≈ 0.41), but multiplicity-adjusted and phylogenetic intervals do not support claiming a unique moisture mechanism. The paper therefore emphasizes **directional consistency and effect sizes**: geographically structured colour variation tends to occur in species with narrower sampled occupied climatic niches than within-population coexistence.

## Repository boundary

Active root-level material should answer one of four questions:

1. How were candidate flower-colour cases found and screened?
2. Why do the 34 frozen classifications enter the paper?
3. How are climatic niches and statistical models reproduced?
4. How are submission-facing outputs validated?

Older phase-theory experiments, alternative sample expansions and abandoned modelling branches are preserved in Git history or explicitly archived when scientifically useful; they are not part of the active paper pipeline.
