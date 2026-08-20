# Frozen 34-species paper pipeline

This document is the canonical map from literature discovery to the Journal of Biogeography analyses. It records **units at every reduction step** so that works, candidate species and final model species are not conflated.

## 1. Scientific question

The paper asks whether species-level occupied climatic niche differs according to the documented spatial organization of natural intraspecific flower-colour variation:

- `within_population`: at least two discrete natural colour variants are explicitly documented as coexisting within a population;
- `among_population`: geographic or among-population colour differentiation is documented without retained evidence of local coexistence;
- `mixed`: both spatial signals are documented;
- `unclear`: retained evidence does not resolve the spatial organization.

Only the first two states enter the frozen binary comparison.

## 2. Verified data-reduction flow that produced the manuscript sample

The repository does **not** support an active, reproducible `180 -> 34` stage. Historical working lists existed, but the preserved QC records support the following exact chain and this is the chain that should be reported.

| Stage | Unit | Verified count | Role |
|---|---|---:|---|
| broad OpenAlex discovery retained after mapping/deduplication | works | 1,075 | literature-discovery pool used by the original evidence pipeline |
| species-level candidates linked to retained literature | species | 664 | high-recall candidate pool across 140 families |
| initial review queue | species | 72 | direct/high-priority evidence before follow-up |
| resolved review queue after targeted follow-up and evidence aggregation | species | 111 | evidence-supported cases carried forward after deferred-candidate rescue |
| final frozen binary climatic model set | species | 34 | 20 within-population + 14 among-population species from 25 families, each with >=20 occupied climate cells |

The initial 72-species queue is **not** a separate inferential sample. Follow-up searches rescued/clarified deferred candidates, producing the resolved 111-species evidence queue. Mixed, unclear, conflicting or otherwise non-binary cases were not forced into the final response variable. The final 34 were frozen independently of the final climate-model results.

### Later systematic-search expansion

A later systematic-map infrastructure (15 query blocks, 52 shards; PR #8) retrieved 79,242 deduplicated bibliographic records. This substantially broadens search coverage and is preserved as **search-completeness/provenance infrastructure**, but it was developed after the original 34-species evidence path. It must not be described as if the final 34 were literally selected by a single `79,242 -> ... -> 34` deterministic chain.

Exploratory expanded sets generated downstream of that newer corpus remained unreviewed and are not part of the current paper.

## 3. Evidence classification

The resolved evidence pipeline preserves source identifiers, titles and evidence passages. Spatial labels are rule-derived and source-traceable. A label is not assigned merely because the word `polymorphism` appears or because a study sampled multiple sites.

Binary inclusion requires an unambiguous retained spatial signal. If both within- and among-population evidence occur, the species is `mixed`; if neither can be established, it is `unclear`. Mixed and unclear cases are excluded from the binary manuscript analysis rather than coerced into one class.

The current repository does **not** document completed independent blinded human review for all 34 species. Until completed reviewer sheets exist, the manuscript must say `source-traceable, rule-derived classifications`.

## 4. Occurrence and climate construction

For the frozen species, coordinate-bearing GBIF records are linked to WorldClim 2.1 bioclimatic variables. Records with identical nine-variable climate vectors within a species are deduplicated to occupied climate cells. The primary frozen paper comparison requires at least 20 occupied climate cells per species.

The five symmetric species-level niche summaries are:

1. temperature breadth;
2. moisture breadth;
3. climatic heterogeneity;
4. PCA dispersion;
5. PCA hull area.

These are realised occupied-climate summaries. They are not morph-specific climatic tolerance or physiological niche estimates.

## 5. Primary statistical analysis

Each climatic metric is fitted separately with the same frozen 34 species:

`among ~ metric_z + effort_z`

where `among = 1` for geographically structured cases and `0` for within-population cases, `metric_z` is the standardized climatic metric, and `effort_z` is standardized `log1p(n_climate_cells)`.

Primary uncertainty and robustness:

- family-clustered sandwich covariance;
- 9,999 label permutations;
- leave-one-family-out refits;
- Holm multiplicity context across the five metrics;
- predictor-correlation, VIF and condition-number diagnostics.

## 6. Phylogenetic and finite-sample sensitivity

Phylogenetic sensitivity is reported in two complementary ways:

- Open Tree induced topology, 30 uniquely matched species, 100 polytomy resolutions, Grafen branch lengths, `phyloglm(logistic_MPLE)`;
- time-scaled `GBOTB.extended.LCVP` / V.PhyloMaker2 trees, all 34 species, placement scenarios S1-S3.

Finite-sample diagnostics include:

- CR2 family-cluster correction with Satterthwaite degrees of freedom;
- design-based power/precision simulation using the observed 34-species predictor and effort distributions.

These diagnostics quantify uncertainty at n=34; they are not used to declare the sample adequate post hoc.

## 7. Final result hierarchy

The five non-phylogenetic point estimates are all negative. Moisture breadth is the strongest observed contrast (OR about 0.41), but its Holm-adjusted and phylogenetic uncertainty does not support a unique moisture mechanism. The defensible conclusion is a directionally coherent comparative pattern: geographically structured flower-colour variation tends to occur toward the narrower end of sampled occupied climatic niche breadth, with the largest observed contrast along the moisture axis.

## 8. Canonical active code

Active paper code is organized by purpose:

- `fcp_pipeline/` — shared constants, model functions and hard validation gates;
- `scripts/run_34species_models.py` — five symmetric GLMs, permutations, LOFO and multiplicity;
- `analysis_34species_environmental_niche_phylogenetic.R` — collinearity + phylogenetic sensitivities;
- `scripts/run_34species_power_precision.py` — design-based finite-sample simulation;
- `analysis_34species_cluster_small_sample.R` — CR2/Satterthwaite sensitivity;
- `.github/workflows/34species-paper.yml` — canonical paper workflow.

Literature acquisition/search provenance is retained separately under `literature/` and the relevant source-data/QC files under `data/`. Historical theory and exploratory expanded-set analyses are not part of the active paper pipeline.

## 9. Hard invariants

The final workflow must fail if any of the following change without an explicit new freeze:

- species != 34;
- families != 25;
- within/among != 20/14;
- any of the five niche metrics is absent;
- fewer than 9,999 valid permutations are produced for a main model;
- an unreviewed expanded-set classification enters the primary dataset.
