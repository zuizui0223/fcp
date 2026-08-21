# Frozen 34-species paper pipeline

This is the canonical map from literature discovery to the Journal of Biogeography analysis. Every reduction step names its **unit** so works, candidate species, review queues and final model species are not conflated.

## 1. Scientific question

Does species-level occupied climatic niche differ according to the documented spatial organization of natural intraspecific flower-colour variation?

- `within_population`: at least two discrete natural colour variants are explicitly documented as coexisting within a population;
- `among_population`: geographic/among-population colour differentiation is documented without retained evidence of local coexistence;
- `mixed`: both signals are documented;
- `unclear`: retained evidence does not resolve spatial organization.

Only the first two enter the frozen binary comparison.

## 2. Verified reduction from literature to 34 species

The repository does **not** preserve an unambiguous reproducible `~180 -> 34` stage. The source-backed chain is:

| Stage | Unit | Verified count | Role |
|---|---|---:|---|
| broad OpenAlex discovery retained after mapping/deduplication | works | 1,075 | original literature-discovery pool |
| species linked to retained literature | species | 664 | high-recall candidate pool across 140 families |
| initial review queue | species | 72 | direct/high-priority evidence before follow-up |
| resolved review queue after targeted follow-up/evidence aggregation | species | 111 | resolved evidence queue |
| frozen binary climatic model set | species | 34 | 20 within + 14 among; 25 families; >=20 occupied climate cells |

The 72- and 111-species stages are screening/evidence layers, not separate inferential datasets. Mixed, unclear, conflicting and non-binary cases were not coerced into the final response. Climate-model results were not used to assign the frozen spatial label.

### Later systematic-search expansion

PR #8 added a broader systematic-map infrastructure: 15 query blocks, 52 shards and 79,242 deduplicated bibliographic records. It strengthens search-completeness provenance but was developed after the original 34-species evidence path. It must not be described as a literal `79,242 -> ... -> 34` deterministic chain. Unreviewed expanded species sets generated from it are not current primary data.

## 3. Evidence classification

The resolved evidence layer preserves source identifiers and evidence passages. Active rule definitions and normalization live in `fcp_pipeline/evidence.py`.

- evidence supporting local coexistence only -> `within_population`;
- evidence supporting geographic differentiation only -> `among_population`;
- both -> `mixed`;
- neither resolved -> `unclear`.

The current labels are **source-traceable, rule-derived classifications**. The repository contains review scaffolding, but completed independent blinded human review must not be claimed unless completed reviewer files actually exist.

## 4. Occurrence and climate construction

Coordinate-bearing GBIF records were linked to WorldClim 2.1 bioclimatic variables. Within species, repeated identical nine-variable climate vectors were deduplicated to occupied climate cells. At least 20 occupied climate cells are required.

Five species-level summaries are evaluated symmetrically:

1. temperature breadth;
2. moisture breadth;
3. climatic heterogeneity;
4. PCA dispersion;
5. PCA hull area.

They are realised occupied-climate summaries, not morph-specific physiological tolerance.

## 5. Durable production freeze

The production statistical input is committed at:

`data/frozen/frozen_34species_five_metric_dataset.csv`

SHA-256:

`bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed`

`data/frozen/freeze_manifest.json` records its recovery provenance. The file contains only the columns required to reproduce the current comparative analysis: species, family, binary spatial label, classification source, climate-cell effort and all five climatic-niche metrics. Source passages remain in the evidence/classification provenance tables instead of being duplicated in the model matrix.

The canonical row order is `canonical_name` ascending. `fcp_pipeline/validation.py` locks both the SHA-256 and biological counts. The analysis helper also canonicalizes species order before finite Monte Carlo permutations, so results no longer depend on caller row ordering.

The freeze was recovered while the historical Actions artifact was still available and then committed durably. The historical short-retention artifact has since expired and is **not a runtime dependency**.

## 6. Primary statistical analysis

Each metric is fitted separately to the same 34 species:

`among ~ metric_z + effort_z`

where `among = 1` for geographically structured cases, `metric_z` is the standardized climatic metric and `effort_z = z(log1p(n_climate_cells))`.

Primary uncertainty/robustness:

- family-clustered sandwich covariance;
- 9,999 label permutations, seed `20260719`;
- leave-one-family-out refits;
- Holm multiplicity context across five metrics;
- predictor correlation, VIF and condition number.

The durable canonical permutation p-values are:

| metric | permutation p |
|---|---:|
| temperature breadth | 0.6131 |
| moisture breadth | 0.0423 |
| climatic heterogeneity | 0.3567 |
| PCA dispersion | 0.3859 |
| PCA hull area | 0.2372 |

The corresponding moisture Holm-adjusted permutation p-value is 0.2115. The shift from earlier finite-draw Monte Carlo values reflects explicit canonicalization of row order; coefficient estimates, clustered uncertainty, LOFO and biological interpretation are unchanged.

## 7. Phylogenetic and finite-sample sensitivity

Phylogenetic sensitivity uses:

- Open Tree induced topology, 30 uniquely matched species, 100 random polytomy resolutions, Grafen branch lengths, `phyloglm(logistic_MPLE)`;
- time-scaled `GBOTB.extended.LCVP` / V.PhyloMaker2 trees, all 34 species, scenarios S1-S3.

Finite-sample diagnostics use:

- CR2 family-cluster correction with Satterthwaite degrees of freedom;
- design-based power/precision simulation retaining the observed 34-species predictor and effort distributions.

These quantify uncertainty at n=34; they are not post-hoc criteria for declaring sample adequacy.

## 8. Result hierarchy

All five non-phylogenetic point estimates are negative. Moisture breadth has the largest observed contrast (OR about 0.41), but Holm-adjusted and phylogenetic uncertainty does not justify a unique moisture mechanism. The defensible result is a directionally coherent comparative pattern: geographically structured flower-colour variation tends to occur toward the narrower end of sampled occupied climatic niche breadth, with the largest observed contrast along moisture.

## 9. Active code by methodological layer

- `fcp_pipeline/` — shared constants, evidence rules, model functions and hard validation;
- `scripts/literature/` — literature acquisition/follow-up provenance;
- `scripts/occurrence/` — GBIF/WorldClim construction utilities;
- `scripts/run_34species_models.py` — five symmetric GLMs, permutation, LOFO and multiplicity;
- `scripts/run_34species_phylogenetic.R` — collinearity + phylogenetic sensitivities;
- `scripts/run_34species_power_precision.py` — design-based finite-sample simulation;
- `scripts/run_34species_cr2.R` — CR2/Satterthwaite sensitivity;
- `scripts/submission/` — retained provenance/citation utilities;
- `.github/workflows/34species-paper.yml` — canonical paper workflow.

Historical phase theory, matched controls, fragmentation/turnover experiments and unreviewed expanded-set ecology are not part of the active pipeline.

## 10. Hard invariants

The production workflow fails if any of the following change without an explicit new freeze:

- dataset SHA-256 changes;
- species != 34;
- families != 25;
- within/among != 20/14;
- `classification_source` differs from `baseline_unambiguous`;
- any of the five metrics is absent;
- any species has <20 climate cells;
- fewer than 9,999 valid permutations are produced;
- numerical regression of the five ORs or canonical Monte Carlo p-values drifts.
