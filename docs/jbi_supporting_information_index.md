# Journal of Biogeography Supporting Information index

This index describes the **active frozen 34-species paper**. Historical S1–S19 numbering from earlier development stages is no longer used because it mixed superseded analyses, transient workflow artifacts and submission-facing provenance files.

## 1. Frozen analysis input

- `data/frozen/frozen_34species_five_metric_dataset.csv`
  - canonical downstream input for the paper;
  - 34 species, 25 families, 20 within-population / 14 geographically structured cases;
  - all species have at least 20 occupied climate cells;
  - contains all five climatic-niche metrics used in the main comparison;
  - SHA-256: `bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed`.
- `data/frozen/freeze_manifest.json`
  - records the freeze rule, checksum, counts and historical provenance of the durable input.

The committed frozen CSV, rather than a temporary GitHub Actions artifact, is the source consumed by the canonical manuscript workflow.

## 2. Classification and evidence provenance

- `docs/supporting/frozen_classification_manifest.csv`
  - 34-species classification manifest used to document the frozen binary states.
- `docs/supporting/classification_correction_log.csv`
  - records corrections to the frozen classification layer; the current file contains no additional post-freeze corrections.
- `docs/supporting/blinded_classification_review.csv`
  - blinded 34-species review sheet prepared for independent human review; blank reviewer fields must not be represented as completed review.
- `docs/supporting/rule_classification_key.csv`
  - separate rule-derived comparison key for use after blinded review.
- `docs/jbi_classification_review_protocol.md`
  - operational instructions for blinded review and adjudication.
- `docs/jbi_classification_rule_audit.md`
  - audit of rule wording and code parity.
- `docs/jbi_literature_search_provenance.md`
  - reconstructs the original literature-discovery and targeted follow-up path used to reach the frozen evidence set.

The manuscript therefore describes the current labels as **source-traceable, rule-derived classifications**. Completed independent blinded review is not claimed unless reviewer fields are actually completed and adjudicated.

## 3. Phylogenetic sensitivity inputs

- `docs/supporting/jbi_opentree_induced_topology.tre`
  - induced Open Tree topology used for the 30-species topology-based sensitivity analysis.
- `docs/supporting/jbi_dated_phylogeny_s1.tre`
- `docs/supporting/jbi_dated_phylogeny_s2.tre`
- `docs/supporting/jbi_dated_phylogeny_s3.tre`
  - time-scaled V.PhyloMaker2 trees retaining all 34 species under placement scenarios S1–S3.
- `docs/supporting/jbi_dated_phylogeny_manifest.json`
  - dated-tree generation and placement provenance.

The Open Tree analysis uses Grafen branch lengths and 100 polytomy resolutions. The dated-tree analysis uses the `GBOTB.extended.LCVP` backbone. These are sensitivity analyses, not alternative primary datasets.

## 4. GBIF and occurrence provenance

- `docs/supporting/gbif_taxon_resolution_audit.csv`
  - species-level GBIF taxonomic-resolution audit retained for the frozen taxa.
- `docs/supporting/jbi_gbif_doi_bundle/`
  - preparatory material for external GBIF occurrence archiving and citation.
- `docs/jbi_gbif_doi_protocol.md`
  - procedure for obtaining the final citable GBIF occurrence identifier / Derived Dataset registration.

The repository still contains historical paginated-retrieval QC files for provenance. They are **not part of the active five-metric model chain** and should not be cited as if the current frozen climatic metrics were regenerated from the paginated sensitivity dataset.

A permanent citable GBIF identifier remains **Not verified** until authenticated external registration is completed.

## 5. Analysis documentation and submission-facing summaries

- `docs/supporting/jbi_environmental_niche_comprehensive_34species.md`
  - human-readable record of the five symmetric models, collinearity diagnostics and phylogenetic sensitivities.
- `docs/supporting/cr2_satterthwaite_summary.csv`
  - verified five-metric CR2/Satterthwaite summary copied from the canonical workflow output; retains exact odds ratios, confidence intervals, Satterthwaite degrees of freedom and p-values.
- `.github/workflows/34species-paper.yml`
  - canonical executable workflow.
- `scripts/run_34species_models.py`
  - five family-clustered GLMs, 9,999 row-order-invariant label permutations, Holm context and leave-one-family-out refits.
- `scripts/run_34species_phylogenetic.R`
  - collinearity diagnostics and Open Tree / dated phylogenetic sensitivity.
- `scripts/run_34species_cr2.R`
  - CR2/Satterthwaite finite-cluster sensitivity.
- `scripts/run_34species_power_precision.py`
  - design-based power/precision diagnostics.

The workflow validates the frozen-data checksum, 34/25/20/14 counts, five metrics, 9,999 valid permutations for every main model, numerical regression of the primary effect estimates, phylogenetic outputs, CR2/Satterthwaite outputs and power/precision diagnostics. Workflow-run and artifact identifiers are intentionally **not** treated as durable Supporting Information metadata because they change across valid reruns and artifacts expire.


## 6. Canonical figures

- `docs/figures/figure1_geographic_context.png` / `.pdf`
  - main-text geographic context for the 34 focal species; broader exact GBIF occurrence subset is shown only as context/QC.
- `docs/figures/figure2_five_metric_forest.png` / `.pdf`
  - central five-metric family-clustered effect-size figure generated directly from the frozen 34-species dataset.
- `docs/figures/figure3_raw_species_metrics.png` / `.pdf`
  - all 34 species shown across standardized versions of the five frozen climatic metrics; visualization only.
- `docs/figures/figureS1_34_species_distribution_context.png` / `.pdf`
  - species-specific geographic occurrence maps using the broader exact GBIF subset as supporting distribution context/QC.
- `scripts/make_paper_figures.py`
  - canonical figure generator.
- `docs/FIGURE_PLAN.md`
  - explains why these figures were selected from the pipeline and why robustness diagnostics remain supporting material.

Figure 1 and Supporting Figure S1 must retain the distinction between the broader exact GBIF citation subset and the exact occurrence sampling path that produced the frozen primary climatic summaries.

## Submission boundary

Material from matched-control, unreviewed expanded-set, range-fragmentation and environmental-turnover experiments is recoverable from Git history and closed development PRs but is not part of the active submission analysis. The Supporting Information for the paper should be generated from the durable frozen dataset and the canonical workflow above, so historical exploratory outputs are not mixed with the final 34-species inference.