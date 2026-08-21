# Journal of Biogeography analysis protocol — frozen 34-species design

## Central question

Do documented cases of natural intraspecific flower-colour variation differ in species-level occupied climatic niche according to whether colour variants coexist within populations or are geographically structured among populations?

The analysis is organized by **methodological layer**, not by the chronology in which exploratory scripts were developed.

## Layer 1 — Evidence discovery and screening

The manuscript sample comes from a broad literature-derived evidence pipeline. The preserved original chain is:

`1,075 retained works -> 664 candidate species -> 72-species initial review queue -> 111-species resolved queue -> 34 frozen binary species`.

The later 79,242-record systematic-map search is additional search-completeness infrastructure and must not be described as the direct deterministic parent of the original 34-species freeze.

Detailed provenance and units are fixed in `docs/PIPELINE_34SPECIES.md`.

## Layer 2 — Spatial classification

The umbrella term is **intraspecific flower-colour variation**.

- `within_population`: at least two discrete natural colour variants are explicitly documented as coexisting within at least one population.
- `among_population`: geographic or among-population differentiation is documented without retained evidence of local coexistence.
- `mixed`: both signals occur.
- `unclear`: retained evidence does not resolve the spatial organization.

Mixed and unclear species are not forced into the binary response. The active deterministic rules are centralized in `fcp_pipeline/evidence.py`.

The current labels are **source-traceable, rule-derived classifications**. The repository contains a blinded review protocol, but completed independent human review must not be claimed unless completed reviewer files actually exist.

## Layer 3 — Frozen comparative dataset

The final binary analysis is hard-frozen at:

- 34 species;
- 25 plant families;
- 20 within-population cases;
- 14 among-population cases;
- at least 20 occupied climate cells per species;
- `classification_source = baseline_unambiguous` for every row.

The frozen biological sample is validated independently of the climatic model results. A change in any of these counts requires a new explicit freeze rather than a silent update.

## Layer 4 — GBIF occurrence and occupied-climate construction

Coordinate-bearing GBIF records are linked to WorldClim 2.1 bioclimatic values. Within a species, records with identical nine-variable climate vectors are deduplicated to occupied climate cells.

The resulting variables describe **realised occupied climate**. They do not represent physiological tolerance, fundamental niche, morph-specific climate, gene flow or local selective environments.

Five species-level climatic-niche summaries are evaluated symmetrically:

1. temperature breadth;
2. moisture breadth;
3. climatic heterogeneity;
4. PCA dispersion;
5. PCA hull area.

## Layer 5 — Primary five-metric models

Every climatic metric uses the same response, species and covariate structure:

`among ~ metric_z + effort_z`

where:

- `among = 1` for geographically structured cases and `0` for within-population cases;
- `metric_z` is the standardized focal climatic metric;
- `effort_z = z(log1p(n_climate_cells))`.

The five metrics are fitted **separately**, not entered jointly. Primary uncertainty uses family-clustered sandwich covariance.

The production implementation is `scripts/run_34species_models.py`, with reusable model logic in `fcp_pipeline/models.py`.

## Layer 6 — Robustness and multiplicity

For every one of the five primary models:

- 9,999 spatial-label permutations are performed;
- leave-one-family-out refits quantify concentration in a represented family;
- Holm-adjusted Wald and permutation p-values are reported across the five metrics;
- predictor correlation, VIF and exact condition number diagnose collinearity with occurrence effort.

Multiplicity is treated as interpretive context. The paper emphasizes effect sizes, confidence intervals and directional stability rather than defining discovery by one p-value threshold.

## Layer 7 — Phylogenetic sensitivity

Two complementary phylogenetic analyses use the same model structure:

### Open Tree topology

- 30 uniquely matched species;
- 100 random polytomy resolutions;
- Grafen branch lengths;
- `phylolm::phyloglm(method = "logistic_MPLE")`.

### Time-scaled V.PhyloMaker2 trees

- all 34 species;
- `GBOTB.extended.LCVP` backbone;
- placement scenarios S1, S2 and S3;
- same `among ~ metric_z + effort_z` formula.

These analyses are phylogenetic **sensitivity analyses**, not a claim that the ecological process itself follows the assumed tree model.

Production entry point: `scripts/run_34species_phylogenetic.R`.

## Layer 8 — Finite-sample diagnostics

Because the frozen sample contains 34 species in 25 families, two additional diagnostics quantify rather than conceal finite-sample uncertainty:

1. CR2 family-cluster covariance with Satterthwaite degrees of freedom (`scripts/run_34species_cr2.R`);
2. design-based power/precision simulation retaining the observed 34-species predictor and effort distributions (`scripts/run_34species_power_precision.py`).

These diagnostics are not used post hoc to declare the sample “adequate.” They show how much directional and threshold-significance information the actual design contains.

## Layer 9 — Result hierarchy

The current frozen five-metric point estimates are all negative. Moisture breadth shows the largest observed contrast (OR about 0.41), but:

- Holm adjustment weakens threshold-based support;
- CR2/Satterthwaite expands uncertainty;
- OpenTree and dated phylogenetic confidence intervals include one.

The paper therefore makes the following bounded inference:

> geographically structured flower-colour variation tends to occur toward the narrower end of sampled occupied climatic niche breadth than within-population coexistence, with the largest observed contrast along the moisture axis.

It does **not** claim that moisture uniquely causes geographic differentiation.

## Interpretation rules

- Use **association**, not causation.
- Do not infer local adaptation without morph-labelled locality evidence.
- Do not equate species-level climatic breadth with morph-specific tolerance.
- Do not call geographically separated variants within-population flower-colour polymorphism.
- Do not treat non-significance as proof of absence.
- Do not describe the 34 species as a random sample of all angiosperms.
- Do not call the study a conventional effect-size meta-analysis; it is a literature-derived cross-species comparative evidence synthesis.

## Canonical reproduction

The active paper analysis is reproduced by `.github/workflows/34species-paper.yml`, which hard-fails if the frozen sample counts, five-metric symmetry or 9,999-permutation requirement are violated.
