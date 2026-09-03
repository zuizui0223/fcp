# JBI Chapter 1 spatial execution status

## Final Chapter-1 decision

The held-out six-species spatial analysis is complete and frozen.

1. **Stage A supported:** continuous flower colour is more locally organized within species than expected under species-conditioned random labelling.
2. **Stage B not confirmed by the frozen primary analysis:** independent species did not show confirmatory concentration of their strongest continuous colour transitions in the same global cells.

The manuscript therefore centers on **repeated within-species spatial organization without a demonstrated universal shared boundary**.

The separate prospective 200-species / 60,000-photo terminal scale-out has also ended. It returned **`not_evaluable_scaleout_measurement_completeness` before species-conditioned spatial organization**, so it does not alter or replicate the six-species biological result. See `docs/JBI_CHAPTER1_TERMINAL_SCALEUP_RESULT.md` and `docs/JBI_IMAGE_FIRST_ATLAS_STATUS.md`.

## Frozen empirical sample

- acquired photographs: **1,200 = 6 species × 200**;
- calibration split: **480 = 80/species**;
- held-out evaluation split: **720 = 120/species**;
- deterministic, outcome-blind, hash-frozen split;
- evaluation records joined to frozen coordinates by `photo_id`;
- complete colour vectors permuted only within species.

The operational phenotype was frozen before evaluation as a **species-specific continuous colour vector** from Florence-localized flower regions, standardized using calibration-only means and population standard deviations. No evaluation-derived feature selection, scaling rule or localization rule was introduced after opening.

## Evaluation integrity

Workflow run `33281907575` completed all 36 deterministic extraction shards:

- evaluation features: **720/720**;
- per species: **120/120**;
- `feature_status = ok`: **720/720**;
- localization failures: **0**;
- `final_label = false`: **720/720**;
- post-opening rule tuning: **none**.

## Stage A — supported local organization

The frozen test used colour-blind spherical within-species nearest-neighbour graphs, RMS discontinuity between complete standardized vectors, an equal-species global statistic and 9,999 complete-vector permutations strictly within species.

| k | observed equal-species Q | null mean | standardized clustering deficit | lower-tail p |
|---:|---:|---:|---:|---:|
| 3 | 1.37661 | 1.42967 | 2.5523 | 0.0066 |
| **5 primary** | **1.39114** | **1.42943** | **2.3113** | **0.0113** |
| 8 | 1.39556 | 1.42978 | 2.6133 | 0.0065 |

Primary species-level lower-tail probabilities were 0.2635 (*Antirrhinum majus*), 0.2560 (*Dactylorhiza sambucina*), 0.7865 (*Gentiana lutea*), 0.2923 (*Ipomoea purpurea*), 0.0023 (*Lysimachia arvensis*) and 0.0080 (*Raphanus sativus*). These describe heterogeneity; the confirmatory claim is the predeclared equal-species global test.

## Stage B — primary shared concentration not confirmed

Geometry-only selection chose the first passing frozen support: **500-km edge cap + 36×18 equal-area grid**. Twenty-four cells had opportunity `A ≥ 2`, nine had `A ≥ 3`, four had `A ≥ 4`, and maximum opportunity was all six species.

Primary result:

- observed concentration: **0.0082315**;
- null mean: **0.0056757**;
- null SD: **0.0017763**;
- standardized concentration excess: **1.4389**;
- upper-tail Monte Carlo `p = 0.0906`;
- descriptive two-sided `p = 0.1372`.

One predeclared coarser sensitivity, 500-km / 24×12, had `p = 0.0445`; the other seven sensitivity supports were not nominally supported. That isolated sensitivity cannot replace the frozen primary support.

## Terminal scale-out — measurement NE, not biological absence

The later prospective atlas extension used eight disjoint 25-species cohorts with 300 fixed observations per species: **200 species / 60,000 photographs**. The exact terminal run `33592829701` completed the frozen location-blind measurement/reassembly layer as **256 partitions -> 16 semantic shards -> 60,000 unique terminal records**.

The frozen measurement gate then found only **58 / 200 measurement-evaluable species**, with cohort counts **3, 4, 7, 8, 8, 12, 7, 9**. All eight cohorts were `not_evaluable`. The binding state was:

- `coordinate_join_permitted = false`;
- `coordinates_opened = false`.

Accordingly, the terminal cascade stopped before species-conditioned spatial organization. Shared-transition, environmental-concordance and pollinator-concordance branches were not run. No threshold, seed, denominator, estimator, source, metric, cohort requirement or branch order was changed after the result.

This scale-out is a **measurement-transfer/evaluability result**, not evidence against spatial organization or any downstream biological hypothesis. The 58 evaluable species cannot be promoted as a favourable post-result subset.

## Claim boundary

Supported:

> Across the six-species held-out sample, continuous flower colour is spatially organized within species more strongly than expected under species-conditioned random labelling.

Not supported as a confirmatory claim:

> Independent species repeatedly place their strongest flower-colour transitions along one common global geographic boundary.

Not evaluated by the terminal 200-species experiment:

- species-conditioned spatial organization;
- shared-transition concordance;
- climatic or historical concordance;
- pollinator biogeographic concordance.

## Integration state

The Chapter-1 result figures, numerical tables and manuscript narrative are already generated and versioned. The terminal scale-out has been incorporated into the manuscript only as a prospectively stopped measurement-scaling extension and limitation; it does not change the Abstract's primary six-species biological conclusion.

Current manuscript/reproducibility entry points:

- `docs/JBI_CHAPTER1_MANUSCRIPT.md`;
- `docs/JBI_CHAPTER1_RESULTS.md`;
- `docs/JBI_CHAPTER1_FIGURE_PLAN.md`;
- `docs/JBI_CHAPTER1_TERMINAL_SCALEUP_RESULT.md`;
- `docs/JBI_IMAGE_FIRST_ATLAS_STATUS.md`;
- `docs/supporting/jbi_atlas_terminal_measurement_v5_receipt.json`.

No additional confirmatory coordinate-opening branch is authorized from the terminal experiment. The next work is submission packaging, reference/claim audit and reproducibility cleanup rather than further atlas inference.
