# JBI Chapter 1 — manuscript-ready spatial results

## Answer first

The held-out six-species analysis supports **within-species spatial organization of continuous flower colour**, but it does not provide confirmatory evidence for **one shared global geography of the strongest colour transitions**.

This yields a two-level result:

- neighbouring observations are more colour-similar than expected after conditioning on species distributions;
- the locations of relatively strong transitions are not consistently concentrated in the same global cells across the six species at the frozen primary scale.

## Study sample and evaluation integrity

The frozen design contained 1,200 community photographs from six species, with 80 calibration and 120 evaluation photographs per species. The colour representation was fixed before evaluation as a species-specific standardized continuous vector derived from the frozen Florence-localized flower region and calibration-only scaling parameters.

All 720 held-out photographs were processed successfully. Each species contributed 120 unique evaluation records; all 720 records had `feature_status = ok`, no localization failure occurred, and no final discrete biological label was emitted. The evaluation path was therefore complete without post-opening rule modification.

## Stage A: continuous colour is locally organized within species

### Confirmatory analysis

Within each species, a spherical five-nearest-neighbour graph was constructed from coordinates without using flower-colour values. RMS Euclidean discontinuity was calculated between standardized continuous colour vectors on graph edges. Complete vectors were then permuted strictly within species while graph geometry remained fixed. Species-specific mean edge discontinuities were calculated first and averaged with equal species weight. The primary test used 9,999 permutations; `k = 3` and `k = 8` were predeclared sensitivities.

### Result

The primary equal-species mean discontinuity was **1.39114**, below the null mean of **1.42943**. The standardized clustering deficit was **2.311**, and the lower-tail Monte Carlo p-value was **0.0113**. Thus geographically neighbouring observations were more similar in continuous flower colour than expected under species-conditioned random labelling.

The result was stable to graph scale:

| Graph k | Observed global Q | Null mean | Standardized clustering deficit | Lower-tail p |
|---:|---:|---:|---:|---:|
| 3 | 1.37661 | 1.42967 | 2.552 | 0.0066 |
| **5, primary** | **1.39114** | **1.42943** | **2.311** | **0.0113** |
| 8 | 1.39556 | 1.42978 | 2.613 | 0.0065 |

The primary species-specific tests were heterogeneous:

| Species | Lower-tail p | Stage-A reading |
|---|---:|---|
| *Antirrhinum majus* | 0.2635 | no individually resolved clustering |
| *Dactylorhiza sambucina* | 0.2560 | no individually resolved clustering |
| *Gentiana lutea* | 0.7865 | no individually resolved clustering |
| *Ipomoea purpurea* | 0.2923 | no individually resolved clustering |
| *Lysimachia arvensis* | 0.0023 | strong local organization |
| *Raphanus sativus* | 0.0080 | strong local organization |

The confirmatory inference is the prospectively defined equal-species global test, not a post-hoc count of individually significant species. The species results show that the common global direction is heterogeneous and is most strongly expressed by *Lysimachia arvensis* and *Raphanus sativus* in the primary graph.

## Stage B: no confirmatory common global transition concentration

### Confirmatory analysis

Because Stage A passed its gate, the analysis proceeded to cross-species transition concentration. A spatial configuration was selected before any observed colour edge score was computed, using only species identities, coordinates and graph distances. The first configuration satisfying the frozen support criteria was a **500-km edge cap with a 36×18 equal-area longitude–sin(latitude) grid**.

A species/cell was detectable only when at least two retained geometry edges had midpoints in that cell. Shared cells required at least two detectable species. For every species, retained edge discontinuities were average-rank transformed to `[0,1]`, then averaged within detectable cells. The shared intensity `S(x)` was the mean of species-cell intensities among detectable species, with `A(x)` retained as the explicit opportunity denominator. Cells with insufficient opportunity were not evaluable rather than assigned zero.

The primary statistic was the `A(x)`-weighted variance of `S(x)`. The complete null recomputed edge discontinuities, within-species ranks, species-cell intensities, the shared surface and concentration after each of 9,999 within-species complete-vector permutations.

### Geometry support

The primary configuration retained **246, 363, 372, 241, 300 and 329** edges for *Antirrhinum*, *Dactylorhiza*, *Gentiana*, *Ipomoea*, *Lysimachia* and *Raphanus*, respectively. Detectable-cell counts were **15, 8, 7, 17, 22 and 6**. Across the global grid, 24 cells had opportunity `A ≥ 2`, nine had `A ≥ 3`, four had `A ≥ 4`, and the maximum opportunity was all six species.

### Result

Observed opportunity-weighted concentration was **0.0082315**, compared with a null mean of **0.0056757** and null SD of **0.0017763**. The standardized concentration excess was **1.4389**. The Monte Carlo upper-tail p-value was **0.0906**; the descriptive two-sided p-value was **0.1372**.

The frozen primary shared-concentration null was therefore **not rejected**.

Predeclared sensitivity analyses were not uniformly supportive:

| Edge cap | Grid | Upper-tail p |
|---:|---:|---:|
| 500 km | 24×12 | 0.0445 |
| 500 km | 18×9 | 0.3415 |
| 1,000 km | 36×18 | 0.2235 |
| 1,000 km | 24×12 | 0.4945 |
| 1,000 km | 18×9 | 0.1920 |
| 2,000 km | 36×18 | 0.4690 |
| 2,000 km | 24×12 | 0.5500 |
| 2,000 km | 18×9 | 0.3495 |

The nominal result at 500 km and 24×12 (`p = 0.0445`) is retained as exploratory scale sensitivity. It cannot replace the first prospectively selected passing configuration, and its lack of consistency across neighbouring supports prevents it from establishing a robust shared global boundary.

## Manuscript result paragraphs

### Results paragraph 1 — held-out measurement

The frozen evaluation workflow produced valid continuous colour measurements for all 720 held-out photographs (120 per species), with no localization failure and no post-opening modification of the representation or scaling rules. Each observation retained a species-specific standardized continuous colour vector rather than a forced universal discrete-morph label.

### Results paragraph 2 — Stage A

Continuous flower colour was more locally organized within species than expected under species-conditioned random labelling. In the primary five-nearest-neighbour graph, the observed equal-species mean edge discontinuity was 1.391, compared with a permutation-null mean of 1.429 (standardized clustering deficit = 2.31; 9,999-permutation lower-tail `p = 0.0113`). The direction was retained at `k = 3` (`p = 0.0066`) and `k = 8` (`p = 0.0065`). Species-specific support was heterogeneous and strongest for *Lysimachia arvensis* (`p = 0.0023`) and *Raphanus sativus* (`p = 0.0080`).

### Results paragraph 3 — Stage B

The subsequent shared-transition test did not reject its frozen primary null. Geometry-only selection chose a 500-km edge cap and 36×18 equal-area grid, yielding 24 globally evaluable cells with at least two detectable species. Opportunity-weighted concentration of shared transition intensity was 0.00823, compared with a null mean of 0.00568 (standardized excess = 1.44; 9,999-permutation upper-tail `p = 0.0906`). One coarser 500-km sensitivity grid was nominally supported (`p = 0.0445`), whereas the other seven sensitivity configurations were not, indicating scale-dependent exploratory overlap rather than a confirmatory universal boundary.

## Discussion-ready interpretation

The two stages distinguish spatial organization from spatial coincidence. Stage A shows that community photographs retain enough structured colour information to detect non-random within-species geographic organization after preserving each species' observed range and sampling geometry. Stage B then asks the stronger question of whether independent species place their strongest relative transitions in the same global cells. Failure to reject that primary null means local organization should not be collapsed into a claim of one shared biogeographic boundary.

The supported biological statement is therefore:

> Independent species can show spatially organized continuous flower-colour variation, but the geography of that organization is heterogeneous rather than demonstrably governed by one universal global boundary in the present six-species sample.

This result is compatible with species-specific combinations of environmental gradients, demographic history, dispersal, biotic interactions and image-scale measurement variation. Those alternatives are not distinguished by the present analysis.

## Claims excluded from the current chapter

The completed analysis does not establish:

- a particular climatic driver;
- a particular historical boundary;
- identical selection mechanisms among species;
- morph-specific fitness differences;
- physiological colour measurements from calibrated spectra.

Because the Stage-B confirmatory gate was not passed, a geographic-reference correspondence search is not advanced as the next confirmatory test. Any later environmental or historical overlay must be described as exploratory and kept separate from the frozen Stage-A and Stage-B claims.

## Versioned evidence

Primary source products:

- `data/evaluation/jbi_ch1_florence_evaluation_features_v1.jsonl`;
- `docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json`;
- `data/evaluation/jbi_ch1_stage_a_primary_null_v1.csv`;
- `docs/supporting/jbi_ch1_stage_b_geometry_audit_v1.json`;
- `docs/supporting/jbi_ch1_stage_b_shared_transition_concentration_v1.json`;
- `data/evaluation/jbi_ch1_stage_b_shared_transition_surface_v1.csv`;
- `data/evaluation/jbi_ch1_stage_b_primary_null_v1.csv`.
