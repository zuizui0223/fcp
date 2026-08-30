# JBI Chapter 1 spatial execution status

## Current decision

Branch: `analysis/jbi-global-colour-boundaries`  
PR: #20, stacked on `audit/upstream-spatial-reaudit` / PR #17.

As of **2026-08-30**, the frozen 720-photograph evaluation has been completed and the two confirmatory spatial stages have been run.

The result is asymmetric:

1. **Stage A is supported.** Continuous flower-colour observations are more locally organized within species than expected under species-conditioned random labelling.
2. **Stage B is not supported by the frozen primary analysis.** The six species do not show confirmatory evidence that their strongest continuous colour transitions concentrate in the same global cells.

The chapter should therefore be framed around **repeated within-species spatial organization without a demonstrated universal shared boundary**, not around a single global transition line.

## Frozen empirical sample

- acquired photographs: **1,200 = 6 species × 200**;
- calibration split: **480 = 80/species**;
- held-out evaluation split: **720 = 120/species**;
- split basis: deterministic, outcome-blind and hash-frozen;
- evaluation records are joined to the frozen coordinates by `photo_id`;
- no colour vector is ever permuted between species.

## Frozen colour representation

The discrete-morph route was not forced where calibration geometry did not support it. Before the evaluation set was opened, the operational representation was frozen as a **species-specific continuous colour vector** derived from Florence-localized flower regions and standardized using calibration-only means and population standard deviations.

The representation is stored in:

- `docs/supporting/jbi_ch1_continuous_colour_representation_v1.json`;
- `docs/supporting/jbi_ch1_evaluation_opening_contract_v1.json`.

Reviewer-2 scoring is no longer a gate for this analysis. It remains optional calibration-QC infrastructure and is not allowed to retroactively change the completed confirmatory representation or evaluation set.

## 720-photo evaluation — complete

Workflow run: `33281907575`.

- deterministic extraction shards: **36/36 successful**;
- evaluation feature records: **720/720**;
- unique `photo_id`: **720/720**;
- per species: **120/120**;
- `feature_status = ok`: **720/720**;
- localization failures: **0**;
- `final_label = false`: **720/720**;
- post-opening rule tuning: **none**.

Committed products:

- `data/evaluation/jbi_ch1_florence_evaluation_features_v1.jsonl`;
- `docs/supporting/jbi_ch1_florence_evaluation_features_v1.json`.

## Stage A — within-species local organization

### Frozen test

For each species:

1. construct a colour-blind spherical k-nearest-neighbour graph from observation coordinates;
2. calculate RMS Euclidean discontinuity between complete standardized colour vectors on every retained edge;
3. keep graph geometry fixed;
4. permute complete colour vectors strictly within species;
5. recompute species-specific discontinuity `Q_i`;
6. average the six `Q_i` values with equal species weight.

Primary geometry: `k = 5`.  
Predeclared sensitivities: `k = 3` and `k = 8`.  
Permutations: **9,999** for each analysis.

Workflow run: `33283136767`.

### Global result

| k | observed equal-species Q | null mean | standardized clustering deficit | lower-tail p |
|---:|---:|---:|---:|---:|
| 3 | — | — | 2.5523 | 0.0066 |
| **5 primary** | **1.39114** | **1.42943** | **2.3113** | **0.0113** |
| 8 | — | — | 2.6133 | 0.0065 |

The lower observed discontinuity means geographically neighbouring photographs are more similar in continuous flower colour than expected after conditioning on each species' locations and colour-vector distribution.

### Species-specific primary results

| Species | lower-tail p | interpretation at species level |
|---|---:|---|
| *Antirrhinum majus* | 0.2635 | not individually resolved |
| *Dactylorhiza sambucina* | 0.2560 | not individually resolved |
| *Gentiana lutea* | 0.7865 | not individually resolved |
| *Ipomoea purpurea* | 0.2923 | not individually resolved |
| *Lysimachia arvensis* | 0.0023 | strong local organization |
| *Raphanus sativus* | 0.0080 | strong local organization |

The confirmatory claim is the equal-species global result. The species table identifies heterogeneous contributors; it does not redefine the global test after inspection.

Committed products:

- `docs/supporting/jbi_ch1_stage_a_continuous_graph_contract_v1.json`;
- `docs/supporting/jbi_ch1_stage_a_continuous_graph_v1.json`;
- `data/evaluation/jbi_ch1_stage_a_primary_null_v1.csv`.

## Stage B — cross-species shared-transition concentration

### Frozen geometry selection

Before observed colour discontinuities were scored, the primary spatial support was selected using species identities, coordinates and graph distances only.

Candidate priority:

- edge caps: 500, 1,000, 2,000 km;
- equal-area longitude–sin(latitude) grids: 36×18, 24×12, 18×9;
- species/cell detectable only with at least two retained geometry edges;
- shared cell evaluable only when at least two species are detectable.

All nine candidate configurations passed the frozen geometry-support criteria. The first passing configuration was therefore the primary:

**500-km edge cap + 36×18 equal-area grid**.

Primary support:

- retained edges by species: 396, 416, 396, 387, 397 and 368;
- detectable cells by species: 40, 38, 37, 49, 34 and 2;
- cells with `A ≥ 2`: 24;
- cells with `A ≥ 3`: 10;
- cells with `A ≥ 4`: 5;
- maximum opportunity `A`: 5 species.

### Frozen concentration test

Within each species and configuration, edge discontinuities were average-rank transformed to transition intensities in `[0,1]`. Species-cell intensities were averaged only across detectable species. Cells with insufficient opportunity were `NaN`, never zero.

The primary statistic was the opportunity-weighted variance of shared transition intensity across cells with `A ≥ 2`. Its complete null pipeline used 9,999 within-species vector permutations and recomputed edge scores, ranks, cell intensities, the shared surface and the concentration statistic each time.

Workflow run: `33284194283`.

### Primary result

- observed concentration: **0.0060573**;
- null mean: **0.0045443**;
- null SD: **0.0011504**;
- standardized concentration excess: **1.3153**;
- Monte Carlo upper-tail p: **0.0906**;
- descriptive two-sided p: **0.1708**;
- reject shared-concentration null at 0.05: **no**.

### Predeclared sensitivity configurations

| Edge cap | Grid | upper-tail p |
|---:|---:|---:|
| 500 km | 24×12 | 0.0445 |
| 500 km | 18×9 | 0.3415 |
| 1,000 km | 36×18 | 0.2235 |
| 1,000 km | 24×12 | 0.4945 |
| 1,000 km | 18×9 | 0.1920 |
| 2,000 km | 36×18 | 0.4690 |
| 2,000 km | 24×12 | 0.5500 |
| 2,000 km | 18×9 | 0.3495 |

The nominal `p = 0.0445` result for the coarser 500-km/24×12 grid is exploratory sensitivity evidence only. It does not replace the prospectively selected primary configuration, and the sensitivity set is not uniformly supportive.

Committed products:

- `docs/supporting/jbi_ch1_stage_b_shared_transition_contract_v1.json`;
- `docs/supporting/jbi_ch1_stage_b_geometry_audit_v1.json`;
- `docs/supporting/jbi_ch1_stage_b_shared_transition_concentration_v1.json`;
- `data/evaluation/jbi_ch1_stage_b_shared_transition_surface_v1.csv`;
- `data/evaluation/jbi_ch1_stage_b_primary_null_v1.csv`.

## Claim boundary

Supported:

> Across six independently sampled species, continuous flower colour is spatially organized within species more strongly than expected under species-conditioned random labelling.

Not supported as a confirmatory claim:

> Independent species repeatedly place their strongest flower-colour transitions along one common global geographic boundary.

Not tested by the completed stages:

- a climatic cause;
- a historical cause;
- identical biological mechanisms among species;
- morph-specific fitness or adaptation.

Because Stage B did not pass its primary gate, a post-discovery geographic-reference correspondence analysis is **not promoted to the confirmatory main line**. Running such a correspondence now would be explicitly exploratory and cannot rescue the Stage B hypothesis.

## Active next work

1. generate the Chapter 1 result figures and manuscript-ready tables from the committed Stage A/B outputs;
2. revise the chapter narrative from “shared global boundary” to “repeated local organization with scale-dependent but non-confirmatory cross-species overlap”;
3. retain the full Stage B surface and sensitivity audit in Supporting Information;
4. keep reviewer-2 and geographic-cause development separated as optional later QC/exploration rather than reopening the frozen evaluation.
