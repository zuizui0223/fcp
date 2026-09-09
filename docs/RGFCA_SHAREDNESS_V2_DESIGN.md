# RGFCA sharedness v2 — current prospective design ladder

Updated 8 September 2026. This document records the current forward design only. It does not reopen or reclassify the completed discovery/reserve outcomes.

## Why v2 exists

The completed RGFCA work established two distinct limitations.

1. The original G1 statistic is not a sharedness-specific statistic. In synthetic worlds it can reject when every species has its own independent spatial boundary, so G1 cannot support a cross-species common-boundary claim.
2. The species-disjoint sharedness-specific predictive statistic fixes that estimand problem, but its qualification on the old geometry failed. The old implementation used 100 training species plus 100 disjoint evaluation species and only 20 photographs per species. It had low power even when synthetic shared boundaries were deliberately implanted.

The v2 response is therefore **not** to relax the sharedness question to generic within-species spatial structure. The response is to prospectively change the data design so that each species is sampled deeply enough to span its own geography before any new colour outcome is opened.

The independent reserve result remains unchanged: the weak distance–flower association replicated, while the prespecified matched flower-minus-background differential did not pass (`p = 0.087`; `flower_specific_robust_replication=false`).

## Gate 1 — metadata capacity

The complete metadata-only parent census contains 42,111 species and no remaining request errors after the frozen 429 transport recovery. Its original request ceiling was 200 observations per species, so rows at `raw_results == 200` are right-censored for higher-depth design questions.

A frozen re-summary at `docs/supporting/rgfca_sharedness_v2_capacity_summary_v1.json` gives observer-capped capacities:

| minimum metadata-eligible photos after observer cap | species |
| ---: | ---: |
| 100 | 4,730 |
| 125 | 3,934 |
| 150 | 3,082 |
| 175 | 1,727 |
| 200 | 21 |

There are 3,440 right-censored species at the old 200-result API ceiling. Of those, 2,866 already retain at least 150 rows and 1,700 retain at least 175 rows after the frozen observer cap. These counts establish substantial high-depth potential but **cannot** establish 300- or 400-photo availability.

A prospective high-depth metadata-only pilot is therefore frozen and running under `rgfca-sharedness-v2-high-depth-capacity-pilot-v1`:

- fixed pool: the 3,440 right-censored species;
- fixed pilot: 500 species chosen by SHA-256 rank of iNaturalist taxon ID, without replacement;
- same research-grade / flowering / geolocation / positional-accuracy / licence / prior-ID-exclusion rules as the parent census;
- observer cap = 2;
- deterministic pagination, at most five 200-row pages per species;
- fixed targets: 300 and 400 retained metadata rows after observer cap;
- primary feasibility PASS: at least 150 of the fixed 500 species reach 300;
- request-error ceiling: 5%; no replacement species;
- page-cap misses are indeterminate, not declared biologically or observationally insufficient.

No image pixel, flower colour, climate value, G1/G3 outcome, or reserve outcome is used by this gate.

## Gate 2 — measurement representation

The v2 measurement output is already frozen before any v2 candidate pixel opening in `rgfca-sharedness-v2-continuous-visible-colour-measurement-v1`.

The design does not invent a new colour extractor. It preserves the continuous visible-colour summaries that the exact frozen ROI-v4 runtime already computes before palette compression:

- flower and matched-background CIELAB `L`, `a`, `b` mean;
- per-channel SD and q10/q50/q90 retained as diagnostics;
- the same horizontal-flip ROI and colour-stability checks;
- the old 12-anchor palette counts retained only as a compatibility/secondary representation.

The future continuous sharedness primary representation is fixed as the three-dimensional visible `L_mean, a_mean, b_mean` vector. Ordinary RGB photographs cannot support UV-bullseye or pollinator-vision claims.

Matched background continuous summaries are mandatory so that a future observed sharedness test cannot ignore the photography/background control that stopped the reserve flower-specific claim.

## Gate 3 — sharedness-specific synthetic identifiability

If and only if Gate 1 passes, the next scientific gate is a new synthetic qualification before any v2 colour outcome is opened.

The target estimand remains cross-species shared geographic structure, not generic within-species autocorrelation. The v1 species-disjoint predictive principle is retained: learn the shared structure from one species set and score genuinely disjoint species against a structured null in which each species may possess its own strong spatial boundary.

The v2 qualification must differ from v1 in the data geometry and colour representation, not by searching for a favourable test statistic after the fact:

- use the prospectively selected high-depth species/photo geometry;
- use species-disjoint train/evaluation splits;
- preserve independent species-specific boundaries in the null;
- simulate continuous visible CIELAB-like responses, with matched nuisance/background structure, rather than relying only on the old binary synthetic label;
- calibrate rejection from structured nuisance worlds before evaluating positive worlds;
- freeze hard power/error gates before the synthetic outcomes are generated;
- do not use observed v2 flower colour in method selection.

A synthetic PASS does not itself authorize an ecological claim. It only permits the later, separately frozen new-image measurement/inference stage.

## Pollination guilds

Pollination guild remains a possible **prospective secondary stratification**, not a current primary gate. Existing FCP literature shortlists are unsuitable because they are outcome-/topic-selected review assets rather than an external complete species trait frame. BIOLFLOR and recent plant–pollinator network compilations provide useful pollination information but are geographically concentrated, particularly in Europe. A global guild stratum will be added only if a separately audited external source has adequate taxonomic and geographic coverage; otherwise it will not be used to define the primary v2 species sample.

## Hard boundaries

- The six-species and 34-species analyses are not part of this main line.
- The old 369-species threshold is not a v2 biological target.
- G1 `p=.070` is not promoted as sharedness evidence.
- The reserve differential `p=.087` is not rescued or reclassified.
- No threshold may be selected by opening several v2 colour outcomes and choosing the favourable one.
- No new image acquisition or v2 pixel opening occurs before metadata capacity and synthetic identifiability gates are resolved.
