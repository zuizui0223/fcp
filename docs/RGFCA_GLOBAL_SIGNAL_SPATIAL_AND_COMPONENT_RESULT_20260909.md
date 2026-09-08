# RGFCA recurrent signal geography and flower/background decomposition — exploratory result

Date: 2026-09-09 JST

Upstream frozen protocols/results:

- signal recurrence protocol: `958eb931b07d768bef41a6fe08a770bf8cef1793`;
- spatial recurrence protocol: `b42a8625f7cd59837301f1d1070cfb73fac103d7`;
- flower/background decomposition protocol: `36f2d1b6d7af31e683adcad93ebf436ea20799e7`.

## Spatial recurrence

The fixed 500-km reporting lattice and >=10-species/realization support rule yield **38 main-map cells** that are supported in >=80% of the 200 realizations. Coverage is not globally uniform: these 38 cells fall in the currently well-supported North American, European, and Australasian parts of the reserve. The map must therefore be read as a support-limited recurrence atlas, not a globally exhaustive map.

| Mode | main-map cells | R_sign >=0.80 | R_sign >=0.90 | cells with observed-minus-null R_sign >=0.10 | maximum R_sign | maximum excess over null |
|---|---:|---:|---:|---:|---:|---:|
| M1 | 38 | 10 | 2 | 27 | 1.000 | 0.465 |
| M2 | 38 | 10 | 3 | 24 | 0.940 | 0.430 |
| M3 | 38 | 11 | 5 | 30 | 0.965 | 0.430 |
| M4 | 38 | 9 | 0 | 29 | 0.890 | 0.365 |
| M5 | 38 | 7 | 3 | 20 | 0.979 | 0.440 |

The species-conditioned coordinate-permutation null is usually near chance sign recurrence in the strongest cells, while observed recurrence can approach 1.0. Thus the retained signal modes are not merely stable in colour space; within the covered regions, their **signs recur geographically**.

Examples among M1–M3 include an M1-positive cluster in eastern Australia, M2-positive cells in the southwestern/southeastern United States, and strong M3 sign structure across western and southeastern North America. These location descriptions are descriptive labels on fixed reporting cells, not ecological explanations.

A separate post-hoc adjacency scan (not a prespecified discovery criterion) finds edge-neighbouring, opposite-sign strong cells in several places, including M3 along western North America and the southeastern United States. These are candidate transition zones only; they are not yet biological boundaries.

## Flower/background decomposition

For every retained state and realization, the matched signal was decomposed exactly as

`S_D = S_F - S_B`,

where D is flower-minus-background, F is flower-only and B is matched-background, all projected onto the already frozen recurrent mode vector after within-species centering. Maximum numerical identity errors are ~1e-16.

### Global component diagnostic

| Mode | median corr(D,F) | median corr(D,-B) | median mean |D| | median mean |F| | median mean |B| | fraction |F|>|B| | fraction |B|>|F| |
|---|---:|---:|---:|---:|---:|---:|---:|
| M1 | 0.584 | 0.799 | 0.1190 | 0.0661 | 0.0973 | 0.368 | 0.632 |
| M2 | 0.702 | 0.586 | 0.0974 | 0.0776 | 0.0725 | 0.503 | 0.497 |
| M3 | 0.470 | 0.861 | 0.0758 | 0.0364 | 0.0666 | 0.315 | 0.685 |

Therefore **M1 and M3 are globally background-heavy in this photographic representation**, while M2 is approximately balanced/slightly more flower-associated. This prevents a blanket interpretation of recurrent D modes as floral adaptation or pollinator signalling.

## Flower-side candidate cells

For transparent exploratory triage only, a post-result candidate table is retained using the descriptive rule:

- retained main-map M1–M3 cell;
- D `R_sign >= 0.80`;
- D observed-minus-null `R_sign >= 0.20`;
- flower-only `R_sign >= 0.80`;
- flower-only recurrence greater than background-only recurrence.

This leaves **13 candidate cells**. The rule is explicitly post hoc and is not a significance or confirmation threshold.

Notable flower-side candidates include:

- M1: Pacific Northwest / western Canada cells and a southeastern-Florida cell;
- M2: Florida, central North America, southeastern Australia, and western-Canada cells;
- M3: adjacent western-North-America cells with opposite signs, a southwestern-US cell, and eastern Australia.

These are the first locations where the recurrent matched signal remains strongly recurrent on the **flower side itself**, but even these remain photo-derived colour structure. They do not yet establish animal-perceived contrast, attraction, adaptation, or a causal environmental driver.

## Verification

- The spatial analysis reproduced all upstream 200-realization mode statistics to floating precision (maximum difference `4.44e-16`).
- The retained D spatial recurrence in the component diagnostic reproduced the prior spatial result exactly (`max difference = 0`).
- Re-running the spatial analysis produced byte-identical science outputs; the only omitted item was a separate lon/lat convenience table generated after the main run.

## Next interpretation boundary

The next stage may overlay ecological covariates on these **already frozen modes/cells**. Island/mainland status, climate, latitude, realm, pollinator guild, reproductive system, and taxonomic composition were not used to construct the recurrent modes or select the spatial cells above. Any ecological association remains exploratory until independently replicated.
