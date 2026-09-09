# RGFCA global flower–background signal recurrence — exploratory result

Date: 2026-09-09 JST

Protocol frozen before recurrence computation: `958eb931b07d768bef41a6fe08a770bf8cef1793`.

## Scope

This result is the hypothesis-free recurrence stage defined in `RGFCA_GLOBAL_SIGNAL_RECURRENCE_EXPLORATION_PROTOCOL_20260909.md`. It does not prespecify an island syndrome, whitening, a universal colour boundary, pollinator guild, climate direction, or ecological covariate. The measured quantity is **photo-derived flower-minus-matched-background 12-anchor palette structure**. It is not calibrated animal-perceived colour or demonstrated attraction.

## Frozen execution

- Eligible source: **24,885 photos / 500 species** from the frozen RGFCA reserve.
- Local states: 300-km EPSG:6933 shifted grid; >=3 eligible photos; >=2 unique observers.
- A species enters a realization only with >=2 retained local states separated by >=500 km.
- Local-state means use observer bootstrap with replacement and independent within-observer photo bootstrap; species are equal-weighted in the covariance decomposition.
- 200 realizations, seed `20260909`; first five modes matched to the zero-shift observer-equal reference by maximum total absolute cosine.
- Reference support: **331 species / 1,848 local states**. Across realizations: species 314–335 (median 326); states 1,796–1,950 (median 1,877).

## Recurrence result

| Mode | Ref. variance | median cosine | 2.5–97.5% cosine | P(cos>=0.80) | P(cos>=0.90) | recurrent dominant palette coordinates |
|---|---:|---:|---:|---:|---:|---|
| M1 | 0.320 | 0.992 | 0.962–0.999 | 1.000 | 1.000 | white 1.000, green 1.000 |
| M2 | 0.220 | 0.989 | 0.958–0.997 | 1.000 | 1.000 | white 1.000, green 0.725, brown 0.260 |
| M3 | 0.133 | 0.983 | 0.955–0.996 | 1.000 | 1.000 | brown 1.000, black 0.995 |
| M4 | 0.086 | 0.961 | 0.711–0.992 | 0.905 | 0.785 | yellow 1.000, black 0.825, bronze 0.175 |
| M5 | 0.074 | 0.950 | 0.731–0.986 | 0.905 | 0.760 | bronze 1.000, brown 0.575, yellow 0.215 |

The first three modes are especially stable: every one of the 200 realizations aligns to its reference mode with cosine >=0.90. Modes 4–5 are less individually stable and exchange ranks in 13% of realizations, although the five-dimensional subspace as a whole remains highly stable.

First-five subspace mean singular-value similarity: median **0.995**, 2.5–97.5% **0.983–0.998**. Minimum principal cosine across the five-dimensional subspace: median **0.977**, 2.5% **0.921**.

## What the recurrent modes look like

- **M1:** green +0.708, white -0.641, brown -0.212, black -0.139.
- **M2:** white +0.666, green +0.495, brown -0.428, yellow -0.262.
- **M3:** brown +0.677, black -0.565, yellow -0.402, bronze +0.189.
- **M4:** yellow -0.805, black +0.508, white -0.191, pink +0.173.
- **M5:** bronze -0.879, brown +0.323, green +0.205, pink +0.191.

The strongest reference mode is therefore a **white–green contrast axis** in the photo-derived flower-minus-background palette representation, not a prespecified whitening axis. M3 is dominated by **brown–black** structure. These names are descriptive coordinates only; ecological meaning is intentionally deferred.

## Palette-wise within-species variability

Across the 200 realizations, species-equal mean absolute local-state deviation is largest for:

- green: median 0.0995 (2.5–97.5% 0.0947–0.1041)
- white: median 0.0956 (2.5–97.5% 0.0913–0.0999)
- brown: median 0.0759 (2.5–97.5% 0.0730–0.0791)
- black: median 0.0648 (2.5–97.5% 0.0617–0.0681)
- yellow: median 0.0495 (2.5–97.5% 0.0460–0.0530)
- bronze: median 0.0475 (2.5–97.5% 0.0450–0.0511)

This recurring ordering is a descriptive property of the current palette/measurement representation. It must not be translated directly into pollinator salience without calibrated receptor-space validation.

## Reproducibility and hard stop

The full 200-realization run was executed twice with the same frozen seed and implementation; all eight output files were byte-identical. A second production implementation using an analytic WGS84/EPSG:6933 transform reproduced all seven science-bearing CSV outputs byte-for-byte.

No island/mainland, area, isolation, climate, latitude, realm, pollinator, or reproductive covariate has been used to choose or rotate the recurrent modes in this result. The next stage may overlay those variables only after this recurrence result is retained, and such associations remain exploratory until independent replication.
