# RGFCA recurrent-signal flower/background decomposition diagnostic

Date frozen: 2026-09-09 JST

This diagnostic follows the retained global signal recurrence and spatial recurrence stages. Its purpose is to determine whether a recurrent matched signal score `D = flower - background` is primarily associated with variation on the flower side, the matched-background side, or both.

## Fixed inputs

- Reuse the exact 200 upstream realizations and seed `20260909`.
- Reuse the same shifted 300 km local states, observer/photo bootstrap draws, >=3-photo / >=2-observer state rule, >=500 km within-species separation rule, species centering, and M1–M5 mode matching.
- Reuse the already retained reference mode vectors. No mode is refit using flower-only or background-only values.
- Reuse the fixed 500 km reporting cells and the previously retained spatial support rule when cell-level diagnostics are reported.

## Linear decomposition

For every selected photo retain three 12-anchor fraction vectors:

- `F`: flower palette fractions;
- `B`: matched-background palette fractions;
- `D = F - B`.

Within the same admitted species/local-state set, independently species-center all three vectors. The identity `Z_D = Z_F - Z_B` must hold numerically.

Project each centered vector on each reference-aligned mode `v_m`:

- `S_D = Z_D v_m`;
- `S_F = Z_F v_m`;
- `S_B = Z_B v_m`;

so `S_D = S_F - S_B`.

## Global diagnostic summaries

For modes M1–M3 report across species-equal state weights and 200 realizations:

1. correlation of `S_D` with `S_F`;
2. correlation of `S_D` with `-S_B`;
3. median species-equal mean absolute `S_D`, `S_F`, and `S_B`;
4. fraction of state scores with `abs(S_F) > abs(S_B)`;
5. fraction with `abs(S_B) > abs(S_F)`.

These are descriptive decomposition quantities, not causal variance partitions.

## Spatial diagnostic

For cells already admitted to the retained M1–M3 main spatial map, aggregate `S_F` and `S_B` by the same species-first equal weighting. Report their sign recurrence alongside the already retained `S_D` recurrence. Do not choose new cells from flower/background outcomes.

## Hard boundary

A background-dominated recurrent signal cannot be presented as evidence of flower adaptation or pollinator attraction. A flower-dominated signal is still only photo-derived floral colour structure; calibrated animal visual models and independent measurement validation remain necessary before claims about pollinator-perceived conspicuousness.
