# RGFCA observer-consensus geography — protocol

Date frozen: 2026-09-09 JST

This stage rebuilds geography after the sole pooled-map flower-side candidate failed observer multi-partition robustness, while the global M1–M3 mode basis passed observer-disjoint validation. No pooled hotspot or previous candidate is privileged in this analysis.

## Fixed inputs

- the frozen full-data reference axes M1–M3;
- the 20 deterministic observer partitions `rgfca-observer-multipart-v1|{salt}|{observer_id}`, salts `0..19`;
- the first 20 frozen 300 km grid shifts from the recurrence ledger;
- the fixed 500 km reporting lattice.

## Half-state construction

For every salt × observer half × grid shift:

1. rebuild 300 km states using only that half's observers;
2. require >=3 eligible photos and >=2 distinct observers per half-state;
3. represent the state by the equal-observer mean flower-minus-background 12-colour vector;
4. retain species with >=2 admitted states and >=500 km maximum state-centroid separation;
5. species-centre the state vectors and project them onto the fixed full-data M1–M3 axes;
6. aggregate each reporting cell by species-equal mean mode score.

A half-cell-shift is supported when >=5 species contribute.

## Metadata-only support census

Before any observer-consensus colour outcome was opened, support was checked from coordinates/species/observer metadata only. With the >=5-species threshold, 30 reporting cells are evaluable in both halves for at least 16/20 salts when a half is required to support the cell in >=10/20 shifts. This threshold is therefore fixed.

## Direction recurrence within a half

For each salt × half × reporting cell × mode:

- require the cell to be supported in >=10/20 shifts;
- call the half direction **positive** when >=70% of supported shift scores are >0;
- call it **negative** when <=30% are >0;
- otherwise call it unstable.

## Agreement across observer halves and salts

A salt agrees for a cell × mode only when both observer halves are directional and have the same sign.

A cell × mode is labelled **observer-consensus geographic** only when:

1. at least 16/20 salts are evaluable in both halves for that cell; and
2. at least 80% of evaluable salts agree on the same direction.

All supported cells and M1–M3 are retained in reporting; no cell is selected based on pooled-map strength.

## Claim boundary

This is a stability map, not a significance map. Surviving cells show geographic mode direction that recurs across disjoint observer compositions and grid shifts in the current reserve. They do not establish adaptation, pollinator perception, causal environmental effects, or universal transition boundaries.
