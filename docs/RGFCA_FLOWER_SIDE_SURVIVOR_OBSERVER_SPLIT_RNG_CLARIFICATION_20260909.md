# RGFCA observer-disjoint nuisance audit — RNG clarification

Date frozen: 2026-09-09 JST

This clarification was committed before opening any observer-half M2 outcome. It changes no cell, mode, state threshold, species threshold, or robustness gate from `RGFCA_FLOWER_SIDE_SURVIVOR_OBSERVER_SPLIT_PROTOCOL_20260909.md`.

The protocol fixed the observer bootstrap logic but did not yet specify an order-independent bootstrap RNG stream. For reproducibility, each observer-half local state in each realization receives its own deterministic NumPy seed derived from SHA-256 of the UTF-8 string

`rgfca-observer-split-v1|<half>|<rep>|<species>|<cell_x>|<cell_y>`.

The first 16 hexadecimal SHA-256 digits are interpreted as an unsigned integer seed. The state then uses the already frozen one-photo-per-drawn-observer bootstrap: sample the state's distinct observers with replacement, with draw count equal to the number of distinct observers, and for each drawn observer select one of that observer's photos uniformly.

Because the seed is a deterministic function of metadata-only identifiers, results do not depend on state iteration order and no outcome-dependent seed choice is possible.

All other rules remain exactly as previously frozen: observer SHA-256 parity split, `>=3` photos and `>=2` observers per half-state, `>=500 km` within-species support, frozen full-data reference-aligned M2, fixed reporting cell `(-23,11)`, `MINSP=5`, and both halves requiring at least 100 supported realizations with positive fraction `>=0.70` to pass.