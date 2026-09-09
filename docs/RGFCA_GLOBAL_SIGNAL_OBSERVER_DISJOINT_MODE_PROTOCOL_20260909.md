# RGFCA global recurrent signal modes — observer-disjoint validation protocol

Date frozen: 2026-09-09 JST

This is a post-discovery validation of the already discovered global flower-minus-background modes M1–M3. It does not test local geography, the Canadian Rockies candidate, island effects, latitude effects, or any ecological predictor.

## Fixed observer partitions

Reuse the 20 deterministic observer partitions defined before the observer multi-partition outcome:

`rgfca-observer-multipart-v1|{salt}|{observer_id}`, salts `0..19`, SHA-256 parity.

Each observer belongs wholly to one half within a salt.

## Fixed grid shifts

Use the first 20 grid shifts from the frozen 200-realization recurrence ledger (`seed=20260909`). These shifts are selected by realization order only, not by any mode or ecological outcome.

## Half-state construction

For each salt × half × grid shift:

- rebuild 300 km states from that half only;
- require >=3 eligible photos and >=2 distinct observers in the half-state;
- represent each state by the equal-observer mean of its 12-component flower-minus-background palette vector; no photo bootstrap is used in this validation;
- retain species with >=2 admitted states and maximum state-centroid separation >=500 km;
- construct the species-equal within-species covariance exactly as in the global recurrence analysis.

## Mode alignment

Extract the top five half-data eigenvectors and optimally permute/sign-align them to the frozen full-data reference M1–M5. M1–M3 are the validation targets.

For every salt × half × shift report M1–M3 cosine similarity and species/state support.

For each salt × half × mode report:

- median cosine across the 20 fixed shifts;
- fraction of shifts with cosine >=0.80;
- species-support median/min/max.

A half passes a mode when median cosine >=0.80 and at least 16/20 shifts have cosine >=0.80. A salt passes a mode when both halves pass. A mode is labelled **observer-disjoint robust** when at least 16/20 salts pass.

## Claim boundary

Passing would support the existence of a recurrent global colour-signal axis that is not dependent on one observer composition. It would not validate any local geographic hotspot or biological transition. Failure would mean the current image reserve is insufficient even for observer-invariant recovery of that global mode.
