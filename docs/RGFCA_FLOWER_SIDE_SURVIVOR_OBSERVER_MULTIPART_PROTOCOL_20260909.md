# RGFCA surviving flower-side candidate observer multi-partition stress test — protocol

Date frozen: 2026-09-09 JST

This is a post-discovery nuisance stress test of the sole surviving local flower-side candidate `(-23,11,M2+)`. The candidate and the first single observer split are already outcome-opened. This protocol therefore does **not** create confirmatory evidence; it prevents further cherry-picking of observer partitions by fixing a complete family of 20 new deterministic partitions before their colour outcomes are opened.

## Fixed observer partitions

Use salts `0..19`. For salt `j`, assign every observer wholly to one of two halves using SHA-256 parity of the literal key:

`rgfca-observer-multipart-v1|{j}|{observer_id}`

No observer may occur in both halves within a salt.

## Metadata-only support census completed before colour outcomes

Using only coordinates, species identity, observer identity, state support and the frozen 200 grid shifts, all 20 salts retained enough support for the fixed diagnostic rule below.

- per-half candidate-species counts across all salt × realization combinations: median about 9;
- every salt had at least 159/200 realizations in which **both** halves retained >=5 candidate species;
- all 20 salts therefore enter the outcome analysis; no salt is selected or dropped based on colour.

## Fixed state and score rules

For each salt, half and each of the exact 200 recurrence grid shifts:

1. rebuild 300 km local states using only observers assigned to that half;
2. require >=3 eligible photos and >=2 distinct observers within that half-state;
3. require the represented species to have >=2 admitted half-states with maximum centroid separation >=500 km;
4. project flower-only palette composition onto the already frozen realization-specific M2 loading;
5. species-centre M2 across the admitted half-states;
6. average candidate-cell `(-23,11)` state scores within species, then average species equally.

Within each half-state, observer/photo resampling is deterministic from a SHA-256 seed keyed by salt, half, realization, species and state cell. No result-dependent reruns are allowed.

A half-realization is supported when >=5 species contribute candidate scores. The threshold is deliberately lower than the original >=10-species reporting rule because observer-disjoint splitting sharply reduces support; it is a nuisance stress-test threshold, not a publication-grade ecological estimator.

## Fixed summaries

For each of the 20 salts and each half report:

- supported realizations;
- candidate-species median/min/max;
- candidate-score median and 2.5–97.5% quantiles;
- fraction of supported realizations with positive M2 candidate score.

A **salt passes** when both halves have >=100 supported realizations and both halves have positive fraction >=0.70.

The multi-partition stress test is labelled **observer-multipart robust** only if at least **16/20 salts pass**. Also report the distribution across salts of the weaker-half positive fraction `min(p_A, p_B)`.

## Claim boundary

Passing would show that the current local pattern is not dependent on one arbitrary observer partition. Failure would leave observer/camera/geographic sampling composition as an unresolved nuisance. Neither outcome is independent ecological replication, measurement validation, pollinator perception, adaptation, or proof of a biological transition.
