# RGFCA M3 genus-structure observer-disjoint stress test — protocol

Date frozen: 2026-09-09 JST

This post-discovery stress test evaluates whether the exploratory genus structure of recurrent D-M3 activity persists when observer composition is split completely. It does not search other modes, genera, geography or thresholds.

## Fixed target

The only primary target is flower-minus-background component `D`, mode `M3`, whose genus-structure activity test was previously frozen and yielded an observed genus R² above its permutation null.

## Observer partitions

Reuse exactly the 20 observer SHA256 partitions already frozen for observer-multipart stress tests: salts `0..19`, with observer assignment defined by parity of SHA256 of

`rgfca-observer-multipart-v1|{salt}|{observer_id}`.

No observer appears in both halves for a salt.

## Activity construction

For each salt, observer half and the first 20 fixed recurrence grid shifts:

- retain existing photo QC;
- define 300 km species × local states from that observer half only;
- require >=3 photos and >=2 observers per state within the half;
- retain species with >=2 admitted states and >=500 km maximum centroid separation;
- use observer-equal state means within the half;
- centre D state vectors within species;
- project onto the frozen reference M3 axis;
- define species shift activity as RMS of centred M3 state scores;
- retain a species activity as the median across shifts when available in >=10 of 20 shifts.

## Genus statistic

Within each salt × half, retain genera with >=2 retained species and compute species-equal one-way genus R² on D-M3 activity.

Require at least 20 retained species and at least 8 multi-species genera for a half to be evaluable.

Use exactly 2,000 genus-label permutations per half with deterministic seed derived from `20260913 + 100*salt + half`, preserving genus group sizes. Report observed R² and permutation p-like.

A half is labelled positive when p-like <=0.10 and observed R² exceeds the null median. A salt passes only when both halves are evaluable and positive.

The observer-disjoint genus-structure stress test passes only if at least 16 of 20 salts pass. All salts remain in the denominator.

## Sensitivity

Also report observed R² and p-like for genera with >=3 retained species where at least 6 such genera and 18 species remain. This sensitivity does not alter the primary gate.

## Claim boundary

Passing would show that genus structure in D-M3 geographic activity is not readily attributable to observer composition. Failing would mean the lineage pattern remains unstable to photography/sampling composition. Neither outcome establishes phylogenetic heritability, adaptation or mechanism.