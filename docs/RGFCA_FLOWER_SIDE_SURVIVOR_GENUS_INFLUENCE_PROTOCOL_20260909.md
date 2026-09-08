# RGFCA surviving flower-side candidate genus-influence audit — protocol

Date frozen: 2026-09-09 JST

After applying both fixed species-disjoint and genus-disjoint robustness gates to all 13 previously discovered flower-side candidate cells, only one candidate passed both: `(rx=-23, ry=11)`, mode M2, positive full-data sign. This audit tests whether that surviving pattern is dominated by any single genus.

## Fixed candidate

Use only the already fixed surviving candidate `(rx=-23, ry=11, M2, positive)`; do not search new cells or modes.

## Fixed realizations and scores

Reuse the exact 200 recurrence realizations (`seed=20260909`), grid shifts, observer/photo bootstrap, admitted local states, frozen M1–M3 axes, 500-km reporting cells, and flower-only score `S_F`.

Within each realization and the fixed reporting cell:

1. average multiple states of the same species;
2. require at least 10 species for the full-cell estimate;
3. compute the equal-species mean M2 flower-only score;
4. define genus from the first token of each species string;
5. recompute the equal-species mean after removing each genus in turn, whenever at least 5 species remain.

## Fixed robustness summaries

Across supported realizations report:

- recurrence of the positive full-data sign for the full-cell estimate;
- fraction of realizations in which **every** eligible leave-one-genus-out estimate remains positive;
- median of the minimum leave-one-genus-out score per realization;
- 2.5–97.5% quantiles of that minimum;
- maximum absolute single-genus influence `|LOO - full|` per realization and its median/quantiles;
- for each genus, number of eligible removals, fraction of removals that flip the sign, and median absolute influence.

The candidate is labelled **single-genus robust** if at least 70% of supported realizations have every eligible genus deletion preserve the positive sign. This is a descriptive robustness gate, not a significance test.

## Claim boundary

This is a post-discovery influence audit. Passing it would show that the candidate is not easily removed by deleting one genus, but it would still not constitute independent ecological replication, a phylogenetic test, adaptation evidence, or proof of a biological boundary.
