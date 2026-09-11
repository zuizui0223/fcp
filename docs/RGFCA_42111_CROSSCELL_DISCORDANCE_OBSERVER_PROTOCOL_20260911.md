# RGFCA 42,111 cross-cell colour-discordance observer protocol

Date: 2026-09-11 JST
Status: frozen before any taxon×cell image pixel or taxon×cell flower-colour outcome is opened.

## Question

Among globally discovered plant species represented in at least two frozen equal-area discovery cells, how often do two independently observed, geographically distinct taxon×cell flower-photo anchors differ in the fixed four-state coarse flower-colour classification?

This endpoint is **cross-cell observed-colour discordance**, not an estimator of within-population coexistence and not automatically equal to species-wide Simpson D. The two taxon×cell anchors are deliberately drawn from distinct geographic cells, so the estimand emphasizes geographically expressed within-species colour heterogeneity.

## Frozen source frame

Use only the already frozen unique taxon×cell frame:

`results/rgfca_42111_taxon_cell_anchor_step8e2_20260911/taxon_cell_anchor_frame.csv.gz`

The frame contains 85,337 unique observation/photo identities, 42,111 species and 128 occupied discovery cells. No taxon×cell flower-colour outcome was opened when this protocol was frozen.

## Observer recovery

Recover current iNaturalist `user.id` for the exact frozen observation IDs in batches of at most 200. This is metadata transport only.

- do not replace a frozen observation when observer metadata is unavailable;
- retain unresolved observer identities explicitly;
- do not query alternative observations to improve pair eligibility;
- do not use image pixels or flower colour in observer recovery.

## Pair eligibility

A species is pair-eligible only if the frozen taxon×cell frame contains at least two rows satisfying all of:

1. distinct `cell_id`;
2. distinct `observation_id`;
3. distinct `photo_id`;
4. non-missing recovered observer IDs;
5. distinct observer IDs.

## Frozen pair selection

For every pair-eligible species, enumerate all eligible unordered row pairs and select exactly one by the minimum SHA256 of:

`20260911|crosscell-pair|inat_taxon_id|min(cell1,cell2)|max(cell1,cell2)|min(observation1,observation2)|max(observation1,observation2)|min(photo1,photo2)|max(photo1,photo2)`.

The rule is independent of flower colour, geographic distance magnitude, climate, literature status and later outcomes. No replacement is allowed after colour opening.

## Primary endpoint after taxon×cell measurement

For the selected pair of each eligible species:

`discordant = 1` if both anchors are classifiable into the frozen four states and their states differ; `0` if both are classifiable and equal; otherwise the species is measurement-unresolved for the discordance endpoint.

Report:

- frozen pair-eligible species count;
- both-classifiable pair count and fraction;
- mean cross-cell discordance among both-classifiable pairs with binomial/Wilson interval;
- the same endpoint descriptively by first-token genus size class and by frozen geographic separation quantiles, with no post-outcome threshold tuning;
- unresolved fraction against the full pair-eligible denominator.

## Interpretation ceiling

This design measures repeatable cross-cell observed-colour discordance at species breadth. It does not establish within-population polymorphism, genetic differentiation, adaptation, C*/S*, modal species colour, or unbiased species-wide Simpson D without additional sampling assumptions.
