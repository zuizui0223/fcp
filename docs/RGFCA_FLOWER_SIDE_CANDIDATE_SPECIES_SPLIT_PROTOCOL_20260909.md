# RGFCA flower-side candidate species-split robustness — protocol

Date frozen: 2026-09-09 JST

This is a post-discovery robustness check of the **13 already fixed flower-side candidate reporting cells** from the completed spatial/component analysis. It does not search for new candidate cells and cannot convert the existing exploratory discovery into a prospective confirmation.

## Fixed candidate set

Use exactly the 13 rows retained in `exploratory_flower_side_candidate_cells.csv` from the completed flower/background decomposition. Keep each candidate's fixed reporting cell `(rx, ry)`, mode M1–M3, and full-data sign.

## Species-disjoint split

Assign every reserve species independently of outcomes using SHA-256 of its exact species string:

- split A if the integer SHA-256 digest is even;
- split B if it is odd.

The metadata-only census before opening split outcomes gives 253 species in A and 247 in B.

## Fixed state and score construction

Reuse the exact 200 recurrence realizations (`seed=20260909`), grid shifts, observer/photo bootstrap, admitted local states, and frozen reference-aligned M1–M3 axes. Use the flower-only state score `S_F` from the completed decomposition.

For each realization, candidate cell, mode, and split:

1. retain only admitted states belonging to species in that split;
2. average multiple states of the same species within the reporting cell;
3. require at least 5 distinct species in that split/cell/realization;
4. compute the equal-species mean flower-only score.

## Robustness summaries

For each fixed candidate and split report:

- number of supported realizations;
- fraction of supported realizations whose sign matches the candidate's full-data sign;
- median split-specific flower-only score;
- 2.5–97.5% quantiles.

A candidate is labelled **cross-split robust** only if both A and B have at least 100 supported realizations and each has full-sign recurrence >=0.70. This threshold is a descriptive robustness gate, not a significance test.

Also report whether the two split medians have the same sign as the full-data candidate.

## Claim boundary

Because the candidate cells were discovered using all species before this split check, passing this gate demonstrates taxonomic robustness only. It is not an independent confirmatory replication and does not establish adaptation, a biological boundary, pollinator causation, or generality outside the current reserve.
