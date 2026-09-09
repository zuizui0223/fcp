# RGFCA flower-side candidate genus-split robustness — protocol

Date frozen: 2026-09-09 JST

This is a second post-discovery robustness check of the same **13 fixed flower-side candidate reporting cells**. It tests whether candidate recurrence depends on closely related species from the same genera appearing in both halves.

## Fixed candidate set

Use exactly the same 13 candidate `(rx, ry, mode, full-data sign)` rows used by the species-split check. No new candidate cells are searched.

## Genus-disjoint split

Define genus as the first token of the exact retained species string. Assign entire genera using SHA-256 of that genus string:

- split A if the integer SHA-256 digest is even;
- split B if odd.

The metadata-only census before opening outcomes gives 371 genera total: 186 genera / 255 species in A and 185 genera / 245 species in B.

## Fixed state and score construction

Reuse the exact 200 recurrence realizations (`seed=20260909`), grid shifts, observer/photo bootstrap, admitted local states, frozen reference-aligned M1–M3 axes, reporting cells, and flower-only state score `S_F`.

For each realization, fixed candidate, and genus split:

1. retain states from species whose genus belongs to that split;
2. average multiple states of the same species within the reporting cell;
3. require at least 5 distinct species in the split/cell/realization;
4. compute the equal-species mean flower-only score.

Report support realizations, recurrence of the original full-data sign, median and 2.5–97.5% flower-only score, and median species count.

A candidate is labelled **cross-genus-split robust** only if both halves have at least 100 supported realizations and each has original-sign recurrence >=0.70. Apply this gate to all 13 candidates, not only the three that passed the prior species split.

## Claim boundary

This remains a post-discovery robustness diagnostic. Passing both species-disjoint and genus-disjoint splits strengthens taxonomic stability but is not an independent data replication, phylogenetic analysis, or evidence for adaptation or a biological boundary.
