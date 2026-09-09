# RGFCA recurrent-mode activity genus-structure exploration — protocol

Date frozen: 2026-09-09 JST

This stage follows the failure of local geographic candidates to remain broadly species-composition robust. It asks a different exploratory question: whether the *amount* of within-species geographic variation along the already frozen recurrent M1–M3 signal axes is structured by genus.

No geographic hotspot, island label, climate variable, or local candidate is used to select species or modes.

## Fixed species universe and state construction

Reuse the existing photo QC and the first 20 fixed recurrence grid shifts from the `seed=20260909` recurrence ledger.

For each shift:

- define 300 km equal-area species × local states;
- require >=3 photos and >=2 observers per state;
- compute observer-equal state means;
- retain species with >=2 admitted states and >=500 km maximum centroid separation;
- centre state vectors within species.

Project the centred vectors onto the already frozen reference M1–M3 axes.

## Mode activity estimand

For species `s`, mode `m`, and shift `r`, define activity as the root mean squared centred state score:

`A_smr = sqrt(mean(score_state^2))`.

This is a non-directional measure of how strongly a species varies geographically along a recurrent mode. It does not encode which side of a geographic transition is positive.

For every species × mode, use the median activity across the 20 fixed shifts as the retained species activity. Require activity to be available in >=10 of 20 shifts.

Compute the analysis separately for:

1. flower-minus-background signal `D` (primary, because M1–M3 were discovered from D);
2. flower-only signal `F` (companion);
3. matched background signal `B` (diagnostic).

The exact identity `D = F - B` is retained at state-score level but RMS activities are not additive.

## Genus support gate

Derive genus as the first token of the accepted species name. The metadata-only census before activity outcomes found 331 spatially eligible species in 261 genera, including 49 genera with >=2 species (119 species total) and 13 genera with >=3 species (47 species).

Primary genus-structure analyses use all genera with >=2 retained species for the relevant component/mode. A sensitivity uses genera with >=3 retained species.

## Genus structure statistic

For each component × mode, fit the one-way genus decomposition on species activities and report

`R2_genus = SS_between_genus / SS_total`,

with species as equal-weight observations. This is descriptive variance explained by genus in the retained species set.

## Permutation null

Use random seed `20260912` and exactly 10,000 permutations. Permute genus labels among retained species while preserving the observed multiset of genus labels and therefore all genus group sizes. Recompute `R2_genus` each time.

Report the observed R², null median and 2.5–97.5% range, and

`p_like = (1 + count(R2_null >= R2_observed)) / 10001`.

This is exploratory calibration, not a confirmatory family-wise-error-controlled test.

## Robustness summaries

Also report, for each genus with >=2 species, the median species activity and within-genus range. No genus is selected for inference after outcome opening.

The >=3-species sensitivity repeats the same R²/permutation calculation without changing thresholds.

## Claim boundary

A positive result would indicate that recurrent geographic colour-signal *activity* is partly lineage-structured at genus level in the current reserve. It would not establish phylogenetic heritability, adaptation, common mechanism, or geographic convergence; formal phylogenetic inference would require an external tree and an independent analysis design.