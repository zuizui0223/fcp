# RGFCA South Florida M2 flower-side species-multipart stress test — protocol

Date frozen: 2026-09-09 JST

This is a post-discovery stress test of the sole observer-consensus, flower-supported local signal candidate: South Florida reporting cell `(-16,6)`, mode `M2`, positive flower-only direction. It does not search new cells, modes, salts, thresholds, or ecological predictors from signal outcomes.

## Fixed target and upstream objects

Reuse the frozen observer-consensus geography and component-attribution results. The only target is reporting cell `(-16,6)`, `M2+`, flower-only score `S_F`.

Use the first 20 fixed recurrence grid shifts from the frozen `seed=20260909` recurrence ledger and the frozen reference M2 loading vector. State geometry remains 300 km equal-area cells; reporting geometry remains 500 km cells.

## Twenty species-disjoint partitions

Use exactly 20 salts, integers `0..19`. Assign every species name wholly to half A or B using the parity of SHA256 of the UTF-8 string

`rgfca-species-multipart-v1|{salt}|{species}`.

No species may appear in both halves for a salt. All observers belonging to an admitted species remain within that species half.

## State construction within each species half

For each salt, half, and grid shift:

- keep the existing photo QC: admitted ROI, globally classifiable, positive 12-colour flower/background totals, positional accuracy <=5 km;
- require >=3 photos and >=2 observers per species × 300 km state;
- compute flower-only palette fractions;
- within each state, first average photographs within observer, then average observers equally;
- for each species with >=2 admitted states and >=500 km maximum state-centroid separation, centre its state flower vectors by that species' mean and project onto the frozen reference M2 axis;
- within reporting cell `(-16,6)`, average multiple candidate states of a species, then average species equally.

A half × shift is supported only when the candidate cell contains >=5 informative species.

## Half and salt gates

For each salt and half, across the 20 fixed shifts report the number of supported shifts, median candidate score, positive fraction, and median candidate-species support.

A half passes when:

- at least 10 of 20 shifts are supported; and
- the positive-score fraction among supported shifts is >=0.70.

A salt passes only when both A and B pass. The candidate is labelled `species_multipart_robust=true` only if at least 16 of the 20 salts pass.

Also report the weaker-half positive fraction for every salt and its median/range across salts. All salts remain in the denominator; none may be dropped after outcome opening.

## Claim boundary

Passing would show that the South Florida M2 flower-side pattern is not dependent on one particular species composition split in the current reserve. Failing would mean taxonomic composition remains a material source of instability even though the signal passed observer-consensus and flower/background attribution. Neither outcome establishes adaptation, pollinator perception, causality, or independent biological replication.