# RGFCA surviving flower-side candidate observer-disjoint nuisance audit — protocol

Date frozen: 2026-09-09 JST

This is a post-discovery photographic-nuisance audit of the sole candidate surviving the species-disjoint, genus-disjoint, leave-one-genus, and season/year-matched audits: reporting cell `(-23,11)`, M2 positive. It does not search new cells or modes.

## Fixed observer split

Assign every eligible photo to observer half A or B by SHA-256 parity of its exact `observer_id` string. All photos from one observer remain in one half. The split is therefore observer-disjoint by construction and does not depend on colour outcomes.

## Metadata-only support census

Before opening observer-split M2 outcomes, the exact 200 frozen grid shifts were inspected using coordinates, species and observer identities only. Retaining the original local-state requirements within each observer half (`>=3` photos and `>=2` distinct observers) produced median candidate support of 8 species in A and 11 in B. Both halves had at least 5 candidate species in 191/200 realizations.

The minimum of 5 species per half is not newly tuned: it is the same `MINSP=5` already used by the previously frozen species-disjoint and genus-disjoint candidate audits.

## Fixed state and species rules

For each observer half and each of the exact 200 upstream grid shifts:

1. rebuild 300-km local states using only that half's photos;
2. require `>=3` photos and `>=2` distinct observers within the half;
3. admit a species only if it has at least two admitted states whose centroids are separated by `>=500 km`;
4. compute observer-bootstrap flower-only state vectors using the same one-photo-per-drawn-observer logic as the recurrence analysis;
5. centre flower vectors within species across that half's admitted states.

## Fixed M2 axis

Do not estimate or rotate a new colour axis from either half. For each realization project the observer-half flower vectors onto the already frozen, reference-aligned full-data M2 loading vector for that realization. This keeps the tested signal definition independent of the observer-half outcome.

## Fixed candidate statistic and gate

Within fixed reporting cell `(-23,11)`, average multiple states within species and then average species equally. A half-realization is supported when at least 5 species contribute.

For each half report supported realizations, median/min/max species support, median score, 2.5–97.5% realization quantiles, and fraction with positive score.

A half passes if it has at least 100 supported realizations and positive fraction `>=0.70`. The candidate is labelled **observer-disjoint robust** only if both halves pass. This is a descriptive nuisance-robustness gate, not a significance test.

## Claim boundary

Passing would show that the within-reserve candidate is not carried only by one deterministic set of observers. It would not remove observer-population differences between regions, camera hardware, exposure/white balance, illumination, habitat/background, or other photographic nuisance; nor would it establish biological adaptation or pollinator perception.
