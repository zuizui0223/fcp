# RGFCA surviving flower-side candidate season/year matched nuisance audit — protocol

Date frozen: 2026-09-09 JST

This is a post-discovery nuisance audit of the sole flower-side candidate that survived both species-disjoint and genus-disjoint partitions and the leave-one-genus audit: reporting cell `(-23,11)`, mode M2, positive flower-only direction. It does not search new cells, modes, ecological predictors, or date thresholds from signal outcomes.

## Metadata-only support census

Before opening any season-matched M2 outcome, the exact 200 recurrence grid shifts were inspected using only coordinates, species identity, `observed_on`, state support and the fixed candidate cell. The candidate contained 13–27 species per realization (median 19). Requiring a same-species comparison state at least 500 km away, within 30 circular day-of-year days and within 5 mean observation years retained 12–27 matched species per realization (median 18). This condition is therefore fixed as the primary nuisance match.

## Fixed realizations and score

Reuse the exact 200 recurrence realizations (`seed=20260909`), grid shifts, observer/photo bootstrap, state admission, frozen reference-aligned M1–M3 axes and flower-only state score `S_F`.

Use only fixed candidate `(-23,11,M2+)`.

## State date metadata

For each admitted local state, calculate from **all eligible photos assigned to that state**, not the bootstrap-selected photos:

- circular mean day of year, treating one year as 365.25 days;
- mean observation year.

Date metadata therefore cannot respond to the bootstrapped colour outcome.

## Primary matched comparison

For every candidate state of a represented species, eligible comparison states must:

1. belong to the same species;
2. fall outside the fixed candidate reporting cell;
3. have centroid distance at least 500 km from that candidate state;
4. have circular mean day-of-year distance <=30 days;
5. have absolute mean observation-year difference <=5 years.

If multiple comparison states satisfy these rules, average all of their flower-only M2 scores with equal state weight. For each candidate state compute candidate M2 minus that comparison mean. Average multiple candidate-state contrasts within species, then average species equally.

A realization is supported only if at least 10 species have a valid matched contrast.

## Fixed summaries and robustness gate

Across supported realizations report:

- number of supported realizations;
- matched-species median/min/max;
- median matched M2 contrast and 2.5–97.5% quantiles;
- fraction of supported realizations with positive matched contrast.

The candidate is labelled **season/year-match robust** if at least 100 realizations are supported and the positive fraction is >=0.70. This is a descriptive nuisance-robustness gate, not a significance test.

Two prespecified metadata-only sensitivities are reported without changing the primary gate:

- season-only: <=30 circular days, no year restriction;
- stricter-year: <=30 circular days and <=2 mean observation years.

Each retains the same >=500 km spatial separation and >=10 matched-species realization support.

## Existing matched-background context

The already completed component analysis is not reselected here. For this fixed cell, the outcome-opened full-data recurrences were F positive 0.95, D positive 0.82 and B positive 0.70. The present audit tests date/season dependence only; it does not reinterpret those component results.

## Claim boundary

Passing this audit would show that the surviving within-reserve M2 pattern is not readily removed by coarse matching on flowering season and observation year. It would not eliminate illumination, exposure, camera-processing or other photographic nuisance, and it would not constitute independent ecological replication, pollinator perception, adaptation, or proof of a biological transition.
