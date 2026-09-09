# RGFCA hypothesis-free flower–background signal exploration — endpoint

Date: 2026-09-09 JST

This checkpoint closes the current-reserve hypothesis-free discovery branch. It is an endpoint, not a new candidate-selection stage.

## What survived discovery

The reserve supports recurrent, species-conditioned photo-derived flower–background signal modes without assuming a common biological boundary, island syndrome, or whitening direction in advance.

- The eligible reserve used by the recurrence analysis contains 24,885 admitted/classifiable photos from 500 species.
- Repeated 300-km local-state construction with observer/photo resampling retained roughly 320–330 geographically informative species per realization.
- The first three reference-aligned modes recur very strongly across the 200 frozen realizations; each M1–M3 alignment remained >=0.90 in all 200 realizations.
- The dominant mode compositions are approximately white/green (M1), white/green versus brown/yellow (M2), and brown versus black/yellow (M3).

These results support recurrent structure in the photo-derived colour-signal representation. They do not by themselves identify a biological transition mechanism.

## What did not survive global ecological overlay

Using all admitted M1–M3 states rather than selected cells:

- flower-only mode scores showed no stable global within-species association with absolute latitude;
- flower-only mode scores showed no stable global within-species association with island share;
- direct flower–background separation magnitude (TV and Hellinger) showed no clear monotonic global association with latitude or island share.

Matched flower-minus-background latitude patterns were more strongly attributable to matched-background structure than to flower-only shifts. In particular, background green increased with absolute latitude in the exploratory palette diagnostic, while flower-side palette slopes were weak.

Therefore the reserve does not support a simple global rule such as `island -> duller/whiter flowers` or `higher latitude -> larger flower-background separation`.

## Local candidate sieve

The fixed spatial recurrence analysis initially produced 13 flower-side candidate cells. They were then subjected to increasingly strict, outcome-fixed nuisance and taxonomic audits.

1. species-disjoint split: 3/13 retained the same direction in both halves;
2. genus-disjoint split: 4/13 retained the same direction in both halves;
3. intersection of the species- and genus-disjoint survivors: 1/13;
4. the sole intersection survivor was reporting cell `(-23,11)`, M2 positive, approximately 116.6 W / 51.7 N (Canadian Rockies region);
5. leave-one-genus analysis showed that this candidate was not carried by one genus;
6. season/year matching passed: primary positive fraction 0.93 across 200/200 supported realizations, with median 18 matched species;
7. observer-disjoint audit failed prospectively fixed replication across photographer sets.

## Decisive observer-disjoint result

The observer split, state rules, support floor and seeds were frozen before observer-half flower outcomes at commits `269845cb338b004fdd16a4687e5efa4fa6f506b7` and `897ac0c2a67f64f545c099d3c687345a3df285dd`.

| half | supported realizations | species median (range) | median flower M2 | 2.5–97.5% | positive fraction | gate |
|---|---:|---:|---:|---:|---:|---|
| A | 191 | 8 (5–14) | +0.009494 | -0.042578 to +0.056367 | 0.6178 | FAIL |
| B | 200 | 11 (5–17) | +0.014373 | -0.015158 to +0.048805 | 0.7800 | PASS |

The candidate-level rule required both halves to pass. Therefore `observer_disjoint_robust = false`.

This does not prove that observer or camera differences generated the candidate. Half A has lower support and a broad realization distribution. But the current reserve does not meet the prospectively fixed requirement that the candidate recur across completely disjoint observer sets.

## Current-reserve endpoint

**No flower-side local candidate in the current reserve survives the complete discovery -> taxonomic -> season/year -> observer-disjoint chain.**

Accordingly:

- do not promote the Canadian Rockies M2+ cell as a nuisance-robust biological transition;
- do not reopen the 13-cell set and choose another candidate after seeing these audits;
- do not relax the observer split, state support thresholds, `MINSP`, or positive-fraction gate to rescue the surviving cell;
- do not reinterpret the stable M1–M3 axes themselves as pollinator perception, adaptation, or causal flower evolution.

The defensible retained result is narrower: recurrent photo-derived flower/background signal structure exists, its global ecological overlays are weak, and the strongest flower-side local candidate does not survive a fully observer-disjoint nuisance gate.

## Next admissible stage

Further biological interpretation requires a genuinely independent replication tranche with new observer/photo support and an independently validated flower-region measurement target. The existing reserve is now outcome-opened for all candidate-selection and nuisance-audit decisions described above and cannot serve as its own confirmatory replication.

The previously failed Monarda localization-agreement gate remains relevant: a new ecological replication should not be promoted as biological evidence until the measurement-validity lane establishes an adequate independent target/reference. Until then, new ecological acquisition may be designed prospectively but should not be used to rescue the current-reserve candidate.
