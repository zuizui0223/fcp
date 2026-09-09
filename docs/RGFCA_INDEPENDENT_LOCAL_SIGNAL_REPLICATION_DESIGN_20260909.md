# RGFCA independent local-signal replication — prospective design

Date frozen: 2026-09-09 JST

This document begins the only admissible follow-up to the closed current-reserve exploration endpoint. It does **not** authorize ecological outcome opening. The current reserve is outcome-opened and cannot serve as its own replication.

## Prerequisite measurement gate

No ecological replication result may be opened unless a separately frozen independent measurement-validity lane first establishes an adequate target/reference for the quantity being measured. The failed Monarda localization-agreement gate and the incomplete FlowerMask paired frame remain binding limitations.

Before ecological pixels are measured, the measurement lane must specify and pass, on independent validation data:

- the visible flower-region target actually annotated by the reference;
- source/version/hash and rights provenance;
- frozen model and image transforms;
- denominator, unknown/ambiguous-tissue semantics and technical-failure handling;
- development versus untouched validation membership;
- a prospectively justified pass/fail agreement gate.

If that gate does not pass, the ecological replication remains blocked. Do not lower the gate because the current-reserve M2+ candidate exists.

## Fixed ecological target

The independent ecological target is fixed to the sole current-reserve candidate that survived species-, genus-, leave-one-genus and season/year audits before failing observer-disjoint replication:

- reporting cell: `(-23,11)` on the frozen 500-km reporting lattice;
- approximate region: Canadian Rockies, around 116.6 W / 51.7 N;
- signal mode: **M2**;
- predicted direction: **positive flower-only M2** in the target cell relative to same-species geographically separated comparison support.

No other cell, mode, sign, island hypothesis, latitude trend or candidate transition is eligible for confirmatory promotion in this tranche.

## Independent data identity

Replication photos must be disjoint from the current reserve at all of the following identifiers before any colour outcome is opened:

- `photo_id`;
- `observation_id`;
- exact `observer_id`.

All observers in the replication tranche must therefore be new relative to the discovery reserve. Species are allowed to overlap because the estimand is within-species geographic change; species overlap does not make the photographs or observer population non-independent.

The exact acquisition manifest must be frozen and hashed before measurement. No candidate-region colour, segmentation score or M2 value may be used to select photos, species or observers.

## Metadata-only qualification gate

Qualification uses coordinates, dates, species and observer IDs only. It occurs before opening replication image measurements.

A species qualifies only if it has both:

1. target-cell support in fixed reporting cell `(-23,11)`; and
2. at least one same-species comparison state outside the target cell whose centroid is >=500 km from the target support and whose circular mean day-of-year differs by <=30 days.

For both target and comparison support require, before pixel opening:

- >=6 photographs;
- >=3 distinct observers.

If several comparison states qualify, all are retained and receive equal state weight. No nearest/best-colour comparison is selected.

The ecological tranche is authorized only if **>=20 species** satisfy this full metadata-only rule. This floor is chosen to meet or exceed the roughly 19-species median support of the discovery candidate without using replication outcomes. If fewer than 20 species qualify, record `replication_support_qualified=false`; do not relax the photo, observer, distance, season or species floor after seeing the census.

## Frozen signal definition

Do not refit PCA or discover a new mode on replication data. Project flower-only measurements onto the exact frozen discovery **reference M2 loading vector** retained by the hypothesis-free recurrence analysis.

For each state:

- average photographs within observer first;
- then average observers equally;
- compute the flower-only M2 state score.

For each species:

- average all target states equally;
- average all eligible comparison states equally;
- define `Delta_s = target_M2 - comparison_M2`.

The primary estimator is the equal-species mean of `Delta_s`.

Matched flower-minus-background M2 and matched-background M2 are required companion diagnostics, using the same species/state weights, but they cannot replace the flower-only primary endpoint.

## Confirmatory resampling and null

Use 1,000 prespecified bootstrap realizations after the qualified manifest and measurement model are frozen.

Within each state, resample observers with replacement and then one photograph within each drawn observer. Recompute `Delta_s` and the equal-species mean each realization.

Report:

- median replicated flower-only M2 contrast;
- 2.5–97.5% bootstrap quantiles;
- fraction of bootstrap realizations >0.

For a fixed-direction calibration, use a species-level sign-flip null on the observed species contrasts with 10,000 deterministic random sign assignments, seed `20260914`. Report the one-sided probability of a null equal-species mean >= the observed equal-species mean.

## Observer-disjoint replication gate

The discovery candidate failed the observer-disjoint gate, so observer independence is mandatory in the replication tranche.

Before outcomes, partition the **new** replication observers by SHA-256 parity of exact observer ID. Recompute the fixed M2 contrast independently in each half without refitting the M2 axis.

A half is evaluable only if >=10 species retain both target and comparison support after the observer split. The observer-disjoint replication condition passes only if:

- both halves are evaluable; and
- both half-level equal-species contrasts are positive.

No half-specific threshold tuning or reassignment is allowed.

## Replication success rule

The independent local-signal replication is called supported only if all of the following hold:

1. prerequisite independent measurement-validity gate passed;
2. metadata-only ecological support qualification passed with >=20 species;
3. primary flower-only M2 bootstrap median >0;
4. primary 2.5% bootstrap quantile >0;
5. one-sided 10,000-draw species sign-flip calibration <=0.05;
6. observer-disjoint replication condition passed in both new-observer halves.

Failure of any one item is retained as replication non-support. Do not rescue a failure by reopening other cells, modes, colour axes, date windows, distance thresholds or observer partitions.

## Claim boundary

Even a successful replication would establish only a reproducible photo-derived local flower-signal shift under the validated measurement target and the fixed M2 representation. It would not by itself establish adaptation, a pollinator-perceived colour distance, causal environmental selection, an island syndrome or a universal flower-colour boundary. Those mechanisms require separately designed evidence.
