# P500 prospective execution epoch

Date: 2026-09-15 JST
Epoch: `P500_H2_PROSPECTIVE_20260915_E1`
Base: `analysis/polymorphism-42111-h1-h2-gates-20260912` at `26f7a792f93fb53bc6e8c9a64efb45ccb548b9ea`

## Purpose

Create a new, auditable chronology boundary for the untouched P500 white-axis confirmation without inventing missing historical evidence.

The earlier measurement-control audit correctly concluded that the saved candidate-metadata result is not a complete pre-opening execution receipt. This epoch does not reinterpret that record. Instead, it starts a new prospective state machine whose first durable state is committed before any execution belonging to this epoch may begin.

This document and its JSON receipt do **not** open images, image bytes, flower colour, palette fractions, D, the metadata-colour join, H2 vectors, W, or an H2 p-value.

## Bound scientific design

The epoch inherits the exact fixed P500 science:

- 499 species and 49,900 frozen candidate rows;
- 100 rows per species;
- no species or row replacement;
- location-blind measurement machine from source commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`;
- complete 256-partition terminal census before metadata-colour join;
- >=40 classifiable photos for a measurement-evaluable species;
- >=250 measurement-evaluable species before H2;
- fixed primary 10% and strict 20% H1/H2 admission thresholds;
- fixed `q_white` and W statistic;
- 999 structured-null worlds with seeds 20260915 and 20261015;
- no axis refit, target relaxation, replacement, or outcome-driven rerun.

The support-floor conflict between the 2026-09-13 selection protocol and the 2026-09-15 measurement protocol is resolved in `POLYMORPHISM_H2_P500_GATE_RECONCILIATION_20260915.md`: 20 vectors remain a computational floor, but the earlier **N>=100** floor remains binding for a confirmatory supported/not-supported verdict.

## State machine

Only four ordered states exist:

1. `PREOPEN_FROZEN`
2. `DISPATCH_RECORDED`
3. `MEASUREMENT_COMPLETE`
4. `INFERENCE_OPENED`

No state may be skipped.

### PREOPEN_FROZEN -> DISPATCH_RECORDED

A first-dispatch event must be durably recorded while all of the following are still false:

- image requests started;
- image bytes opened;
- measurement started;
- H2 opened.

The event must bind the canonical SHA-256 of the exact prior epoch receipt, 499 species and 49,900 rows. The returned advanced receipt must be saved durably **before** the acquisition/measurement action begins.

### DISPATCH_RECORDED -> MEASUREMENT_COMPLETE

The completion event is valid only for exactly:

- 49,900 terminal rows;
- 49,900 unique measurement IDs;
- 256 terminal partitions;
- zero duplicate IDs;
- zero missing IDs;
- no early stopping;
- no replacement.

Failures remain terminal reason-coded rows and do not authorize substitution.

### MEASUREMENT_COMPLETE -> INFERENCE_OPENED

The metadata-colour join/H2 surface may open only after the complete terminal census is verified and a separate inference-opening record is durably saved. The record must assert that neither the join nor H2 had already been opened before that record.

## Provenance ceiling

This epoch establishes chronology only **from the committed epoch receipt forward**. It explicitly does not claim that all historical non-access before this freeze has been independently reconstructed or proven.

A successful offline test suite certifies only the state-transition safeguards and immutable design constants. It is not a biological result and does not itself authorize pixel opening.

## Current state

`PREOPEN_FROZEN`.

No dispatch receipt has been created. No measurement workflow is authorized to acquire P500 images from this branch yet. The next admissible action is to create and durably save one exact `first_dispatch_record`; only a subsequent execution step may begin image acquisition.
