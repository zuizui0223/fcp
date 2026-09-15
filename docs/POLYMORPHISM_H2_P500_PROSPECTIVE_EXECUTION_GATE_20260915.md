# P500 prospective execution gate

Date: 2026-09-15 JST

## Purpose

This gate closes the implementation gap between the frozen P500 metadata/H2 design and the first prospective pixel opening. It does **not** reinterpret the earlier measurement-control receipt as a verified historical execution record. Instead it establishes a new forward-only chronology from a clean branch state and requires exact stage-by-stage evidence before H2 can open.

## Frozen scientific inputs

The gate accepts only the already frozen P500 state:

- P500 species-set SHA256: `f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4`;
- candidate metadata SHA256: `a2339a3eba7bec8c29e726edc8c71a64a1a98b5cad4367764f7466c436e9f595`;
- exact-100 species SHA256: `568a64e692dc21f4a5f76c5eadb615a63b07d71f01bc9cd91e8496ccd0fa3f9e`;
- 499 species and 49,900 frozen image rows;
- inherited measurement source commit: `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`;
- primary H2 target: fixed `q_white` / `W` structured-null test;
- no axis refit, row/species replacement, target relaxation, threshold change, or post-outcome rescue.

The gate validates the real candidate-metadata result and the frozen 2026-09-15 prospective authorization before permitting any forward start receipt.

## Why this is separate from the 2026-09-13 white-control freeze

The earlier white-control implementation correctly found that the candidate metadata file is not itself a complete preopening execution record. Missing historical flags are not backfilled. A static declaration cannot prove that an execution stage actually occurred.

This new gate therefore starts a new auditable chronology **before any P500 pixel opening**. It does not claim to certify unrecorded external access before this start point. Its only role is to ensure that the future execution can no longer advance without machine-checkable receipts.

## Forward start receipt

Ordinary push/PR CI is qualification-only and always records `opening_authorized=false`.

A PREOPENING receipt can be created only by an explicit workflow dispatch with `issue_start_receipt=true`. In that run the workflow must first:

1. validate the real frozen candidate metadata and prospective authorization;
2. run the execution-gate test suite;
3. verify that the designated prospective output namespaces do not yet exist;
4. bind the receipt to the exact branch and 40-character `GITHUB_SHA`;
5. record all preopening outcome flags as false.

The resulting artifact is forward-only evidence. A later image-opening workflow must consume and validate that exact receipt rather than synthesize its fields.

## Frozen execution ladder

The execution may advance only one stage at a time:

`PREOPENING -> ACQUISITION_COMPLETE -> MEASUREMENT_COMPLETE -> REASSEMBLY_COMPLETE -> SUPPORT_GATE_COMPLETE -> H2_COMPLETE`

Each transition fails closed.

### ACQUISITION_COMPLETE

- exactly 49,900 frozen rows must have terminal acquisition states;
- no replacement rows or species.

### MEASUREMENT_COMPLETE

- exactly 49,900 terminal measurement rows;
- exactly 256 terminal partitions;
- no replacement rows or species.

### REASSEMBLY_COMPLETE

- exactly 49,900 unique measurement IDs;
- zero duplicate measurement IDs.

### SUPPORT_GATE_COMPLETE

- measurement-evaluable species are counted using the frozen `n_classifiable >= 40` rule;
- `>=250` gives `PASS`;
- `<250` gives `NOT_EVALUABLE` and H2 must remain closed.

### H2_COMPLETE

H2 can open only after support `PASS`.

- fewer than 20 primary H2 vectors must yield `NOT_EVALUABLE_VECTOR_SUPPORT`;
- otherwise the fixed 999-world structured-null p-value is evaluated;
- `p < 0.05` -> `CONFIRMED`;
- `p >= 0.05` -> `NOT_CONFIRMED`.

The strict 20% sensitivity cannot rescue the primary 10% decision.

## Current state after this change

This branch adds and tests the forward execution contract. It does **not** fetch images, decode pixels, run ROI/segmentation, open colour, compute D, or calculate H2. Until an explicit PREOPENING receipt is issued and consumed by a separately qualified execution workflow, the prospective P500 biological result remains unopened.
