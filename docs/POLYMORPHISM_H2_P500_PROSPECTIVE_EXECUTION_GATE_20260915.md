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
4. dry-run the complete 49,900-row location-blind firewall and verify all 256 partitions;
5. bind the receipt to the exact branch and 40-character `GITHUB_SHA`;
6. record all preopening outcome flags as false.

The resulting artifact is forward-only evidence. A later image-opening workflow must consume and validate that exact receipt rather than synthesize its fields.

## Frozen execution ladder

The global execution may advance only one stage at a time:

`PREOPENING -> FIREWALL_FROZEN -> PARTITION_MEASUREMENT_COMPLETE -> REASSEMBLY_COMPLETE -> SUPPORT_GATE_COMPLETE -> H2_COMPLETE`

This ordering matches the safe implementation: each partition acquires its images, immediately measures them location-blind, seals terminal rows, and deletes pixels. The design therefore never requires all 49,900 images to be retained simultaneously.

### FIREWALL_FROZEN

- exactly 49,900 unique blinded measurement IDs;
- two blind batches;
- 32 semantic shards per batch;
- four compute partitions per shard;
- all 256 terminal partitions populated;
- worker surface contains no species, coordinates, source URL, observer identity or biological outcome;
- candidate pixels remain unopened.

### PARTITION_MEASUREMENT_COMPLETE

Across the complete 256-partition set:

- exactly 49,900 frozen rows have terminal acquisition states;
- exactly 49,900 rows have terminal measurement states;
- all 256 terminal partitions are present;
- no image pixels or masks are persisted after partition sealing;
- no replacement rows or species.

### REASSEMBLY_COMPLETE

- exactly 49,900 unique measurement IDs;
- zero duplicate measurement IDs;
- the sealed metadata-colour join opens only after complete partition coverage.

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

## Current qualification result

The branch qualification has already reconstructed the real frozen inputs and dry-run the complete P500 firewall without opening pixels. The dry run retained 499 species / 49,900 rows and populated all 256 partitions under the exact inherited blind-measurement interface.

The H2 executor was additionally qualified by rebuilding the legacy discovery/reserve Delta vectors, matching the committed frozen `result.json`, and requiring exact equality of all 999 structured-null draws between the legacy and prospective executors in all four primary/strict × discovery/reserve cases.

A legacy 34-species figure workflow subsequently moved the PR head without changing any P500 scientific input or implementation. Because PREOPENING evidence is bound to an exact commit SHA, that head move is treated conservatively as invalidating the earlier head-specific qualification. This document update intentionally triggers a fresh P500 qualification on the new post-figure head; no first pixel may open until that fresh qualification succeeds.

This remains a technical qualification result only. No image was fetched or decoded by the qualification workflow, no ROI/segmentation was run, no colour outcome was opened, and no H2 statistic was calculated. The prospective biological result remains unopened until a PREOPENING receipt is explicitly issued and consumed by a separately qualified execution workflow.
