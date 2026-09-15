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

## Forward start receipt

Ordinary push/PR CI is qualification-only and always records `opening_authorized=false`.

A PREOPENING receipt can be created only by an explicit workflow dispatch with `issue_start_receipt=true`. In that run the workflow must first:

1. validate the real frozen candidate metadata and prospective authorization;
2. run the execution-gate test suite;
3. verify that the designated prospective output namespaces do not yet exist;
4. validate the `SUPPORT_GATE_COMPLETE -> H2_COMPLETE` receipt contract;
5. replay the frozen H2 implementation against all four legacy cases;
6. dry-run the complete 49,900-row location-blind firewall and verify all 256 partitions;
7. bind the receipt to the exact branch and 40-character `GITHUB_SHA`;
8. record all preopening outcome flags as false.

The resulting artifact is forward-only evidence. The image-opening workflow consumes and validates that exact receipt rather than synthesizing its fields.

## Frozen execution ladder

`PREOPENING -> FIREWALL_FROZEN -> PARTITION_MEASUREMENT_COMPLETE -> REASSEMBLY_COMPLETE -> SUPPORT_GATE_COMPLETE -> H2_COMPLETE`

The implementation realizes the complete ladder. Each of the 256 blind partitions acquires its frozen image rows, measures colour without species/coordinate/source-URL fields on the worker surface, seals terminal rows, and destroys image pixels. Reassembly opens the metadata-colour join only after all 49,900 frozen rows and all 256 partition receipts are present. H2 opens only after the frozen measurement-support gate passes.

### SUPPORT_GATE_COMPLETE

- measurement-evaluable species use the frozen `n_classifiable >= 40` rule;
- `>=250` gives `PASS`;
- `<250` gives `NOT_EVALUABLE`, and H2 remains unopened.

### H2_COMPLETE

- fewer than 20 primary H2 vectors gives `NOT_EVALUABLE_VECTOR_SUPPORT`;
- otherwise the fixed 999-world structured-null test is evaluated;
- primary `p < 0.05` gives `CONFIRMED`;
- primary `p >= 0.05` gives `NOT_CONFIRMED`;
- the strict 20% analysis is sensitivity only and cannot rescue the 10% primary result.

## Qualification requirements

A head is first-pixel-ready only when the P500 execution-gate run on that exact SHA succeeds. The qualification must pass all of the following without opening a P500 pixel:

- execution-gate unit tests;
- both workflow YAML contracts, including the terminal `h2-prospective-test` job;
- frozen candidate metadata and authorization checks;
- unopened future-output namespace checks;
- exact inherited blind-measurement interface materialization;
- prospective execution-script compilation;
- synthetic PASS and NOT_EVALUABLE support-receipt -> H2 contract tests;
- four-case legacy H2 replay with exact 999-draw legacy/prospective structured-null equality;
- the 499-species / 49,900-row / 256-partition firewall dry run;
- qualification receipt and artifact creation.

The qualification artifact, not a hand-written SHA in this document, is authoritative for the exact qualified head. This prevents documentation edits or unrelated branch movement from being silently treated as preopening evidence.

The four frozen replay anchors are:

| threshold | cohort | N | W | structured-null p |
|---|---|---:|---:|---:|
| 0.10 | discovery | 152 | 0.5146245742375285 | 0.001 |
| 0.10 | reserve | 129 | 0.5145857778624213 | 0.001 |
| 0.20 | discovery | 75 | 0.5423549779816461 | 0.001 |
| 0.20 | reserve | 65 | 0.5105174254190324 | 0.008 |

For all four cases the rebuilt legacy numerical summary must match the committed frozen result and the full 999-value legacy/prospective null vectors must match exactly.

## Head-stability guard

The legacy 34-species figure workflow previously had permission to auto-commit generated figures on any pull-request branch and therefore moved the P500 head after qualification. Its `generate-figures` job is now restricted to branches beginning `analysis/34species-paper-`. It cannot auto-commit on the P500 execution branch.

## Claim boundary

Qualification CI is technical evidence only. It does not open P500 images, run P500 ROI/segmentation, measure P500 colour, calculate prospective D, or test prospective H2. No prospective biological conclusion exists until a successful PREOPENING receipt is issued on the exact qualified head and consumed by the location-blind measurement workflow.
