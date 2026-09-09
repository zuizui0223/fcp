# FlowerMask: fixed 300-ID paired-metadata feasibility gate

Prospective plan: **8 September 2026**, before any new paired-folder request.
This implements the next gate in the
[measurement-target design](RGFCA_MEASUREMENT_TARGET_AND_VALIDATION_V2.md).
It authorizes a bounded public **metadata** inspection after preflight, not
images, masks, annotation bodies, model execution or a new ecological test.

## Frozen scope

- Frame: the 300 annotation IDs already exposed in the earlier public metadata
  audit, with six schema-opened IDs explicitly retained. No selection by colour,
  location, model performance, size or successful download.
- Correspondence rule: within each exact provider group, map `N.json` to the
  `N_dataset` child of its masking container. All 300 mappings are uniquely
  reconstructible from the saved 319-folder response. This is **documentary
  filename correspondence**, not proven image/mask alignment.
- Inputs: exact original annotation plan SHA-256
  `48436a74160b4700e2d9b2bbc803aeac7fcfe8a1f09885140a1eade034c5aade`,
  folder response `fd30ad1ba248c43eca04c1a53cec31e87617dbd815b5a98b752eb461dec2b892`,
  and cached Butterfly Pea folder metadata
  `a227790321584118998d390850cca1b66ab4f7ac54276ebc0c0262f566c77485`.
  The latter two exact-byte public records are retained in
  `data/validation/flowermask_pair_metadata_inputs_v1/`; no image payload is there.
- New [plan](supporting/rgfca_flowermask_pair_metadata_plan_v1.json) SHA-256:
  `2ce2a578db6c470a9021189241e2e5c9fae561c2435896f64e61d50e0ccd4d83`.
  It contains all 300 rows: one cached response and **at most 299 new requests**.

## Access and preservation rules

Request only the exact version-1 public file-list endpoints in the plan. Start
requests at least one second apart; use a 30-second timeout and a 250,000-byte
response cap. No retries, redirects, pagination follow-ups, alternate URLs,
author contact, authentication, paid access or changed rights. Never follow
returned download/view URLs. These endpoints return file metadata; the program
rejects non-JSON and non-list responses and an unexpected annotation/imageData
schema. Only valid metadata response bodies are retained. Failed responses
retain byte count/hash, when available, and error type, not opaque payloads.

Before the first request, exclusively claim a new output directory and durably
write all 300 pending rows. Checkpoint request starts and every processed row
with atomic, flushed replacement. On the first transport/schema failure stop
requests, retain that failure and mark all remaining rows not attempted after
failure. No automatic resume or replacement. A killed execution retains its
last checkpoint; it is not counted complete or restarted automatically.

Within valid metadata, preserve each filename, file ID, advertised size/hash
and content type. Detect duplicate IDs/names, casefold collisions, wrong folders
and conflicting sizes. The bounded pairing indicator requires exactly named
`img.png`, `label.png`, `label_names.txt` with PNG/PNG/text content types.
Missing members yield an unresolved pair indicator, not a corrected or excluded
row. Advertised payload hashes are **not verified image hashes**. A returned
listing does not establish full pagination or complete 3,600-image access.

## Execution sequence

1. Commit this protocol, the exact plan/input receipts, code, artificial tests
   and workflow before acquisition. Push preflight reconstructs the whole plan
   and tests failures without requesting any provider resource.
2. After that preflight succeeds and its artifact/source identity is checked,
   dispatch the same qualified commit with the exact plan hash. A read-only
   GitHub history gate requires one first-attempt dispatch for this v1 workflow;
   prior/concurrent dispatches, retries or incomplete history stop execution.
   Fixed workflow concurrency prevents parallel v1 acquisition. The saved
   result's later presence closes execution at subsequent heads.
3. Preserve all terminal rows and valid response bytes, even for non-support or
   source failure. Recompute metadata summaries from saved bytes; do not
   reacquire for verification. Record code/plan/run identity and actual request
   counts. CI success means faithful metadata inspection, not model validation.

Implementation: `scripts/analysis/audit_rgfca_flowermask_pair_metadata.py`,
`tests/test_rgfca_flowermask_pair_metadata.py`, and
`.github/workflows/rgfca-flowermask-pair-metadata.yml`.
The raw response grammar and all scope/rate/denominator rules are fixed before
these requests. A source/API incompatibility is retained, not silently repaired
inside this acquisition.

## Claim ceiling and following gate

The resulting table can identify publicly advertised candidate image/mask
pairs for these 300 IDs only. It cannot establish annotation completeness,
anatomical inclusion rules, focal-taxon attribution, pixel alignment, colour
accuracy, observer/event independence, model-training independence or an
ecological signal. No benchmark is admitted by a positive pairing indicator.

If documentary pairing is complete, the next possible gate is separately
planned reference correspondence/annotation-target inspection with a declared
payload budget, development/validation status and rights scope. If pairing or
access fails, retain the precise unresolved outcome and do not score a
convenient subset. The closed six-document FlowerMask inspection, Monarda
non-support and all completed FCP ecological tests remain unchanged.

## Completed result: metadata inspected, paired frame incomplete

The one authorized [run 34201981494](https://github.com/zuizui0223/fcp/actions/runs/34201981494)
completed on **8 September 2026 at 08:03:35 UTC**. Its acquisition code and plan
were fixed at `49d24a1b3ca060608d629577a01ce337b2bb45fe`, after successful
push preflight `34201856801` and exact-byte verification of its nine qualification
files. The metadata job itself ran for 5 minutes 14 seconds. It inspected all
300 rows, reused the one cached response and started exactly 299 new requests.
No transport/schema failure occurred. All six previously schema-opened IDs and
all six provider groups of 50 rows are retained, with provider spellings unchanged.

**299/300 returned listings advertise the required three pair members;
`all_300_pair_metadata_present=false`.** Butterfly Pea `33.json`, annotation ID
`bf1604ee-fe76-4ec2-b2b7-51c50a732090`, maps to `33_dataset`, folder
`d88759cb-dfb0-4996-8b33-d1002cdc33eb`. That returned listing has `label.png`,
`label_names.txt` and `label_viz.png`, but no `img.png` entry. Its exact response
is 2,251 bytes, SHA-256
`a3246f20b3eab5436a461a7e9d7d6609f2f299e74a3390ad741fae2f2f75a83f`.
This establishes an absence **in the inspected listing**, not that an original
image does not exist elsewhere. Pagination and alternate source paths were not
followed. The visualization is not substituted for the missing original.

The fixed inspection completed successfully, but its full-pairing indicator
did not pass. The unresolved row is neither dropped nor replaced. There is no
299-image convenience benchmark, retry, hidden endpoint recovery, new model run
or ecological conclusion from this outcome. The next reference-acquisition or
correspondence design must be specified separately; it is not authorized here.
Annotation scope, rights, actual alignment and development/validation status
would still require evidence even if every filename were present.

### Retained evidence and offline reconstruction

The [exact 303-file bundle](../data/validation/flowermask_pair_metadata_result_v1/)
preserves the 300 valid metadata response bodies, full terminal result,
execution-start receipt and dispatch history from artifact `10046355184`.
It totals **1,550,546 bytes** and contains no image or mask payload. These are
native byte copies, not regenerated responses. The result SHA-256 is
`27c79f60b4088c4c2942778c43c9089509ac3dcce5ec825166d97804166359b3`;
execution-start SHA-256 is
`b133ba62b29a0ed5539b28ed54218b5e74da0069d669aec8598b7f101aa08793`;
dispatch-history SHA-256 is
`cc1ba4ef73ec7f7a57feffe4907cc591ab0782f57e94040436993e70999f3da0`.
The provider reports artifact-ZIP SHA-256
`a26dabf38c4a51274303dad0191fb2284df65b0806832f9c988210018df950e6`;
the compressed ZIP digest was not locally checked. Extracted response bytes,
identities, summaries and receipt hashes were checked directly.

The [retained summary](supporting/rgfca_flowermask_pair_metadata_result_v1.json)
also closes subsequent acquisition heads. The additional retained-result
verifier pins the completed run's receipts, requires the exact 303-file census,
checks the first-attempt dispatch and reconstructs all 300 summaries using the
original frozen metadata parser. It does not rerun acquisition or modify that
parser. Corruption tests and a separate read-only retained-result CI preserve
the incomplete-pairing outcome as a valid, reproducible result.

```text
python -m scripts.analysis.verify_rgfca_flowermask_pair_result data/validation/flowermask_pair_metadata_result_v1 --verify-summary docs/supporting/rgfca_flowermask_pair_metadata_result_v1.json
```

CI success certifies faithful preservation, not pairing completeness or an
admitted measurement benchmark. Model execution, decoded image pixels and
coordinate joins all remain false for this gate.
