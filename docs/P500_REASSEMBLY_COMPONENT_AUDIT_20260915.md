# P500 reassembly: synthetic component audit

Executed successfully on 2026-09-15 using only 49,900 artificial IDs and
499 artificial species. No P500 result artifacts or pixels were opened.

Reproduce from a checkout containing the recorded historical Git object:

```
python scripts/analysis/audit_p500_reassembly_synthetic.py
```

The script refuses a different source blob, performs no fetch, and loads
`fcp_pipeline/photo_first_measurement_execution.py` from measurement source
`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`, blob
`b6e76efb4db8e7529b5ce5a91c60b6ae68637fb5`. This is the source materialized by
the active P500 workflow, not a replacement implementation.

## Observed behavior

- All 49,900 artificial acquisition failures in 256 partitions survive joining
  as mixed_uncertain; zero are counted as classified biological states.
- Missing partition, missing terminal ID, duplicate terminal ID, foreign
  terminal ID, duplicate worker ID, duplicate metadata ID and missing metadata
  ID each raise the expected error before returning a joined result.
- An unknown terminal status is rejected by the partition sealing validator.
- The aggregate helper itself accepts an unknown status when that prior
  validator is bypassed. This characterizes reliance on the earlier sealing
  stage, not evidence that the active run contains invalid statuses.

## Remaining verification

This exercise tests the inherited helper and sealing validator, not the full
P500 wrapper, artifact transport, receipt filenames, source hashes at runtime,
partition membership, palette validity, support gate or H2 statistics.
Synthetic success is not evidence that real measurement has completed.

Source inspection of the P500 wrapper at
`f403606b5611d895988a6a4f15d9ca1e281ce848` shows exact partition-key set checks
and an aggregate receipt-row total, but not a per-receipt CSV-content digest
comparison. The sealing receipt records counts rather than a CSV hash.
At completion, independently verify artifact archive identities, each
receipt's row counts against its CSV, deterministic membership, allowed
status/morph pairs, and the full unique frozen-ID universe. Do not use a
successful aggregate job alone as evidence for all of these properties.

No active workflow, frozen source, result or admission rule was changed.
The [protocol compatibility limits](P500_PROTOCOL_COMPATIBILITY_AUDIT_20260915.md)
remain separate from these data-integrity checks.
