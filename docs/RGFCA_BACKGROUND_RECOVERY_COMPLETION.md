# Discovery matched-background recovery: completed, not evaluable

Verified 7 September 2026. This is an execution and evidence-quality audit of the
closed **discovery** recovery, not a new ecological analysis. Reserve outcomes
were not read, images were not reacquired, and the original recovery/inference
code, masks, thresholds, seeds and matched denominator were not changed.

## Result and claim ceiling

All 128 recovery partitions and finalization succeeded in
[run 34094784607](https://github.com/zuizui0223/fcp/actions/runs/34094784607),
executed at `3ba0c11a66b1b0c6aef09483ebc7eb0666f244a3`.
The workflow completed successfully because its predetermined failure outcome
was correctly recorded; green CI does not mean the scientific control passed.
The [immutable summary](https://github.com/zuizui0223/fcp/blob/2f00847ccc4d15be8637a1e3589239ce2ab5c114/docs/supporting/global_rgfca_background_falsification_result_v2.json)
reports `not_evaluable_incomplete_exact_background_recovery`.

| Saved reproduction check | Rows failing / 21,424 |
|---|---:|
| Original-image SHA agreement | 0 |
| ROI admission reproduced | 0 |
| Flower-mask pixel count | 79 |
| Full flower-palette counts | 79 |
| Background pixel count | 6 |

The 79 flower-count failures are the same records; the six background-count
failures are separate. Thus **21,339 exact + 85 failed = 21,424 retained records**.
No background-adjusted statistic, null distribution or p-value was computed.
No replacement, threshold relaxation, imputation or successful-subset inference
is permitted. The failure is not a negative ecological result and does not show
that background confounding is absent or present. The original flower-only
rho = 0.0270213, p = 0.001 is unchanged but remains an exploratory photo-derived
association, without completed flower-specific validation.

These checks inspect saved worker flags and byte/ID provenance, not newly
decoded images. The saved fields identify which equality checks failed but do
not identify the root cause or retain mismatch magnitudes. Equal mask pixel
totals and palettes would not prove bitwise mask identity or correct focal-taxon
attribution. No such claim is made here.

## Artifact and complete-census verification

- Final artifact: `10017470032`, `global-rgfca-background-falsification-full-v2`,
  807,701 bytes; GitHub-reported archive digest
  `1a63d4ab22855d614d20d4569ed9f65cca66f845d4e3d4e1a4cb2b6eebe4e4d4`.
- The downloaded final artifact contains only the summary and full terminal CSV,
  not species statistics or null outputs. Summary bytes equal result commit
  `2f00847ccc4d15be8637a1e3589239ce2ab5c114`.
- All 128 named partition CSVs and their manifests were downloaded from the
  completed run. Their rows, in finalizer filename order, exactly reproduce the
  final CSV, including every failure. No missing or duplicate partition is used.
- All 21,424 unique measurement IDs match the original frozen eligible
  369-species frame. Original contract, timing amendment and measured-table
  hashes match the summary's lineage using exact Git bytes, not platform-newline
  substitutions.
- Terminal CSV SHA-256:
  `4ac8b8edbaaf3c784ee7648367f50b994e1cde7114f9cc9d4ad5434ee9d56ae4`.
- Machine-readable verification:
  [completion audit](supporting/rgfca_background_recovery_completion_audit_v1.json).

Reproduce without images or statistical inference:

```powershell
gh run download 34094784607 --repo zuizui0223/fcp --name global-rgfca-background-falsification-full-v2 --dir .artifacts/background-34094784607/final
gh run download 34094784607 --repo zuizui0223/fcp --pattern 'global-rgfca-background-recovery-*' --dir .artifacts/background-34094784607/partitions
python scripts/analysis/audit_rgfca_background_recovery_completion.py --full-dir .artifacts/background-34094784607/final --partitions-dir .artifacts/background-34094784607/partitions
```

The focused regression tests execute the exact frozen finalizer's stop branch
with artificial status-only inputs: one failure or 85 failures must stop before
opening the statistical frame; missing partitions and duplicate IDs must fail.
These are software checks, not synthetic ecological evidence.

## Remaining work

The original integer-parser failure in run `34091174522` and its
[technical continuation receipt](supporting/global_rgfca_background_control_numeric_recovery_v2b.json)
remain preserved. The new terminal reproduction failure is distinct from that
earlier exception; the integer fix is not repeated or blamed without evidence.

The separate prospective reserve run `34091091640` measures flower and background
on the first decode. It continues unchanged and is not a repaired version of this
discovery recovery. All 256 partitions and 50,000 terminal records must be
verified before its fixed four tests. Their outcomes, target-domain measurement
validity and the final submission package remain outstanding.
