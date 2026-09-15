# P500 execution and measurement-control compatibility

Scope: source audit, not a colour-result analysis. Inspected execution revision:
`f403606b5611d895988a6a4f15d9ca1e281ce848`, run
[34919485994](https://github.com/zuizui0223/fcp/actions/runs/34919485994).
No partition contents, measured palettes or prospective statistics were read.
The run remains active; final compliance is not established by this audit.

## Decision

The executing location-blind protocol does **not establish satisfaction of the
older white measurement-control contract**. Its eventual `CONFIRMED` label is
an H2 statistical decision, not `CLEAR` under that control contract. Keep the
two decisions separate. No biological result or computed artifact verdict is
issued here. The original control receipt is unchanged.

| Requirement | Inspected execution | Interpretation |
| --- | --- | --- |
| Fixed population, no replacements | Authorization checks 499 species / 49,900 rows and frozen metadata identities | Design alignment; terminal completeness still needs artifact verification |
| Fixed ROI, palette and H2 axis | Materializer fixes the inherited measurement source; runner fixes q_white, thresholds, seeds and 999 null worlds | No evidence here of outcome-driven retuning; not evidence of focal-petal validity |
| Species/location hidden during measurement | Worker surface is restricted; acquisition surfaces are removed before measurement | Location-blind design, not response-blind technical diagnostics |
| All terminal records before metadata join | Workflow requires all 256 partitions and invokes reassembly before H2 | Verify receipts, unique IDs and failure accounting at completion |
| Technical highlight table before response join | No separate technical-table construction/freeze stage in the inspected workflow | Required historical stage is not established |
| Freeze high_clip membership before colour-response join | No high_clip freeze stage in the inspected workflow | Cannot reconstruct its historical ordering after the join |
| Within-species white/near-clip coupling and OR interval | H2 runner performs q_white/W tests, not the conditional coupling analysis | H2 significance does not answer the artifact-coupling question |
| High-clipping exclusion H2 sensitivity with >=90% retention | Runner compares 0.10 and 0.20 morph thresholds only | The strict morph threshold is not the high-clipping sensitivity |
| Control CLEAR / FLAGGED / INDETERMINATE before interpretation | No control evaluator is invoked before the H2 decision | Do not translate the H2 label into a control clearance |

## Source evidence and boundaries

The [older control protocol](POLYMORPHISM_H2_P500_WHITE_MEASUREMENT_CONTROL_PROTOCOL_20260913.md)
sections 4, 6 and 8 require the response-blind exclusion set, retained-species
sensitivity and ordered freeze. The actual
[workflow](https://github.com/zuizui0223/fcp/blob/f403606b5611d895988a6a4f15d9ca1e281ce848/.github/workflows/polymorphism-h2-p500-prospective-measurement.yml)
runs acquisition, ROI/palette measurement, partition sealing, reassembly and
the [prospective H2 runner](https://github.com/zuizui0223/fcp/blob/f403606b5611d895988a6a4f15d9ca1e281ce848/scripts/analysis/run_polymorphism_h2_p500_prospective_white_axis_20260915.py).
The runner imports the existing geometry/structured-null helpers, computes
primary and strict morph-threshold tests, and emits its own H2 verdict.
It does not invoke the publication branch's snapshot or coupling helpers.

This inspection does not prove that every possible technical diagnostic is
absent from every inherited output. Even if a diagnostic column is found in
the complete aggregate later, its existence alone would not prove that a
response-blind exclusion membership was frozen in the required order.
The workflow also deletes ephemeral image pixels and masks; do not assume
missing diagnostics can be recovered without a new measurement operation.

## Handling the eventual aggregate

1. Verify the complete run and immutable aggregate identities; distinguish an
   artifact-upload or Git-commit failure from a scientific non-support result.
2. Verify 256 terminal receipts, the full 49,900-ID universe, uniqueness and
   retained non-classified states before scientific interpretation.
3. Read support gates before interpreting any `primary_confirmed=false` flag:
   measurement/vector insufficiency means not evaluable, not a negative test.
4. Preserve the exact primary and strict results. A strict result cannot rescue
   primary non-support. Neither threshold is the missing highlight sensitivity.
5. Report any support as conditional on the frozen photo-derived measurement
   and admission design, with unresolved highlight control and the failed
   Monarda agreement gate. Do not claim artifact-cleared biological replication.

No restart, pixel reacquisition, retroactive pre-opening receipt, new exclusion
set or changed test is authorized here. A later exploratory diagnostic would
require explicit labeling and cannot repair historical pre-opening compliance.
For the actual opening chronology see the
[Actions correction](FCP_H123_P500_ACTIONS_CORRECTION_20260915.md).
