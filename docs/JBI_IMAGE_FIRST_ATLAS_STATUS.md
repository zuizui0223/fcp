# FCP image-first atlas execution status

## Final terminal decision

The prospective 200-species / 60,000-photo image-first atlas has reached its frozen terminal stopping rule.

**Final decision: `not_evaluable`.**

The exact terminal run is GitHub Actions run `33592829701`, attempt 2, at frozen execution head `aea19a4eff9585f501aa6a833ad44c80080eddcf`. The workflow concluded `failure` because the protected coordinate join was deliberately fail-closed after the measurement gate returned `not_evaluable`; this is not a failure of the exact compute-partition reassembly.

The exact terminal artifact is `jbi-atlas-terminal-measurement-v5` (artifact `9857419176`, SHA-256 `6bfff229a90215d016b2dd7e2bcca9446474f6d216efe33782f9faac11e69b55`). A compact immutable receipt is committed at `docs/supporting/jbi_atlas_terminal_measurement_v5_receipt.json`.

## What completed successfully

The location-blind measurement execution completed the frozen scale-out without replacement, favourable replicate selection or early stopping:

- 256 / 256 compute partitions were present and validated;
- the exact union contained 60,000 unique measurement IDs;
- reassembly produced 16 semantic shards;
- no candidate image pixels were persisted in the reassembly evidence;
- measurement did not open coordinates;
- the reassembled ROI-v4 bundle contained exactly 60,000 terminal records;
- estimator identity, trained-weight hash and ROI-contract hash matched the frozen v5 measurement contract.

`compute_partition_reassembly_v1.json` therefore passed as `pass_exact_256_compute_partition_coverage`, and the reassembled bundle passed as `pass_complete_location_blind_roi_v4_measurement_v5_bundle`.

## Binding measurement gate

The frozen `measurement_gate_v5.json` returned:

`not_evaluable_scaleout_measurement_completeness`

Only **58 / 200 species** met the frozen measurement-evaluable rule. Evaluable-species counts by the eight fixed cohorts were:

| Cohort | Evaluable species |
|---|---:|
| C01 | 3 |
| C02 | 4 |
| C03 | 7 |
| C04 | 8 |
| C05 | 8 |
| C06 | 12 |
| C07 | 7 |
| C08 | 9 |

Every cohort was therefore `not_evaluable` under the already-frozen cohort requirement.

The binding state is:

- `coordinate_join_permitted = false`;
- `coordinates_opened = false`;
- `superseded_v3_ordered_inference_used = false`;
- claim ceiling: measurement completeness only.

## Cascade closure

The confirmatory sequence was frozen as:

`species-conditioned spatial organization -> shared transition -> environmental concordance -> pollinator biogeographic concordance`

Because `not_evaluable` never advances the cascade, the terminal experiment stops **before species-conditioned spatial organization**. The protected coordinate-colour join is not opened. Shared-transition, environmental-concordance and pollinator-concordance analyses are not run.

No threshold, seed, metric, estimator, source, denominator, cohort requirement, branch order or `not_evaluable` stopping rule is altered after observing the terminal measurement result. The 60,000-record denominator remains fixed.

## Immutable earlier evidence

The terminal result does not reopen or overwrite completed evidence:

- six-species Chapter 1: 1,200 photographs; Stage A `p = 0.0113`, Stage B `p = 0.0906`;
- 34-species / 25-family literature comparison: separate comparative method and supplementary context;
- three-species automated negative validation: 717 locked photographs, 306 admitted encounters, zero supported species;
- 50-species sentinel geometry: 20,200 observations, retained as an unopened precursor and never pooled into the terminal experiment.

The terminal atlas is therefore **not a biological replication result**. It is a prospectively executed measurement-scale validation that became non-evaluable under its own frozen completeness rule.

### Historical automated-colour pilot STOP

For audit compatibility, the earlier automated-colour pilot checkpoint is preserved verbatim as a historical state: **bulk atlas image opening is stopped**. That sentence describes the pre-scale-out STOP that was binding before independent ROI-v4 qualification, dated-source/source-role resolution and final preimage authorization subsequently passed. It is not the current terminal state and does not imply that terminal candidate pixels remained unopened; the terminal run later opened pixels location-blind and ended at the measurement-completeness gate described above.

## Manuscript role

For the Journal of Biogeography package, the six-species held-out Stage A/B analysis remains the primary biological result. The terminal 200-species scale-out should be reported as a transparent prospective extension that completed exact 60,000-record location-blind measurement and then stopped at the predeclared completeness gate. It must not be described as evidence for or against spatial organization, shared transitions, environmental concordance or pollinator concordance.

The scientific consequence is narrower but useful: increasing image count and taxonomic breadth did not by itself guarantee an evaluable cross-species colour field under the frozen ROI-v4 completeness requirements. Any future atlas attempt must be a new prospective experiment with its measurement model and qualification rules fixed before new outcome pixels are opened; the present 60,000 records cannot be mined with substituted post-result gates.

## Durable evidence

- terminal receipt: `docs/supporting/jbi_atlas_terminal_measurement_v5_receipt.json`;
- frozen measurement contract: `docs/supporting/jbi_atlas_measurement_execution_contract_v5.json`;
- frozen real-colour inference amendment: `docs/supporting/jbi_atlas_real_colour_inference_amendment_v5.json`;
- frozen terminal follow-up implementation amendment: `docs/supporting/jbi_atlas_terminal_followup_implementation_amendment_v1.json`;
- exact run: `33592829701`;
- exact artifact: `9857419176` / `jbi-atlas-terminal-measurement-v5`.
