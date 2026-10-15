# FCP PR #200 — distributed original-photo production handoff (2026-10-15)

**Status at authoring:** 42,111-taxon source census [Actions 38676167958](https://github.com/zuizui0223/fcp/actions/runs/38676167958) is running, **not finished**. Zero full-study images have passed final assembly. The separate 96-image pilot is not biological full-census input. This file is an operating procedure, NOT an all-species finding.

## Admission before pixels
A complete dated 2026-09-27 source is mandatory. Accept only a retained run artifact whose full_source_execution_receipt.json says COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE with complete_metadata_census=true, 4 pinned original metadata SHA256s and 42,111 complete taxon checkpoints. Its frozen_location_blind_all_photo_firewall/fresh_study_freeze.json must say FRESH_BLIND_METADATA_FROZEN, source_species=42111, worker_packets=128, and carry the full distinct photo-ID worker ledger. Source photo total is unknown until that exact source audit finishes. Do not reuse old 149,900 photo colour outcomes.

## Automated admission after source completion

Download **both** the successful source-status artifact and complete-source
tar.gz artifact from Actions #38676167958 onto durable storage, then unpack
the complete source archive without reselecting taxa. The status artifact
carries the receipt, dated archive pin and 42,111-row photo-count ledger;
the complete-source artifact carries all manifests and the blind firewall.

Run this **non-image** admission check before ANY source-photo requests:

~~~bash
python scripts/acquisition/admit_fcp_complete_photo_source.py \
  --source-package-root /persistent/fcp/unpacked_source \
  --source-status-receipt /persistent/fcp/source_status/full_source_execution_receipt.json \
  --original-pin /persistent/fcp/source_status/fcp_source_archive_pin.json \
  --source-counts /persistent/fcp/source_status/all_42111_species_source_photo_counts.csv \
  --out /persistent/fcp/new_measurement_workplan --batch-size 64
~~~

The tool refuses incomplete source receipts, changed dated archive/table hashes,
missing or duplicate species, corrupted manifests, inconsistent photo counts
and blind-worker hash drift. Its source_admission_receipt.json must read
FULL_SOURCE_PHOTO_WORK_PLAN_READY_NO_PIXELS_YET with **42,111 taxa and 128
source-blind worker partitions** plus the exact all-photo denominator.

Admission is not execution. The receipt reports a theoretical minimum
single-stream duration based solely on 1.1 seconds between new HTTP
attempts, NOT an actual measurement ETA or an endorsement of large-scale
provider downloads. Dedicated durable compute and an appropriate
provider-approved global access plan are still necessary. Never re-download
a 36-GB source archive merely because the first run is still in progress.

## Workers: partition-local output, never concurrent shared output
**Provider rate is a global constraint.** The image loader now enforces ≥1.1 seconds between actual S3 attempts (including retries and across work-unit boundaries) within a process. For multiple cooperating processes on the **same Linux/POSIX-locking storage with synchronized clocks**, set the **same absolute** `FCP_PROVIDER_GLOBAL_RATE_FILE` for every worker; its locked timestamp also covers separate invocations. Without that variable, only process-local pacing is guaranteed. This does **not** globally coordinate geographically separate hosts, third-party downloads or a public-provider-approved quota. Do not start 128 concurrent public S3 downloaders without independently approved aggregate access, pacing and persistent compute/storage. An unavailable/malformed rate lock **or persistent source-wide HTTP 429/5xx/401/403 or connection outage** aborts the current work unit without producing sealed terminal records. This is not a reason to mark the frozen source photos permanently unclassifiable. Only a verified per-photo terminal condition (e.g. original source HTTP 404/410, zero/oversize bytes under the declared 12MiB cap, or undecodable returned photo bytes) can receive an explicit unclassified terminal without replacing its photo ID. **Unexpected exceptions, invalid provider URLs, and local photo-write/disk errors now abort the work unit unsealed rather than masquerade as photographic missingness.** The original GitHub Actions source runner is ephemeral; an artifact must be moved to durable storage before heavy ROI work.

For one partition, after materializing the original source + detector/EfficientSAM SHA-verified model assets, invoke:

~~~bash
F=/persistent/fcp/frozen_location_blind_all_photo_firewall
ROOT=/persistent/fcp/partition_outputs
# Same absolute path in every cooperating worker process; POSIX lock required.
export FCP_PROVIDER_GLOBAL_RATE_FILE=/persistent/fcp/provider_s3_rate.lock
bash scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --full

python scripts/acquisition/execute_fcp_full_photo_work_queue.py \
  --firewall "$F" --output "$ROOT/partition_000" \
  --partition 0 --batch-size 64 --max-work-units 32 \
  --original-measure-script scripts/analysis/measure_random_photo_first_blind_partition.py \
  --original-detector data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt \
  --original-efficient-sam /persistent/models/fcp-original-efficient-sam
~~~

Repeat with identical source/method bytes and partition-local output, honoring one globally coordinated request budget. The 32-work-unit invocation budget and 64-photo compute batch size never cap biological inclusion. **If a provider outage interrupts a unit**, its staging folder remains under `unsealed_inflight_batches`; the next execution deliberately rejects that folder. Inspect and explicitly quarantine the unsealed attempt before restarting the same frozen partition/batch, retaining its diagnostic evidence. Never move incomplete CSV/JSON into `sealed_terminal_batches`, relabel 429 photos as missing, or discard source IDs.

## Safely join distributed sealed results
Only while partition workers are quiescent:

~~~bash
python scripts/acquisition/collect_fcp_partition_measurements.py \
  --firewall "$F" --partition-output-root "$ROOT" \
  --output /persistent/fcp/collected_all_partitions --batch-size 64
~~~

The collector rejects wrong partition file names, method drift, orphan CSV/JSON receipts, conflicting copies or source hashes, and verifies every source-bound partition/work unit with the existing uncapped planner. It never opens images, geography or historical measurement outcomes. It may be rerun to pick up newly sealed results from all 128 disjoint partition outputs.

Its partition_collection_progress.json distinguishes:
- PARTIAL_ORIGINAL_PHOTO_TERMINALS_NOT_SCIENTIFIC_RESULT: work units remain
- ALL_ORIGINAL_PHOTO_TERMINALS_SEALED_READY_FOR_REASSEMBLY_NOT_H2: 0 pending; independent assembly still required

At zero pending, run full source/photo/method reconciliation:

~~~bash
python -m fcp_pipeline.fcp_fresh_study reassemble \
  --firewall "$F" \
  --sealed-results /persistent/fcp/collected_all_partitions/sealed_terminal_batches \
  --method-attestation /persistent/fcp/collected_all_partitions/sealed_terminal_batches/frozen_method_attestation.json \
  --outdir /persistent/fcp/full_new_photo_results
~~~

Only fresh_complete_receipt.json with status COMPLETE_FRESH_MEASURED_PHOTO_CENSUS and identical source/photo terminal denominators supports the statement that original-photo measurement is finished. All original photographs receive a terminal status (including strict QC failures), while taxa with fewer than 40 classifiable photos remain insufficiently supported rather than silently excluded.

Fixed W/H2 coarse-state-preserving 999-world inference and direct per-photo climate/soil comparison are separate later gates. A full terminal census is not an H2/ecological conclusion.
