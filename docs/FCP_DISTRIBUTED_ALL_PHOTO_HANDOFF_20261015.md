# FCP PR #200 — distributed original-photo production handoff (2026-10-15)

**Execution status (15 October 2026):** The first 42,111-taxon source job [38676167958](https://github.com/zuizui0223/fcp/actions/runs/38676167958) **failed** while writing full photo manifests (`OSError: [Errno 28] No space left on device`) after successfully verifying the 36-GB dated archive and all four metadata tables. The corrected disk-reclamation recovery [38699975166](https://github.com/zuizui0223/fcp/actions/runs/38699975166) began at 07:37 UTC and remained **in progress at 09:50 UTC**; no complete source artifact or new full-study ROI-v4 measurement was available at that check. Note that a later indexed SQLite streaming optimization was committed **after** recovery started, so its improvement is not part of this currently running recovery. The separate 96-image pilot is not biological full-census input. This file is an operating procedure, NOT a completed all-species finding.

## Admission before pixels
A complete dated 2026-09-27 source is mandatory. Accept only a retained run artifact whose full_source_execution_receipt.json says COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE with complete_metadata_census=true, 4 pinned original metadata SHA256s and 42,111 complete taxon checkpoints. Its frozen_location_blind_all_photo_firewall/fresh_study_freeze.json must say FRESH_BLIND_METADATA_FROZEN, source_species=42111, worker_packets=128, and carry the full distinct photo-ID worker ledger. Source photo total is unknown until that exact source audit finishes. Do not reuse old 149,900 photo colour outcomes.

## One-shot recovery and real first-photo auto-handoff (15 October 2026)

The first full snapshot census failed while writing species photo manifests, after the
immutable 36GB archive and four source tables were SHA-verified, due to exhausted
GitHub Actions scratch. The disk-reclaimed second source run is
[38699975166](https://github.com/zuizui0223/fcp/actions/runs/38699975166).
Do not treat its running state or CI test passes as a real source receipt.

There are **two disjoint, source-gated paths** already wired on `main`:

- If recovery #38699975166 succeeds, the [first original-photo automation](https://github.com/zuizui0223/fcp/blob/main/.github/workflows/fcp-pr200-first-real-photos-when-source-completes.yml) accepts the exact run SHA and both source artifacts, then executes one 64-photo genuine ROI-v4 work unit.
- If recovery #38699975166 **fails**, the [indexed source fallback](https://github.com/zuizui0223/fcp/blob/main/.github/workflows/fcp-pr200-fallback-indexed-source.yml) downloads its status artifact *first* and refuses another 36GB transfer unless the archive/metadata hashes are correct, all 42,111 original taxa were targeted, and the failure is in the source photo-manifest scan after all four metadata tables were extracted. The fallback uses pinned source code `324e5b0926001272e8864d4517fc198c639e0e24`, which adds the indexed SQLite streaming scan absent from the already running second attempt. Only successful fallback artifacts are eligible for the **same** first real photo automation.

Neither path autoruns the remaining uncapped millions of source photos without
durable capacity and coordinated lawful provider access. No all-species
measurement, fixed W/H2 test or direct soil/climate association may be claimed
from a metadata-only receipt or a 64-photo first work unit.

## Split the SOURCE census itself, not just the image model

The global 42,111-taxon study no longer has to re-read the original ~36GB dated
iNaturalist archive for every source photo-manifest shard. In the **new opt-in
durable-storage path**, scan all four SHA-pinned provider tables into a
disk-backed, globally deduplicated SQLite source **exactly once**, preserve
that database with an immutable SHA256 receipt, and create up to 128 complete,
disjoint metadata/photograph-ID shards in small resumable batches.

On a **persistent disk** with sufficient real free capacity, and a local
trusted original archive (or one deliberately authorised download), run:

~~~bash
ROOT=/persistent/fcp/complete_42111_durable_source
FRAME=/persistent/fcp/42111_species.csv
ARCHIVE=/persistent/fcp/inaturalist-open-data-20260927.tar.gz

# A. Global observation/photo ID join once, preserve exact SQLite bytes.
python scripts/acquisition/execute_fcp_42111_census.py \
  --frame "$FRAME" --archive "$ARCHIVE" --work "$ROOT" \
  --expected-species 42111 --shards 128 --worker-partitions 128 \
  --stage-only --release-verified-temporary-tables

# B. SHA256/integrity-check the huge SQLite source ONCE, emit and seal
#    128 independent source shards in this same process. Each shard gets its
#    own complete receipt; failure leaves previously completed shards intact.
python scripts/acquisition/emit_fcp_staged_source_shards.py \
  --frame "$FRAME" \
  --staged-metadata "$ROOT/complete_42111_metadata_source" \
  --shards 128 --expected-species 42111 \
  --shard-start 0 --shard-stop 128

# An interrupted durable run can be resumed by requesting the same range.
# Each previously sealed shard is independently hash-checked and not rewritten.
# A partially written unsealed shard must be inspected/quarantined first.

# C. Audit all shards, only then freeze all 128 blind photo workers and
#    pack the portable complete metadata source WITHOUT the huge scratch DB.
python scripts/acquisition/finalize_fcp_staged_source.py \
  --frame "$FRAME" --work "$ROOT" \
  --expected-species 42111 --shards 128 --worker-partitions 128
~~~

**Important:** This is genuine *source-side* splitting; original image
measurement is a further, independently partitioned phase. A completed source
shard is never a complete flower-colour result. A failed partially written
shard is quarantined for review rather than overwritten, while previously
sealed shards remain intact. The original SQLite source requires durable,
sufficiently large storage: this change does **not** make the initial full
provider scan fit on an undersized ephemeral CI disk and does **not** authorize
multiple bulk downloads or unconstrained image requests.

To reduce very large SQLite rechecks, the live *conditional* main fallback now runs all 128 independently sealing shards inside a single verified SQLite session. This avoids 16 redundant full-database SHA256 and integrity scans while retaining 128 per-shard receipts. Source metadata processing still has no photo cap, and these engineering gains do not guarantee sufficient disk or within-job time. The separate fixed-run indexed fallback on [main](https://github.com/zuizui0223/fcp/blob/main/.github/workflows/fcp-pr200-fallback-indexed-source.yml)
is now configured to use this single-source, 8-shards-at-a-time path if the
still-running second bulk census fails its strict hash/phase gate. That
fallback still uses an ephemeral GitHub-hosted runner, so its checkpoints may
be lost if the *whole* job fails; a persistent machine is necessary for
durable **cross-job** resumption. The staging SQLite is excluded from the
portable source artifact only after all 128 shards have verified complete.

## Automated admission after source completion

Download **both** the successful source-status artifact and complete-source
tar.gz artifact from the **successful 42,111-species recovery run**, currently
[Actions #38699975166](https://github.com/zuizui0223/fcp/actions/runs/38699975166), onto durable storage, then unpack
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
bash scripts/acquisition/bootstrap_fcp_original_roi_v4.sh --full --sam-dir /persistent/models/fcp-original-efficient-sam
# Install frozen method dependencies in the intended production environment.
python -m pip install -r requirements-atlas-roi-v4.txt

python scripts/acquisition/execute_fcp_full_photo_work_queue.py \
  --firewall "$F" --output "$ROOT/partition_000" \
  --partition 0 --batch-size 64 --max-work-units 32 \
  --original-measure-script scripts/analysis/measure_random_photo_first_blind_partition.py \
  --original-detector data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt \
  --original-efficient-sam /persistent/models/fcp-original-efficient-sam
~~~

Repeat with identical source/method bytes and partition-local output, honoring one globally coordinated request budget. The selected `--partition` invocation now validates and plans **only that immutable source partition**, avoiding a scan of all 128 original-photo manifests twice per invocation. Its `n_pending_work_units=0` means **that partition alone** has finished; it is explicitly forbidden to set the full-study-complete flag. Only the independent all-128 collector/reassembler below can certify global coverage. The 32-work-unit invocation budget and 64-photo compute batch size never cap biological inclusion. **If a provider outage interrupts a unit**, its staging folder remains under `unsealed_inflight_batches`; the next execution deliberately rejects that folder. Inspect and explicitly quarantine the unsealed attempt before restarting the same frozen partition/batch, retaining its diagnostic evidence. Never move incomplete CSV/JSON into `sealed_terminal_batches`, relabel 429 photos as missing, or discard source IDs.

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
