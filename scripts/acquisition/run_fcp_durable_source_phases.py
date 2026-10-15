#!/usr/bin/env python3
"""Run immutable 42,111-taxon FCP photo-ID source phases ONLY on durable storage.

The original 93,209,351 count-only evidence is NOT the photo-ID universe.
This runner performs the missing 4-table -> SQLite -> 128 shards -> firewall
freeze sequence, then independently reconciles every source taxon's photo count.
It never downloads original photograph pixels or reads old colour outcomes.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
ARCHIVE_SHA = "bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
COUNTS_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
AUDIT_SHA = "ed62aee0f9730ad17347d25180575af1cd101e0fc35c581e470c4a17866ebf56"
SOURCE_PHOTOS = 93_209_351
SPECIES = 42_111
SHARDS = 128
DATE = "2026-09-27"
SOURCE_STAGED = "SOURCE_SQLITE_STAGED_FOR_INDEPENDENT_SHARDS_NOT_FULL_CENSUS"
FULL_SOURCE = "COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE"
ADMITTED = "REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS"
SOURCE_ARCHIVE = "inaturalist-open-data-20260927.tar.gz"
MIN_BYTES = {"preflight": 200_000_000_000, "stage": 200_000_000_000,
             "emit": 60_000_000_000, "finalize": 60_000_000_000}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read(path: Path) -> dict:
    x = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(x, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return x


def validate_root(root: Path, runner_temp: str | None = None) -> Path:
    if not root.is_absolute():
        raise ValueError("FCP_DURABLE_ROOT must be an absolute path")
    if root.is_symlink():
        raise ValueError("FCP_DURABLE_ROOT must not be a symbolic link")
    root = root.resolve()
    forbidden = [Path("/"), Path("/tmp"), Path("/var/tmp"), Path("/dev/shm"),
                 Path("/home/runner/work"), Path("/home/runner/_work")]
    if runner_temp:
        forbidden.append(Path(runner_temp).resolve())
    if any(root == x or x in root.parents for x in forbidden):
        raise ValueError("Ephemeral/unsafe FCP source storage; use a persistent mounted volume")
    if not root.is_dir() or root.is_symlink():
        raise ValueError("Persistent FCP root directory must already exist and not be a symlink")
    # This explicit sentinel must be provisioned by the runner operator.
    # Its presence does not prove hardware durability; the operator owns that check.
    marker = root / "FCP_DURABLE_STORAGE_APPROVED"
    if not marker.is_file() or marker.is_symlink() or marker.read_text(encoding="utf-8").strip() != (
            "I confirm this is persistent FCP source storage, not hosted-runner scratch."):
        raise ValueError("No explicit persistent-storage operator sentinel")
    return root


def validate_shards(start: int, stop: int) -> None:
    if not 0 <= start < stop <= SHARDS:
        raise ValueError("Shard range must be half-open within [0,128)")
    if stop - start > 4:
        raise ValueError("One invocation may emit at most four original source shards; no biological cap")


def validate_reference(root: Path) -> dict:
    frame = root / "inputs" / "42111_species.csv"
    counts = root / "references" / "all_42111_species_photo_opportunity_counts.csv"
    audit = root / "references" / "all_42111_photo_opportunity_audit.json"
    expected = ((frame, FRAME_SHA), (counts, COUNTS_SHA), (audit, AUDIT_SHA))
    for path, checksum in expected:
        if not path.is_file() or path.is_symlink() or sha256(path) != checksum:
            raise ValueError(f"Missing/changed fixed original 42111 source reference: {path.name}")
    result = _read(audit)
    if (result.get("status") != "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
        or result.get("snapshot_date") != DATE
        or result.get("source_frame_species") != SPECIES
        or result.get("source_frame_sha256") != FRAME_SHA
        or result.get("source_archive_sha256") != ARCHIVE_SHA
        or result.get("count_ledger_sha256") != COUNTS_SHA
        or result.get("all_source_photo_opportunity_ids") != SOURCE_PHOTOS
        or result.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or result.get("all_source_photo_id_manifests_frozen") is not False
        or result.get("all_original_image_pixels_measured") is not False):
        raise ValueError("Source-reference audit failed original 42111 species and no-colour gates")
    return {"species": SPECIES, "source_photo_opportunities": SOURCE_PHOTOS,
            "original_date": DATE, "legacy_colour_rows_used": 0}


def run_phase(root: Path, phase: str, start: int = 0, stop: int = 4,
              *, repo: Path | None = None) -> dict:
    if phase not in MIN_BYTES:
        raise ValueError("Unknown FCP durable phase")
    if phase == "emit":
        validate_shards(start, stop)
    root = validate_root(root, os.environ.get("RUNNER_TEMP"))
    repo = (repo or Path(__file__).resolve().parents[2]).resolve()
    if not (repo / "scripts/acquisition/execute_fcp_42111_census.py").is_file():
        raise ValueError("Run from the exact reviewed FCP repository checkout")
    root_status = root / "status"
    root_status.mkdir(exist_ok=True)
    lock = root / ".fcp_source_phase.lock"
    with lock.open("a+") as lock_handle:
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another process is using this same persistent FCP source") from exc

        summary = validate_reference(root)
        work = root / "source_20260927"
        frame = root / "inputs/42111_species.csv"
        archive = root / "snapshots" / SOURCE_ARCHIVE
        staged = work / "complete_42111_metadata_source"
        receipt = work / "full_source_execution_receipt.json"
        if phase in {"preflight", "stage"}:
            if (not archive.is_file() or archive.is_symlink()
                or archive.stat().st_size != 35_966_549_880
                or sha256(archive) != ARCHIVE_SHA):
                raise ValueError("Exact original 36GB dated archive absent or SHA256 mismatch")
        if phase == "stage" and work.exists() and any(work.iterdir()):
            raise ValueError("Existing staged-source work requires inspection; never overwrite a partial 93M source")
        if phase == "emit" and (not receipt.is_file() or _read(receipt).get("status") != SOURCE_STAGED):
            raise ValueError("No completely staged, source SHA-verified original SQLite photo-ID universe")
        if phase == "finalize" and (not receipt.is_file() or _read(receipt).get("status") != SOURCE_STAGED):
            raise ValueError("All 128 source shards must be emitted from immutable SQLite before finalization")
        free = shutil.disk_usage(root).free
        if free < MIN_BYTES[phase]:
            raise ValueError(f"Persistent disk free={free} bytes; phase {phase} requires >= {MIN_BYTES[phase]}")

        commands = {
            "stage": [
                sys.executable, "scripts/acquisition/execute_fcp_42111_census.py",
                "--frame", str(frame), "--archive", str(archive), "--work", str(work),
                "--expected-species", str(SPECIES), "--shards", str(SHARDS),
                "--worker-partitions", str(SHARDS), "--stage-only",
                "--release-verified-temporary-tables",
            ],
            "emit": [
                sys.executable, "scripts/acquisition/emit_fcp_staged_source_shards.py",
                "--frame", str(frame), "--staged-metadata", str(staged),
                "--expected-species", str(SPECIES), "--shards", str(SHARDS),
                "--shard-start", str(start), "--shard-stop", str(stop),
            ],
            "finalize": [
                sys.executable, "scripts/acquisition/finalize_fcp_staged_source.py",
                "--frame", str(frame), "--work", str(work), "--expected-species", str(SPECIES),
                "--shards", str(SHARDS), "--worker-partitions", str(SHARDS),
                "--skip-package",
            ],
        }
        if phase != "preflight":
            subprocess.run(commands[phase], check=True, cwd=repo)
            if phase == "stage":
                doc = _read(receipt)
                if (doc.get("status") != SOURCE_STAGED or
                    doc.get("staged_original_source_photo_ids") != SOURCE_PHOTOS or
                    doc.get("complete_metadata_census") is not False):
                    raise ValueError("Original SQLite stage did not preserve exact 93,209,351 photo IDs")
            if phase == "finalize":
                admitted = work / "independent_original_all_ID_source_admission.json"
                subprocess.run([
                    sys.executable, "scripts/acquisition/verify_fcp_complete_photo_ids_vs_original_counts.py",
                    "--baseline-counts", str(root / "references/all_42111_species_photo_opportunity_counts.csv"),
                    "--independent-audit", str(root / "references/all_42111_photo_opportunity_audit.json"),
                    "--complete-source-counts", str(work / "all_42111_species_source_photo_counts.csv"),
                    "--full-source-receipt", str(receipt),
                    "--firewall-freeze", str(work / "frozen_location_blind_all_photo_firewall/fresh_study_freeze.json"),
                    "--output", str(admitted),
                ], check=True, cwd=repo)
                doc = _read(admitted)
                if (doc.get("status") != ADMITTED or
                    doc.get("distinct_original_photo_ids") != SOURCE_PHOTOS or
                    doc.get("full_study_original_pixels_classified") is not False):
                    raise ValueError("Complete photo-ID admission fails original provider census")
                summary["full_source_photo_ID_universe_frozen"] = True

        summary.update(phase=phase, start=start if phase == "emit" else None,
                       stop=stop if phase == "emit" else None,
                       original_photo_pixels_measured=False, flower_colour_W_estimated=False,
                       status=("SOURCE_PHOTO_IDS_FROZEN_RECONCILED_NO_COLOURS"
                               if phase == "finalize" else "ORIGINAL_SOURCE_STAGE_ONLY_NO_COLOURS"),
                       completed_utc=datetime.now(timezone.utc).isoformat())
        target = root_status / f"phase_{phase}_{start:03}_{stop:03}.json"
        pending = target.with_suffix(".json.tmp")
        pending.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(pending, target)
        return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--phase", required=True, choices=list(MIN_BYTES))
    p.add_argument("--shard-start", type=int, default=0)
    p.add_argument("--shard-stop", type=int, default=4)
    args = p.parse_args(argv)
    print("FCP_DURABLE_PHASE", json.dumps(run_phase(args.root, args.phase,
          args.shard_start, args.shard_stop), sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
