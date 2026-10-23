#!/usr/bin/env python3
"""Collect disjoint partition-local FCP original-photo results, fail closed.

Run AFTER workers are quiescent. Workers must use *separate* output directories
(e.g. partition_000/...); they must not write to one shared progress.json.
Never opens image pixels, geography, or historical flower-colour measurements.
Only a complete 42,111-source firewall permits production collection.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import validate_method_attestation
from scripts.acquisition.plan_fcp_full_measurements import plan_jobs


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _copy_immutable(source: Path, target: Path) -> None:
    """Resume an interrupted pair copy only if its original bytes still agree."""
    digest = sha256_file(source)
    if target.is_file():
        if sha256_file(target) != digest:
            raise ValueError(f"Previously collected FCP terminal changed: {target.name}")
        return
    if target.exists():
        raise ValueError(f"Terminal destination is not a regular file: {target.name}")
    temporary = target.with_name("." + target.name + ".copying")
    if temporary.exists():
        raise ValueError(f"Unverified interrupted copy needs inspection: {temporary.name}")
    shutil.copyfile(source, temporary)
    if sha256_file(temporary) != digest:
        raise ValueError(f"Source FCP terminal changed during copying: {source.name}")
    os.replace(temporary, target)


def collect(firewall: Path, partition_root: Path, output: Path, *,
            batch_size: int = 64, _test_only_small_source: bool = False) -> dict:
    freeze = _read(firewall / "fresh_study_freeze.json")
    if (freeze.get("status") != "FRESH_BLIND_METADATA_FROZEN"
            or freeze.get("historical_measured_CSV_read") is not False
            or freeze.get("new_measurements_completed") != 0):
        raise ValueError("Complete, untouched fresh-source photo firewall required")
    if not _test_only_small_source and (
            freeze.get("source_species") != 42111 or freeze.get("worker_packets") != 128):
        raise ValueError("Cannot collect production results without all 42111 taxa and 128 partitions")
    if not 1 <= batch_size <= 1024:
        raise ValueError("Work-unit size outside frozen execution bounds")
    expected_names = {f"partition_{p:03}" for p in range(int(freeze["worker_packets"]))}
    if partition_root.exists():
        unknown = {p.name for p in partition_root.iterdir()} - expected_names
        if unknown:
            raise ValueError(f"Unknown measurement partitions: {sorted(unknown)[:3]}")
    output.mkdir(parents=True, exist_ok=True)
    dest = output / "sealed_terminal_batches"
    dest.mkdir(exist_ok=True)
    progress_file = output / "partition_collection_progress.json"
    previous = _read(progress_file) if progress_file.exists() else None
    if previous is not None and previous.get("source_firewall_sha256") != freeze["blind_firewall_sha256"]:
        raise ValueError("Cannot mix different frozen source-photo studies in collection")
    found = 0
    contributor_partitions = 0
    for part in range(int(freeze["worker_packets"])):
        folder = partition_root / f"partition_{part:03}" / "sealed_terminal_batches"
        if not folder.exists():
            continue
        if not folder.is_dir():
            raise ValueError(f"Partition {part}: sealed output is not a directory")
        files = list(folder.iterdir())
        pattern = re.compile(rf"terminal_p{part:03}_b[0-9]{{7}}\.(?:csv|json)$")
        other = [p.name for p in files if p.name != "frozen_method_attestation.json"
                 and not pattern.fullmatch(p.name)]
        if other:
            raise ValueError(f"Unexpected or wrong-partition ROI-v4 terminal: {other[:3]}")
        csvs = sorted(p for p in files if p.suffix == ".csv")
        jsons = {p.stem for p in files if p.name.startswith("terminal_") and p.suffix == ".json"}
        if {p.stem for p in csvs} != jsons:
            raise ValueError(f"Partition {part}: unsealed CSV/receipt pair")
        if not csvs:
            continue
        method_path = folder / "frozen_method_attestation.json"
        if not method_path.is_file():
            raise ValueError(f"Partition {part}: missing frozen ROI-v4 method attestation")
        method = validate_method_attestation(method_path)
        _copy_immutable(method_path, dest / "frozen_method_attestation.json")
        contributor_partitions += 1
        for src in csvs:
            receipt_file = src.with_suffix(".json")
            receipt = _read(receipt_file)
            batch = int(src.stem.split("_b")[1])
            if (receipt.get("status") != "COMPLETE_WORK_UNIT"
                    or receipt.get("partition") != part or receipt.get("batch") != batch
                    or receipt.get("terminal_csv_sha256") != sha256_file(src)
                    or receipt.get("method") != method
                    or receipt.get("reused_historical_measurement_csv") is not False):
                raise ValueError(f"Partition {part}: invalid source-bound photo receipt {src.name}")
            _copy_immutable(src, dest / src.name)
            _copy_immutable(receipt_file, dest / receipt_file.name)
            found += 1
    with tempfile.TemporaryDirectory(prefix="fcp_partition_collect_audit_") as work:
        audit = plan_jobs(firewall, dest, Path(work) / "plan",
                          batch_size=batch_size,
                          _test_only_small_source=_test_only_small_source)
    if previous is not None and audit["n_completed_work_units"] < previous["n_completed_work_units"]:
        raise ValueError("Production collection cannot regress previously sealed ROI-v4 units")
    complete = audit["n_pending_work_units"] == 0
    result = {
        "schema": "fcp_partition_local_photo_terminal_collection_v1",
        "status": ("ALL_ORIGINAL_PHOTO_TERMINALS_SEALED_READY_FOR_REASSEMBLY_NOT_H2"
                   if complete else "PARTIAL_ORIGINAL_PHOTO_TERMINALS_NOT_SCIENTIFIC_RESULT"),
        "source_species": audit["source_species"],
        "source_firewall_sha256": audit["source_firewall_sha256"],
        "source_photo_ids": audit["source_photo_ids"],
        "n_work_units": audit["work_units"],
        "n_completed_work_units": audit["n_completed_work_units"],
        "n_pending_work_units": audit["n_pending_work_units"],
        "partitions_with_sealed_work": contributor_partitions,
        "n_source_batches_seen": found,
        "complete_source_photo_terminal_coverage": complete,
        "fresh_reassembly_completed": False,
        "W_H2_environment_computed": False,
        "old_flower_colour_csv_reused": False,
        "synthetic_test_only": _test_only_small_source,
    }
    temp = progress_file.with_suffix(".json.tmp")
    if temp.exists():
        raise ValueError("Previous interrupted FCP collection progress write needs review")
    temp.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temp, progress_file)
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--firewall", type=Path, required=True)
    p.add_argument("--partition-output-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--batch-size", type=int, default=64)
    args = p.parse_args()
    print(json.dumps(collect(args.firewall, args.partition_output_root,
                             args.output, batch_size=args.batch_size), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
