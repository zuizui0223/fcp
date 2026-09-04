#!/usr/bin/env python3
"""Build the global two-batch measurement firewall without opening image pixels."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    WORKER_FIELDS,
    ACQUISITION_FIELDS,
    semantic_shard,
    compute_partition,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
DEFAULT_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"
DEFAULT_EXECUTION = ROOT / "docs/supporting/global_monte_carlo_measurement_execution_contract_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _blind(salt: str, label: str, value: object, *, length: int = 24) -> str:
    payload = f"{salt}\x1f{label}\x1f{value}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()[:length]


def measurement_id(photo_id: object, *, salt: str) -> str:
    return "FCPG-" + _blind(salt, "photo", photo_id)


def measurement_batch(measurement: str, batches: int = 2) -> int:
    if int(batches) < 1:
        raise ValueError("batches must be positive")
    digest = hashlib.sha256(
        f"fcp-global-measurement-batch\x1f{measurement}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big") % int(batches)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATE)
    p.add_argument("--candidate-manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    candidate = pd.read_csv(args.candidate_csv)
    manifest = json.loads(args.candidate_manifest.read_text(encoding="utf-8"))
    execution = json.loads(args.execution_contract.read_text(encoding="utf-8"))

    if execution.get("status") != "frozen_before_capacity_outcome_before_candidate_acquisition_outcome_and_before_global_candidate_pixels":
        raise RuntimeError("global measurement execution contract is not frozen pre-pixel")
    if manifest.get("status") != execution["pixel_opening"]["requires_candidate_status"]:
        raise RuntimeError("global candidate acquisition status does not authorize firewall construction")
    gate = manifest.get("premeasurement_gate", {})
    if gate.get("pass") is not True:
        raise RuntimeError("global candidate premeasurement gate did not pass")
    if manifest.get("candidate_image_pixels_opened") is not False or manifest.get("flower_colour_used") is not False:
        raise RuntimeError("global candidate stage already opened forbidden outcomes")
    if manifest.get("measurement_authorized") is not False:
        raise RuntimeError("candidate manifest unexpectedly self-authorized measurement")
    expected_sha = str(manifest.get("lineage", {}).get("candidate_photos_sha256") or "")
    if not expected_sha or sha256_file(args.candidate_csv) != expected_sha:
        raise RuntimeError("candidate photo pool SHA-256 differs from frozen manifest")

    target = int(manifest["capacity_selected_raw_photo_target"])
    full_species = int(manifest["full_target_species"])
    expected_rows = full_species * target
    if len(candidate) != expected_rows or int(manifest["candidate_rows"]) != expected_rows:
        raise RuntimeError("global candidate denominator does not equal full_species * target")
    required = {"species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large", "photo_license", "latitude", "longitude"}
    missing = sorted(required - set(candidate.columns))
    if missing:
        raise RuntimeError(f"candidate table missing firewall fields: {missing}")
    if candidate["observation_id"].nunique() != len(candidate) or candidate["photo_id"].nunique() != len(candidate):
        raise RuntimeError("global candidate observation/photo IDs are not unique")
    counts = candidate.groupby("inat_taxon_id", observed=True).size()
    if len(counts) != full_species or not (counts.astype(int) == target).all():
        raise RuntimeError("global candidate per-species denominator drifted")

    partition = execution["partitioning"]
    batches = int(partition["blind_batches"])
    validated = partition["validated_partitions_per_batch"]
    semantic_shards = int(validated["semantic_shards"])
    compute_partitions = int(validated["compute_partitions_per_semantic_shard"])
    if batches != 2 or semantic_shards != 32 or compute_partitions != 4:
        raise RuntimeError("global measurement batch/validated partition geometry drifted")
    if batches * semantic_shards * compute_partitions != int(partition["total_terminal_partitions"]):
        raise RuntimeError("global measurement total terminal partition count drifted")

    metadata = candidate.reset_index(drop=True).copy()
    salt = str(execution["blinding"]["measurement_id_salt"])
    metadata.insert(0, "measurement_id", [measurement_id(v, salt=salt) for v in metadata["photo_id"]])
    metadata.insert(1, "measurement_batch", [measurement_batch(v, batches) for v in metadata["measurement_id"].astype(str)])
    if metadata["measurement_id"].nunique() != len(metadata):
        raise RuntimeError("global blinded measurement IDs are not unique")

    worker_all = pd.DataFrame({
        "measurement_id": metadata["measurement_id"].astype(str),
        "image_filename": metadata["measurement_id"].astype(str) + ".jpg",
        "photo_license": metadata["photo_license"].astype(str),
    })
    acquisition_all = pd.DataFrame({
        "measurement_id": worker_all["measurement_id"],
        "image_filename": worker_all["image_filename"],
        "photo_url_large": metadata["photo_url_large"].astype(str),
        "photo_license": worker_all["photo_license"],
    })
    if tuple(worker_all.columns) != WORKER_FIELDS or tuple(acquisition_all.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("validated blind interface fields drifted")
    metadata_join = metadata.drop(columns=["photo_url_large"]).copy()

    out = args.output_dir
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    metadata_join.to_csv(out / "sealed_keys/metadata_join_key.csv", index=False, lineterminator="\n")
    assignments = []
    batch_counts = {}
    max_partition_rows = 0
    nonempty = 0
    for batch in range(batches):
        keep = metadata["measurement_batch"].astype(int).eq(batch).to_numpy()
        worker = worker_all.loc[keep].reset_index(drop=True)
        acquisition = acquisition_all.loc[keep].reset_index(drop=True)
        batch_dir = out / f"batch_{batch}"
        (batch_dir / "worker_packet").mkdir(parents=True, exist_ok=True)
        (batch_dir / "sealed_keys").mkdir(parents=True, exist_ok=True)
        worker.to_csv(batch_dir / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
        acquisition.to_csv(batch_dir / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
        batch_counts[str(batch)] = int(len(worker))
        local = pd.DataFrame({
            "measurement_id": worker["measurement_id"].astype(str),
            "measurement_batch": batch,
            "semantic_shard": [semantic_shard(x, semantic_shards) for x in worker["measurement_id"].astype(str)],
            "compute_partition": [compute_partition(x, compute_partitions) for x in worker["measurement_id"].astype(str)],
        })
        assignments.append(local)
        c = local.groupby(["semantic_shard", "compute_partition"], observed=True).size()
        nonempty += int(len(c))
        if len(c):
            max_partition_rows = max(max_partition_rows, int(c.max()))
    assignment = pd.concat(assignments, ignore_index=True).sort_values("measurement_id", kind="mergesort")
    assignment.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")
    if sum(batch_counts.values()) != len(candidate):
        raise RuntimeError("blind batch split changed the candidate denominator")

    firewall = {
        "protocol": execution["protocol"],
        "status": "global_measurement_firewall_frozen_before_pixels",
        "frozen_candidate_rows": int(len(candidate)),
        "frozen_species": full_species,
        "selected_raw_photo_target": target,
        "measurement_ids": int(worker_all["measurement_id"].nunique()),
        "blind_batches": batches,
        "rows_per_blind_batch": batch_counts,
        "validated_semantic_shards_per_batch": semantic_shards,
        "validated_compute_partitions_per_semantic_shard": compute_partitions,
        "total_terminal_partitions": int(batches * semantic_shards * compute_partitions),
        "nonempty_terminal_partitions": nonempty,
        "maximum_rows_in_terminal_partition": max_partition_rows,
        "worker_fields": list(worker_all.columns),
        "acquisition_fields": list(acquisition_all.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
        "candidate_photos_sha256": sha256_file(args.candidate_csv),
        "execution_contract_sha256": sha256_file(args.execution_contract),
    }
    (out / "measurement_firewall_manifest.json").write_text(json.dumps(firewall, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(firewall, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
