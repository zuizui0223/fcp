#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    WORKER_FIELDS,
    ACQUISITION_FIELDS,
    measurement_id,
    semantic_shard,
    compute_partition,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv"
DEFAULT_MANIFEST = ROOT / "docs/supporting/random_photo_first_h9_fresh_metadata_manifest_v1.json"
DEFAULT_H9_EXECUTION = ROOT / "docs/supporting/random_photo_first_h9_measurement_execution_contract_v1.json"
DEFAULT_OLD_EXECUTION = ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--candidate-manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--h9-execution-contract", type=Path, default=DEFAULT_H9_EXECUTION)
    parser.add_argument("--old-execution-contract", type=Path, default=DEFAULT_OLD_EXECUTION)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    candidate = pd.read_csv(args.candidate_csv)
    manifest = json.loads(args.candidate_manifest.read_text())
    h9 = json.loads(args.h9_execution_contract.read_text())
    old = json.loads(args.old_execution_contract.read_text())

    gate = manifest["premeasurement_gate"]
    h9_gate = h9["pixel_gate"]
    if manifest["status"] != h9_gate["required_status"]:
        raise RuntimeError("H9 fresh metadata status mismatch")
    if gate["pass"] is not True or gate["decision"] != h9_gate["required_decision"]:
        raise RuntimeError("H9 measurement pixel gate did not pass")
    if len(candidate) != int(h9_gate["required_rows"]):
        raise RuntimeError("H9 measurement denominator row count mismatch")
    if sha256_file(args.candidate_csv) != str(h9_gate["required_observations_sha256"]):
        raise RuntimeError("H9 fresh metadata SHA-256 mismatch")
    if int(gate["full_species"]) != int(h9_gate["required_full_species"]):
        raise RuntimeError("H9 full-species count mismatch")
    if candidate["observation_id"].nunique() != len(candidate) or candidate["photo_id"].nunique() != len(candidate):
        raise RuntimeError("H9 frozen metadata IDs are not unique")
    if int(candidate["inat_taxon_id"].nunique()) != 38:
        raise RuntimeError("H9 candidate species denominator drifted")
    if not (candidate.groupby("inat_taxon_id").size() == 60).all():
        raise RuntimeError("H9 raw photo denominator is not exactly 60 per species")
    if manifest["candidate_image_pixels_opened"] is not False or manifest["colour_used_for_selection"] is not False:
        raise RuntimeError("H9 pre-pixel outcome firewall already opened")

    if old["protocol"] != "random-photo-first-measurement-execution-v1":
        raise RuntimeError("validated measurement execution protocol mismatch")
    salt = str(old["blinding"]["salt"])
    if salt != str(h9["blinding"]["measurement_id_salt"]):
        raise RuntimeError("H9 measurement-ID salt drifted")
    part = h9["partitioning"]
    if int(part["semantic_shards"]) != 32 or int(part["compute_partitions_per_semantic_shard"]) != 4:
        raise RuntimeError("H9 partition geometry drifted")

    metadata = candidate.reset_index(drop=True).copy()
    metadata.insert(0, "measurement_id", [measurement_id(v, salt=salt) for v in metadata["photo_id"]])
    if metadata["measurement_id"].nunique() != len(metadata):
        raise RuntimeError("H9 measurement IDs are not unique")
    worker = pd.DataFrame({
        "measurement_id": metadata["measurement_id"].astype(str),
        "image_filename": metadata["measurement_id"].astype(str) + ".jpg",
        "photo_license": metadata["photo_license"].astype(str),
    })
    acquisition = pd.DataFrame({
        "measurement_id": worker["measurement_id"],
        "image_filename": worker["image_filename"],
        "photo_url_large": metadata["photo_url_large"].astype(str),
        "photo_license": worker["photo_license"],
    })
    metadata_join = metadata.drop(columns=["photo_url_large"]).copy()
    if tuple(worker.columns) != WORKER_FIELDS or tuple(acquisition.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("H9 blind interface fields drifted")

    assignments = pd.DataFrame({
        "measurement_id": worker["measurement_id"].astype(str),
        "semantic_shard": [semantic_shard(x, 32) for x in worker["measurement_id"].astype(str)],
        "compute_partition": [compute_partition(x, 4) for x in worker["measurement_id"].astype(str)],
    })
    counts = assignments.groupby(["semantic_shard", "compute_partition"]).size()

    out = args.output_dir
    (out / "worker_packet").mkdir(parents=True, exist_ok=True)
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    worker.to_csv(out / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
    acquisition.to_csv(out / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
    metadata_join.to_csv(out / "sealed_keys/metadata_join_key.csv", index=False, lineterminator="\n")
    assignments.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")

    firewall = {
        "protocol": h9["protocol"],
        "status": "h9_measurement_firewall_frozen_before_pixels",
        "frozen_candidate_rows": int(len(candidate)),
        "frozen_species": int(candidate["inat_taxon_id"].nunique()),
        "measurement_ids": int(worker["measurement_id"].nunique()),
        "worker_fields": list(worker.columns),
        "acquisition_fields": list(acquisition.columns),
        "semantic_shards": 32,
        "compute_partitions_per_semantic_shard": 4,
        "total_compute_partitions": 128,
        "nonempty_compute_partitions": int(len(counts)),
        "maximum_rows_in_compute_partition": int(counts.max()) if len(counts) else 0,
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
        "fresh_metadata_sha256": sha256_file(args.candidate_csv),
    }
    (out / "measurement_firewall_manifest.json").write_text(json.dumps(firewall, indent=2) + "\n")
    print(json.dumps(firewall, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
