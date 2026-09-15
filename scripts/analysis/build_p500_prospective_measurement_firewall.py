#!/usr/bin/env python3
"""Build the prospective P500 two-batch location-blind firewall without pixels."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.p500_prospective_execution_gate import (
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    FULL100_SHA256,
    METADATA_SHA256,
    validate_authorization,
    validate_candidate_metadata,
)
from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS,
    WORKER_FIELDS,
    compute_partition,
    measurement_id,
    semantic_shard,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CANDIDATE = ROOT / "data/frozen/polymorphism_h2_p500_candidate_metadata_v1.csv.gz"
DEFAULT_CANDIDATE_RESULT = ROOT / "results/polymorphism_h2_p500_candidate_metadata_20260913/result.json"
DEFAULT_FULL100 = ROOT / "results/polymorphism_h2_p500_candidate_metadata_20260913/full100_species.csv"
DEFAULT_AUTH = ROOT / "docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_AUTHORIZATION_20260915.json"
SALT = "FCP_H2_P500_PROSPECTIVE_20260915_V1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def measurement_batch(value: str, batches: int = 2) -> int:
    digest = hashlib.sha256(f"fcp-p500-batch\x1f{value}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % int(batches)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-csv", type=Path, default=DEFAULT_CANDIDATE)
    p.add_argument("--candidate-result", type=Path, default=DEFAULT_CANDIDATE_RESULT)
    p.add_argument("--full100", type=Path, default=DEFAULT_FULL100)
    p.add_argument("--authorization", type=Path, default=DEFAULT_AUTH)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    candidate_result = json.loads(args.candidate_result.read_text(encoding="utf-8"))
    authorization = json.loads(args.authorization.read_text(encoding="utf-8"))
    validate_candidate_metadata(candidate_result)
    validate_authorization(authorization)

    if sha256_file(args.candidate_csv) != METADATA_SHA256:
        raise RuntimeError("P500 candidate metadata byte hash mismatch")
    if sha256_file(args.full100) != FULL100_SHA256:
        raise RuntimeError("P500 full100 species byte hash mismatch")

    candidate = pd.read_csv(args.candidate_csv, low_memory=False)
    full100 = pd.read_csv(args.full100)
    if len(full100) != EXPECTED_SPECIES:
        raise RuntimeError("full100 species denominator mismatch")
    taxon_ids = set(pd.to_numeric(full100["inat_taxon_id"], errors="raise").astype(int))
    if len(taxon_ids) != EXPECTED_SPECIES:
        raise RuntimeError("full100 taxon IDs are not unique")

    required = {
        "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
        "photo_license", "latitude", "longitude",
    }
    missing = sorted(required - set(candidate.columns))
    if missing:
        raise RuntimeError(f"candidate metadata missing fields: {missing}")
    species_column = "query_species" if "query_species" in candidate.columns else "species"
    if species_column not in candidate.columns:
        raise RuntimeError("candidate metadata lacks species label")

    candidate["inat_taxon_id"] = pd.to_numeric(candidate["inat_taxon_id"], errors="raise").astype(int)
    measured = candidate.loc[candidate["inat_taxon_id"].isin(taxon_ids)].copy()
    if len(measured) != EXPECTED_ROWS:
        raise RuntimeError(f"P500 measured denominator mismatch: {len(measured)}")
    if measured["inat_taxon_id"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("P500 measured species denominator mismatch")
    counts = measured.groupby("inat_taxon_id", observed=True).size().astype(int)
    if not (counts == 100).all():
        raise RuntimeError("P500 per-species row target differs from 100")
    if measured["observation_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate observation IDs in P500 measurement denominator")
    if measured["photo_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate photo IDs in P500 measurement denominator")

    measured = measured.reset_index(drop=True)
    mids = [measurement_id(v, salt=SALT) for v in measured["photo_id"]]
    if len(set(mids)) != EXPECTED_ROWS:
        raise RuntimeError("P500 blinded measurement IDs are not unique")
    measured.insert(0, "measurement_id", mids)
    measured.insert(1, "measurement_batch", [measurement_batch(x) for x in mids])

    worker_all = pd.DataFrame({
        "measurement_id": measured["measurement_id"].astype(str),
        "image_filename": measured["measurement_id"].astype(str) + ".jpg",
        "photo_license": measured["photo_license"].astype(str),
    })
    acquisition_all = pd.DataFrame({
        "measurement_id": worker_all["measurement_id"],
        "image_filename": worker_all["image_filename"],
        "photo_url_large": measured["photo_url_large"].astype(str),
        "photo_license": worker_all["photo_license"],
    })
    if tuple(worker_all.columns) != WORKER_FIELDS or tuple(acquisition_all.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("location-blind interface fields drifted")

    metadata_join = measured.drop(columns=["photo_url_large"]).copy()
    if species_column != "species":
        metadata_join["species"] = metadata_join[species_column].astype(str)

    out = args.output_dir
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    metadata_join.to_csv(out / "sealed_keys/metadata_join_key.csv", index=False, lineterminator="\n")

    assignments = []
    batch_counts: dict[str, int] = {}
    partition_counts: dict[str, int] = {}
    for batch in (0, 1):
        keep = measured["measurement_batch"].astype(int).eq(batch).to_numpy()
        worker = worker_all.loc[keep].reset_index(drop=True)
        acquisition = acquisition_all.loc[keep].reset_index(drop=True)
        batch_dir = out / f"batch_{batch}"
        (batch_dir / "worker_packet").mkdir(parents=True, exist_ok=True)
        (batch_dir / "sealed_keys").mkdir(parents=True, exist_ok=True)
        worker.to_csv(batch_dir / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
        acquisition.to_csv(batch_dir / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
        batch_counts[str(batch)] = int(len(worker))

        local_rows = []
        for mid in worker["measurement_id"].astype(str):
            s = semantic_shard(mid, 32)
            p = compute_partition(mid, 4)
            key = f"b{batch}-s{s:02d}-p{p:02d}"
            partition_counts[key] = partition_counts.get(key, 0) + 1
            local_rows.append({
                "measurement_id": mid,
                "measurement_batch": batch,
                "semantic_shard": s,
                "compute_partition": p,
            })
        assignments.append(pd.DataFrame(local_rows))

    if sum(batch_counts.values()) != EXPECTED_ROWS:
        raise RuntimeError("P500 blind-batch denominator changed")
    if len(partition_counts) != 256:
        raise RuntimeError(f"not all 256 P500 terminal partitions are populated: {len(partition_counts)}")

    assignment = pd.concat(assignments, ignore_index=True).sort_values("measurement_id", kind="mergesort")
    assignment.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")

    manifest = {
        "schema": "p500_prospective_measurement_firewall_v1",
        "status": "p500_measurement_firewall_frozen_before_pixels",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "measurement_ids": int(worker_all["measurement_id"].nunique()),
        "measurement_id_salt": SALT,
        "blind_batches": 2,
        "rows_per_batch": batch_counts,
        "semantic_shards_per_batch": 32,
        "compute_partitions_per_semantic_shard": 4,
        "terminal_partitions": 256,
        "nonempty_terminal_partitions": len(partition_counts),
        "maximum_partition_rows": max(partition_counts.values()),
        "worker_fields": list(worker_all.columns),
        "acquisition_fields": list(acquisition_all.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
        "candidate_metadata_sha256": sha256_file(args.candidate_csv),
        "full100_species_sha256": sha256_file(args.full100),
        "authorization_sha256": sha256_file(args.authorization),
    }
    (out / "measurement_firewall_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
