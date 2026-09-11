#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS,
    WORKER_FIELDS,
    compute_partition,
    semantic_shard,
)

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_42111_breadth_measurement_contract_v1.json"
STEP8D_RESULT = ROOT / "results/rgfca_42111_anchor_resolution_full_step8d_20260911/result.json"
STEP8D_TABLE = ROOT / "results/rgfca_42111_anchor_resolution_full_step8d_20260911/anchor_metadata_42111.csv.gz"
STEP8E_GATE = ROOT / "results/rgfca_42111_breadth_prepixel_step8e_20260911/result.json"
OUT_DEFAULT = ROOT / "results/rgfca_42111_breadth_measurement_firewall_step8f_20260911"

BIOLOGICAL = {"white", "yellow_orange", "red_pink", "blue_purple"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def measurement_id(photo_id: int, salt: str) -> str:
    payload = f"{salt}\x1fphoto\x1f{int(photo_id)}".encode("utf-8")
    return "FCPB-" + hashlib.sha256(payload).hexdigest().upper()[:24]


def unresolved_id(taxon_id: int, salt: str) -> str:
    payload = f"{salt}\x1funresolved\x1f{int(taxon_id)}".encode("utf-8")
    return "FCPB-U-" + hashlib.sha256(payload).hexdigest().upper()[:22]


def main(output_dir: Path = OUT_DEFAULT) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    d_result = json.loads(STEP8D_RESULT.read_text(encoding="utf-8"))
    gate = json.loads(STEP8E_GATE.read_text(encoding="utf-8"))

    if contract.get("status") != "authorize_exactly_one_breadth_measurement_after_transport_gate_before_pixels":
        raise RuntimeError("breadth measurement contract is not in the frozen authorization state")
    if d_result.get("status") != "complete_metadata_only_full_anchor_resolution":
        raise RuntimeError("Step 8D anchor resolution incomplete")
    if gate.get("status") != contract["transport_gate"]["required_status"] or gate.get("transport_evaluable") is not True:
        raise RuntimeError("Step 8E transport gate did not pass")
    if int(gate.get("species_denominator", -1)) != int(contract["species_universe"]):
        raise RuntimeError("species denominator drift")
    if int(gate.get("resolved_exact_anchor_species", -1)) != int(contract["resolved_anchor_species"]):
        raise RuntimeError("resolved-anchor denominator drift")
    if float(gate.get("resolved_fraction", 0.0)) < float(contract["transport_gate"]["minimum_resolved_fraction"]):
        raise RuntimeError("transport resolution fraction below frozen floor")
    if gate.get("image_pixels_opened") is not False or gate.get("flower_colour_used") is not False:
        raise RuntimeError("pre-pixel gate already opened forbidden outcomes")
    if sha256_file(STEP8D_TABLE) != contract["lineage"]["required_step8d_resolved_table_sha256"]:
        raise RuntimeError("Step 8D resolved table hash drift")
    if gate.get("lineage", {}).get("denominator_sha256") != contract["lineage"]["required_step8e_denominator_sha256"]:
        raise RuntimeError("Step 8E denominator hash drift")
    if gate.get("lineage", {}).get("source_manifest_sha256") != contract["lineage"]["required_step8e_source_manifest_sha256"]:
        raise RuntimeError("Step 8E source-manifest hash drift")

    anchors = pd.read_csv(STEP8D_TABLE)
    if len(anchors) != 42111 or anchors["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("Step 8D table is not the exact 42,111-species universe")
    resolved = anchors.loc[anchors["status"].astype(str).eq("resolved")].copy()
    unresolved = anchors.loc[~anchors["status"].astype(str).eq("resolved")].copy()
    if len(resolved) != 42110 or len(unresolved) != 1:
        raise RuntimeError("resolved/unresolved counts differ from frozen contract")
    if resolved["photo_url_large"].fillna("").astype(str).eq("").any():
        raise RuntimeError("resolved anchor lacks source URL")
    if resolved["photo_id"].nunique() != len(resolved) or resolved["observation_id"].nunique() != len(resolved):
        raise RuntimeError("resolved species anchors are not unique identities")

    salt = str(contract["blinding"]["measurement_id_salt"])
    resolved["measurement_id"] = [measurement_id(int(v), salt) for v in resolved["photo_id"]]
    if resolved["measurement_id"].nunique() != len(resolved):
        raise RuntimeError("breadth measurement IDs are not unique")
    resolved["image_filename"] = resolved["measurement_id"].astype(str) + ".jpg"

    worker = resolved[["measurement_id", "image_filename", "photo_license"]].copy()
    acquisition = resolved[["measurement_id", "image_filename", "photo_url_large", "photo_license"]].copy()
    if tuple(worker.columns) != WORKER_FIELDS or tuple(acquisition.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("validated blind worker/acquisition interface drifted")

    metadata_cols = [
        "measurement_id", "species", "inat_taxon_id", "breadth_rank", "observation_id", "photo_id",
        "anchor_source", "opened_existing_species", "status", "photo_license", "attribution"
    ]
    missing = [c for c in metadata_cols if c not in resolved.columns]
    if missing:
        raise RuntimeError(f"resolved table missing metadata columns: {missing}")
    metadata_join = resolved[metadata_cols].copy()
    if any(c in metadata_join.columns for c in ["photo_url_large", "latitude", "longitude"]):
        raise RuntimeError("metadata join unexpectedly contains source URL or coordinates")

    u = unresolved.iloc[0]
    unresolved_terminal = pd.DataFrame([{
        "measurement_id": unresolved_id(int(u["inat_taxon_id"]), salt),
        "species": str(u["species"]),
        "inat_taxon_id": int(u["inat_taxon_id"]),
        "breadth_rank": int(u["breadth_rank"]),
        "observation_id": int(u["observation_id"]),
        "photo_id": int(u["photo_id"]),
        "anchor_source": str(u["anchor_source"]),
        "opened_existing_species": bool(u["opened_existing_species"]),
        "anchor_resolution_status": str(u["status"]),
        "morph": "mixed_uncertain",
        "measurement_status": "anchor_resolution_failed",
        "roi_status": "not_opened_anchor_resolution_failed",
        "failure_reasons": str(u["status"]),
        "flower_effective_pixels": 0,
    }])

    n_sem = int(contract["partitioning"]["semantic_shards"])
    n_comp = int(contract["partitioning"]["compute_partitions_per_semantic_shard"])
    assignments = pd.DataFrame({
        "measurement_id": resolved["measurement_id"].astype(str),
        "semantic_shard": [semantic_shard(v, n_sem) for v in resolved["measurement_id"].astype(str)],
        "compute_partition": [compute_partition(v, n_comp) for v in resolved["measurement_id"].astype(str)],
    })
    counts = assignments.groupby(["semantic_shard", "compute_partition"], observed=True).size()
    if len(counts) != n_sem * n_comp:
        raise RuntimeError("one or more frozen measurement partitions are empty unexpectedly")

    (output_dir / "worker_packet").mkdir(parents=True, exist_ok=True)
    (output_dir / "sealed_keys").mkdir(parents=True, exist_ok=True)
    worker.to_csv(output_dir / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
    acquisition.to_csv(output_dir / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
    metadata_join.to_csv(output_dir / "sealed_keys/metadata_join_key.csv", index=False, lineterminator="\n")
    unresolved_terminal.to_csv(output_dir / "sealed_keys/unresolved_terminal_species.csv", index=False, lineterminator="\n")
    assignments.to_csv(output_dir / "partition_assignments.csv", index=False, lineterminator="\n")

    manifest = {
        "protocol": contract["protocol"],
        "status": "breadth_measurement_firewall_frozen_before_pixels",
        "species_denominator": 42111,
        "resolved_measurement_rows": 42110,
        "preexisting_unresolved_terminal_rows": 1,
        "measurement_ids": int(worker["measurement_id"].nunique()),
        "semantic_shards": n_sem,
        "compute_partitions_per_semantic_shard": n_comp,
        "terminal_partitions": n_sem * n_comp,
        "minimum_partition_rows": int(counts.min()),
        "maximum_partition_rows": int(counts.max()),
        "worker_fields": list(worker.columns),
        "acquisition_fields": list(acquisition.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "coordinate_or_species_colour_join_opened": False,
        "lineage": {
            "outer_contract_sha256": sha256_file(CONTRACT),
            "step8d_table_sha256": sha256_file(STEP8D_TABLE),
            "step8e_gate_sha256": sha256_file(STEP8E_GATE),
            "worker_sha256": sha256_file(output_dir / "worker_packet/measurement_manifest.csv"),
            "acquisition_key_sha256": sha256_file(output_dir / "sealed_keys/acquisition_key.csv"),
            "metadata_join_sha256": sha256_file(output_dir / "sealed_keys/metadata_join_key.csv"),
            "unresolved_terminal_sha256": sha256_file(output_dir / "sealed_keys/unresolved_terminal_species.csv"),
        }
    }
    (output_dir / "firewall_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
