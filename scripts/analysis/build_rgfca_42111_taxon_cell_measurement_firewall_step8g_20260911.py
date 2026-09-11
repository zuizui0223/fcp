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
CONTRACT = ROOT / "docs/supporting/rgfca_42111_taxon_cell_measurement_contract_v1.json"
MAP_RESULT = ROOT / "results/rgfca_42111_taxon_cell_resolution_step8e3_20260911/result.json"
MAP_TABLE = ROOT / "results/rgfca_42111_taxon_cell_resolution_step8e3_20260911/taxon_cell_anchor_metadata_85337.csv.gz"
BREADTH_RESULT = ROOT / "results/rgfca_42111_breadth_measurement_step8f_20260911/result.json"
BREADTH_TABLE = ROOT / "results/rgfca_42111_breadth_measurement_step8f_20260911/rgfca_42111_species_breadth_measured.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxon_cell_measurement_firewall_step8g_20260911"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def measurement_id(photo_id: int, salt: str) -> str:
    payload = f"{salt}\x1fphoto\x1f{int(photo_id)}".encode("utf-8")
    return "FCPG-" + hashlib.sha256(payload).hexdigest().upper()[:24]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    mr = json.loads(MAP_RESULT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_before_any_taxon_cell_specific_image_pixel":
        raise RuntimeError("taxon-cell measurement contract status drift")
    if mr.get("status") != "complete_metadata_only_taxon_cell_resolution":
        raise RuntimeError("taxon-cell transport is incomplete")
    if mr.get("image_pixels_opened") is not False or mr.get("flower_colour_used") is not False:
        raise RuntimeError("taxon-cell transport opened forbidden outcomes")
    if int(mr.get("taxon_cell_pairs", -1)) != 85337 or int(mr.get("species", -1)) != 42111 or int(mr.get("occupied_cells", -1)) != 128:
        raise RuntimeError("taxon-cell transport denominator drift")
    if sha256_file(MAP_TABLE) != contract["transport"]["required_source_table_sha256"]:
        raise RuntimeError("taxon-cell transport table hash drift")

    # Deliberately fail closed until the one-shot species-breadth measurement is complete.
    if not BREADTH_RESULT.exists() or not BREADTH_TABLE.exists():
        raise RuntimeError("breadth prerequisite absent; taxon-cell-specific pixels remain closed")
    br = json.loads(BREADTH_RESULT.read_text(encoding="utf-8"))
    if br.get("status") != contract["breadth_reuse_prerequisite"]["required_breadth_status"]:
        raise RuntimeError("breadth prerequisite is not complete")
    if int(br.get("final_rows", -1)) != 42111 or int(br.get("terminal_partition_receipts", -1)) != 128:
        raise RuntimeError("breadth terminal denominator incomplete")
    if br.get("species_metadata_join_opened_only_after_complete_measurement") is not True:
        raise RuntimeError("breadth result did not satisfy complete-join firewall")

    m = pd.read_csv(MAP_TABLE)
    b = pd.read_csv(BREADTH_TABLE)
    if len(m) != 85337 or m[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != 85337:
        raise RuntimeError("map table is not exact one row per frozen taxon-cell")
    if len(b) != 42111 or b["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("breadth table is not exact 42,111 species")

    reuse = m.loc[m["resolution_status"].astype(str).eq("reused_species_anchor_resolution")].copy()
    new = m.loc[m["resolution_status"].astype(str).eq("resolved_by_api")].copy()
    unresolved = m.loc[~m["resolution_status"].astype(str).isin(["reused_species_anchor_resolution", "resolved_by_api"])].copy()
    if len(reuse) != 26904 or len(new) != 58431 or len(unresolved) != 2:
        raise RuntimeError(f"frozen accounting drift reuse={len(reuse)} new={len(new)} unresolved={len(unresolved)}")

    # Reuse is keyed only by durable image identity; no reused image may be opened again.
    breadth_cols = [
        "observation_id", "photo_id", "morph", "measurement_status", "roi_status", "failure_reasons", "flower_effective_pixels"
    ]
    missing_b = [c for c in breadth_cols if c not in b.columns]
    if missing_b:
        raise RuntimeError(f"breadth result lacks reusable terminal columns: {missing_b}")
    reusable = b[breadth_cols].drop_duplicates(["observation_id", "photo_id"])
    reuse_join = reuse.merge(reusable, on=["observation_id", "photo_id"], how="left", validate="one_to_one")
    if reuse_join["measurement_status"].isna().any():
        raise RuntimeError("one or more frozen reuse rows lack a complete breadth terminal measurement")

    # New measurements use the same blind interface as the qualified breadth pipeline.
    if new["photo_url_large"].fillna("").astype(str).eq("").any():
        raise RuntimeError("new map measurement row lacks resolved URL")
    if new["photo_id"].nunique() != len(new) or new["observation_id"].nunique() != len(new):
        raise RuntimeError("new map measurement identities are not unique")
    if set(zip(new["observation_id"].astype(int), new["photo_id"].astype(int))) & set(zip(reuse["observation_id"].astype(int), reuse["photo_id"].astype(int))):
        raise RuntimeError("new and reused map image identities overlap")

    salt = str(contract["blinding"]["measurement_id_salt"])
    new["measurement_id"] = [measurement_id(int(v), salt) for v in new["photo_id"]]
    new["image_filename"] = new["measurement_id"].astype(str) + ".jpg"
    if new["measurement_id"].nunique() != len(new):
        raise RuntimeError("new map measurement IDs are not unique")

    worker = new[["measurement_id", "image_filename", "photo_license"]].copy()
    acquisition = new[["measurement_id", "image_filename", "photo_url_large", "photo_license"]].copy()
    if tuple(worker.columns) != WORKER_FIELDS or tuple(acquisition.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("blind worker/acquisition interface drift")

    meta_cols = ["measurement_id", "species", "inat_taxon_id", "cell_id", "observation_id", "photo_id", "discovery_source", "resolution_status", "photo_license", "attribution"]
    missing_m = [c for c in meta_cols if c not in new.columns]
    if missing_m:
        raise RuntimeError(f"new map rows lack metadata columns: {missing_m}")
    metadata_join = new[meta_cols].copy()
    if "photo_url_large" in metadata_join.columns:
        raise RuntimeError("sealed metadata join leaked URL")

    n_sem = int(contract["partitioning"]["semantic_shards"])
    n_comp = int(contract["partitioning"]["compute_partitions_per_semantic_shard"])
    assignments = pd.DataFrame({
        "measurement_id": new["measurement_id"].astype(str),
        "semantic_shard": [semantic_shard(v, n_sem) for v in new["measurement_id"].astype(str)],
        "compute_partition": [compute_partition(v, n_comp) for v in new["measurement_id"].astype(str)],
    })
    counts = assignments.groupby(["semantic_shard", "compute_partition"], observed=True).size()
    if len(counts) != 128:
        raise RuntimeError("unexpected empty map measurement partition")

    (OUT / "worker_packet").mkdir(parents=True, exist_ok=True)
    (OUT / "sealed_keys").mkdir(parents=True, exist_ok=True)
    worker.to_csv(OUT / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
    acquisition.to_csv(OUT / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
    metadata_join.to_csv(OUT / "sealed_keys/metadata_join_key.csv", index=False, lineterminator="\n")
    reuse_join.to_csv(OUT / "sealed_keys/reused_breadth_terminal_rows.csv.gz", index=False, compression="gzip", lineterminator="\n")
    unresolved.to_csv(OUT / "sealed_keys/transport_unresolved_rows.csv", index=False, lineterminator="\n")
    assignments.to_csv(OUT / "partition_assignments.csv", index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "taxon_cell_measurement_firewall_frozen_before_new_pixels",
        "taxon_cell_denominator": 85337,
        "species_universe": 42111,
        "occupied_cells": 128,
        "reused_breadth_terminal_rows": len(reuse_join),
        "new_image_measurement_rows": len(new),
        "transport_unresolved_rows": len(unresolved),
        "new_terminal_partitions": 128,
        "minimum_partition_rows": int(counts.min()),
        "maximum_partition_rows": int(counts.max()),
        "reused_image_pixels_reopened": False,
        "new_image_pixels_opened": False,
        "taxon_cell_colour_join_opened": False,
        "worker_contains_species": False,
        "worker_contains_cell": False,
        "worker_contains_source_url": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "transport_table_sha256": sha256_file(MAP_TABLE),
            "breadth_table_sha256": sha256_file(BREADTH_TABLE),
            "worker_sha256": sha256_file(OUT / "worker_packet/measurement_manifest.csv"),
            "acquisition_sha256": sha256_file(OUT / "sealed_keys/acquisition_key.csv"),
            "metadata_join_sha256": sha256_file(OUT / "sealed_keys/metadata_join_key.csv"),
        }
    }
    (OUT / "firewall_manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
