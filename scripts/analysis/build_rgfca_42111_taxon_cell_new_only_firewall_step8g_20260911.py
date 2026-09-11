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
CONTRACT = ROOT / "docs/supporting/rgfca_42111_taxon_cell_measurement_contract_v2.json"
MAP_RESULT = ROOT / "results/rgfca_42111_taxon_cell_resolution_step8e3_20260911/result.json"
MAP_TABLE = ROOT / "results/rgfca_42111_taxon_cell_resolution_step8e3_20260911/taxon_cell_anchor_metadata_85337.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxon_cell_new_firewall_step8g_20260911"


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

    if contract.get("status") != "frozen_before_any_taxon_cell_specific_image_pixel_v2":
        raise RuntimeError("taxon-cell v2 contract status drift")
    if contract.get("supersedes", {}).get("scientific_estimand_changed") is not False:
        raise RuntimeError("v2 contract does not preserve the scientific estimand")
    if mr.get("status") != "complete_metadata_only_taxon_cell_resolution":
        raise RuntimeError("taxon-cell transport is incomplete")
    if mr.get("image_pixels_opened") is not False or mr.get("flower_colour_used") is not False:
        raise RuntimeError("taxon-cell transport opened forbidden outcomes")
    if int(mr.get("taxon_cell_pairs", -1)) != 85337 or int(mr.get("species", -1)) != 42111 or int(mr.get("occupied_cells", -1)) != 128:
        raise RuntimeError("taxon-cell transport denominator drift")
    if sha256_file(MAP_TABLE) != contract["transport"]["required_source_table_sha256"]:
        raise RuntimeError("taxon-cell transport table hash drift")

    m = pd.read_csv(MAP_TABLE)
    if len(m) != 85337 or m[["inat_taxon_id", "cell_id"]].drop_duplicates().shape[0] != 85337:
        raise RuntimeError("map table is not exact one row per frozen taxon-cell")

    reuse = m.loc[m["resolution_status"].astype(str).eq("reused_species_anchor_resolution")].copy()
    new = m.loc[m["resolution_status"].astype(str).eq("resolved_by_api")].copy()
    unresolved = m.loc[~m["resolution_status"].astype(str).isin(["reused_species_anchor_resolution", "resolved_by_api"])].copy()
    if len(reuse) != 26904 or len(new) != 58431 or len(unresolved) != 2:
        raise RuntimeError(f"frozen accounting drift reuse={len(reuse)} new={len(new)} unresolved={len(unresolved)}")

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

    meta_cols = [
        "species", "inat_taxon_id", "cell_id", "observation_id", "photo_id",
        "discovery_source", "resolution_status", "photo_license", "attribution"
    ]
    missing_meta = [c for c in meta_cols if c not in m.columns]
    if missing_meta:
        raise RuntimeError(f"map rows lack required metadata columns: {missing_meta}")

    metadata_join = new[["measurement_id", *meta_cols]].copy()
    reuse_key = reuse[meta_cols].copy()
    unresolved_key = unresolved[meta_cols].copy()
    for frame_name, frame in {
        "metadata_join": metadata_join,
        "reuse_key": reuse_key,
        "unresolved_key": unresolved_key,
    }.items():
        forbidden = {"photo_url_large", "latitude", "longitude"} & set(frame.columns)
        if forbidden:
            raise RuntimeError(f"{frame_name} leaked forbidden transport fields: {sorted(forbidden)}")

    n_sem = int(contract["partitioning"]["semantic_shards"])
    n_comp = int(contract["partitioning"]["compute_partitions_per_semantic_shard"])
    assignments = pd.DataFrame({
        "measurement_id": new["measurement_id"].astype(str),
        "semantic_shard": [semantic_shard(v, n_sem) for v in new["measurement_id"].astype(str)],
        "compute_partition": [compute_partition(v, n_comp) for v in new["measurement_id"].astype(str)],
    })
    counts = assignments.groupby(["semantic_shard", "compute_partition"], observed=True).size()
    if len(counts) != 128:
        raise RuntimeError("unexpected empty new-only taxon-cell partition")
    if int(counts.sum()) != 58431:
        raise RuntimeError("new-only partition denominator drift")

    (OUT / "worker_packet").mkdir(parents=True, exist_ok=True)
    (OUT / "sealed_keys").mkdir(parents=True, exist_ok=True)
    worker.to_csv(OUT / "worker_packet/measurement_manifest.csv", index=False, lineterminator="\n")
    acquisition.to_csv(OUT / "sealed_keys/acquisition_key.csv", index=False, lineterminator="\n")
    metadata_join.to_csv(OUT / "sealed_keys/new_metadata_join_key.csv.gz", index=False, compression="gzip", lineterminator="\n")
    reuse_key.to_csv(OUT / "sealed_keys/reuse_deferred_key.csv.gz", index=False, compression="gzip", lineterminator="\n")
    unresolved_key.to_csv(OUT / "sealed_keys/transport_unresolved_key.csv", index=False, lineterminator="\n")
    assignments.to_csv(OUT / "partition_assignments.csv", index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "taxon_cell_new_only_firewall_frozen_before_pixels_v2",
        "taxon_cell_denominator": 85337,
        "species_universe": 42111,
        "occupied_cells": 128,
        "reused_breadth_rows_deferred": int(len(reuse_key)),
        "new_image_measurement_rows": int(len(new)),
        "transport_unresolved_rows": int(len(unresolved_key)),
        "new_terminal_partitions": 128,
        "minimum_partition_rows": int(counts.min()),
        "maximum_partition_rows": int(counts.max()),
        "breadth_outcomes_required_to_measure_new_images": False,
        "reused_image_pixels_reopened": False,
        "new_image_pixels_opened": False,
        "taxon_cell_colour_join_opened": False,
        "worker_contains_species": False,
        "worker_contains_cell": False,
        "worker_contains_source_url": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "transport_table_sha256": sha256_file(MAP_TABLE),
            "worker_sha256": sha256_file(OUT / "worker_packet/measurement_manifest.csv"),
            "acquisition_sha256": sha256_file(OUT / "sealed_keys/acquisition_key.csv"),
            "new_metadata_join_sha256": sha256_file(OUT / "sealed_keys/new_metadata_join_key.csv.gz"),
            "reuse_deferred_sha256": sha256_file(OUT / "sealed_keys/reuse_deferred_key.csv.gz"),
            "unresolved_sha256": sha256_file(OUT / "sealed_keys/transport_unresolved_key.csv"),
        },
        "claim_boundary": "This firewall authorizes only blind measurement of the 58,431 pre-resolved new taxon-cell photos. No species/cell colour join is opened."
    }
    (OUT / "firewall_manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
