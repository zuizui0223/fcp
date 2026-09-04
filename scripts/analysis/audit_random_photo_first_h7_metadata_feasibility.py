#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h7_metadata_feasibility_contract_v1.json"
POOL = ROOT / "data/frozen/random_photo_first_candidate_pool_v1.csv"
MANIFEST = ROOT / "docs/supporting/random_photo_first_candidate_pool_manifest_v1.json"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h7_metadata_feasibility_result_v1.json"
OUT_CSV = ROOT / "data/derived/random_photo_first_h7_metadata_feasibility_species_v1.csv"

ALLOWED = {
    "species", "inat_taxon_id", "cell_id", "observation_id", "photo_id",
    "observer_id", "latitude", "longitude"
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    contract = json.loads(CONTRACT.read_text())
    manifest = json.loads(MANIFEST.read_text())
    if manifest["candidate_table_sha256"] != contract["source"]["required_candidate_sha256"]:
        raise RuntimeError("candidate pool lineage mismatch")
    if sha256_file(POOL) != contract["source"]["required_candidate_sha256"]:
        raise RuntimeError("candidate pool file hash mismatch")

    df = pd.read_csv(POOL)
    missing = ALLOWED - set(df.columns)
    if missing:
        raise RuntimeError(f"missing metadata columns: {sorted(missing)}")
    # Hard firewall: this audit uses only columns that existed before image pixels opened.
    df = df[sorted(ALLOWED)].copy()
    if len(df) != int(manifest["counts"]["observations"]):
        raise RuntimeError("candidate row count mismatch")
    if df["observation_id"].duplicated().any() or df["photo_id"].duplicated().any():
        raise RuntimeError("candidate IDs are not unique")

    cell_counts = (
        df.groupby(["species", "inat_taxon_id", "cell_id"], observed=True)
        .agg(records=("observation_id", "size"), observers=("observer_id", "nunique"))
        .reset_index()
    )
    species = (
        cell_counts.groupby(["species", "inat_taxon_id"], observed=True)
        .agg(
            total_metadata_records=("records", "sum"),
            occupied_cells=("cell_id", "nunique"),
            cells_ge_2=("records", lambda s: int((s >= 2).sum())),
            cells_ge_3=("records", lambda s: int((s >= 3).sum())),
            max_records_in_cell=("records", "max"),
            median_records_per_occupied_cell=("records", "median"),
        )
        .reset_index()
    )
    species = species.sort_values(
        ["cells_ge_2", "occupied_cells", "total_metadata_records", "inat_taxon_id"],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    species["metadata_rank"] = species.index + 1
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    species.to_csv(OUT_CSV, index=False)

    matrix = []
    for min_records in contract["audit_matrix"]["minimum_old_metadata_records_per_species_cell"]:
        per_species = (
            cell_counts.assign(ok=cell_counts["records"] >= int(min_records))
            .groupby(["species", "inat_taxon_id"], observed=True)["ok"]
            .sum()
        )
        for min_cells in contract["audit_matrix"]["minimum_qualifying_cells_per_species"]:
            matrix.append({
                "minimum_old_metadata_records_per_species_cell": int(min_records),
                "minimum_qualifying_cells_per_species": int(min_cells),
                "qualifying_species": int((per_species >= int(min_cells)).sum()),
            })

    target = contract["prospective_design_target_to_assess"]
    target_feasible_species = next(
        row["qualifying_species"] for row in matrix
        if row["minimum_old_metadata_records_per_species_cell"] == 2
        and row["minimum_qualifying_cells_per_species"] == int(target["fresh_target_cells_per_species"])
    )
    result = {
        "protocol": contract["protocol"],
        "status": "complete_metadata_only_feasibility_audit",
        "source_candidate_rows": int(len(df)),
        "source_species": int(species.shape[0]),
        "matrix": matrix,
        "proposed_target_check": {
            "reference_rule": "at least 2 old metadata records in each of at least 6 cells",
            "qualifying_species": int(target_feasible_species),
            "target_species": int(target["target_species"]),
            "target_is_metadata_feasible": bool(target_feasible_species >= int(target["target_species"])),
        },
        "outcome_firewall": {
            "colour_columns_read": False,
            "measured_table_read": False,
            "h1_h6_species_effects_read": False,
            "fresh_h7_api_queries_opened": False,
            "fresh_h7_image_pixels_opened": False,
        },
        "lineage": {
            "candidate_table_sha256": sha256_file(POOL),
            "candidate_manifest_blob_status": manifest["status"],
        },
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
