#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h9_metadata_feasibility_contract_v1.json"
H7 = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_v1.csv"
H7_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_fresh_metadata_manifest_v1.json"
H8_RESULT = ROOT / "docs/supporting/random_photo_first_h8_design_audit_result_v1.json"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h9_metadata_feasibility_result_v1.json"
OUT_SPECIES = ROOT / "data/frozen/random_photo_first_h9_metadata_species_v1.csv"

EARTH_RADIUS_KM = 6371.0088


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def max_span_km(lat: np.ndarray, lon: np.ndarray) -> float:
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if len(lat) < 2:
        return 0.0
    phi = np.deg2rad(lat)
    lam = np.deg2rad(lon)
    cos_phi = np.cos(phi)
    xyz = np.column_stack([cos_phi * np.cos(lam), cos_phi * np.sin(lam), np.sin(phi)])
    dots = np.clip(xyz @ xyz.T, -1.0, 1.0)
    angular = np.arccos(dots)
    return float(np.nanmax(angular) * EARTH_RADIUS_KM)


def main() -> None:
    contract = json.loads(CONTRACT.read_text())
    h7m = json.loads(H7_MANIFEST.read_text())
    h8 = json.loads(H8_RESULT.read_text())
    if bool(h7m["candidate_image_pixels_opened"]) != bool(contract["source"]["required_h7_pixels_opened"]):
        raise RuntimeError("H7 pixel-opening precondition mismatch")
    if bool(h8["automatic_selection"]["h8_metadata_design_feasible"]) != bool(contract["source"]["required_h8_feasible"]):
        raise RuntimeError("H8 feasibility precondition mismatch")
    expected_h7_sha = str(h7m["files"]["observations"]["sha256"])
    if sha256_file(H7) != expected_h7_sha:
        raise RuntimeError("H7 fresh metadata lineage mismatch")

    df = pd.read_csv(H7)
    required = {"species", "inat_taxon_id", "observation_id", "photo_id", "latitude", "longitude", "cell_id"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"H7 fresh metadata missing columns: {sorted(missing)}")
    if len(df) != int(h7m["retained_fresh_photos"]):
        raise RuntimeError("H7 fresh metadata row-count mismatch")
    if df["observation_id"].duplicated().any() or df["photo_id"].duplicated().any():
        raise RuntimeError("H7 fresh metadata IDs are not unique")

    rows = []
    for (species, taxon_id), group in df.groupby(["species", "inat_taxon_id"], sort=True, observed=True):
        rows.append({
            "species": str(species),
            "inat_taxon_id": int(taxon_id),
            "h7_fresh_records": int(len(group)),
            "h7_occupied_target_cells": int(group["cell_id"].nunique()),
            "h7_maximum_span_km": max_span_km(group["latitude"].to_numpy(), group["longitude"].to_numpy()),
            "h7_unique_observers": int(group["observer_id"].nunique()) if "observer_id" in group.columns else 0,
        })
    species = pd.DataFrame(rows).sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)

    counts = [int(x) for x in contract["audit_matrix"]["minimum_h7_fresh_records_per_species"]]
    spans = [float(x) for x in contract["audit_matrix"]["minimum_h7_fresh_maximum_span_km"]]
    minimum_species = int(contract["audit_matrix"]["minimum_species_required"])
    matrix = []
    for count in counts:
        for span in spans:
            qualifies = (species["h7_fresh_records"] >= count) & (species["h7_maximum_span_km"] >= span)
            n = int(qualifies.sum())
            matrix.append({
                "minimum_h7_fresh_records_per_species": count,
                "minimum_h7_fresh_maximum_span_km": span,
                "qualifying_species": n,
                "meets_minimum_species": bool(n >= minimum_species),
            })

    eligible = [row for row in matrix if row["meets_minimum_species"]]
    selected_rule = None
    selected = species.iloc[0:0].copy()
    if eligible:
        selected_rule = sorted(
            eligible,
            key=lambda row: (
                int(row["minimum_h7_fresh_records_per_species"]),
                float(row["minimum_h7_fresh_maximum_span_km"]),
            ),
            reverse=True,
        )[0]
        selected = species.loc[
            (species["h7_fresh_records"] >= int(selected_rule["minimum_h7_fresh_records_per_species"]))
            & (species["h7_maximum_span_km"] >= float(selected_rule["minimum_h7_fresh_maximum_span_km"]))
        ].copy()
    selected["h9_selected_from_metadata_only"] = True
    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(OUT_SPECIES, index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "complete_h9_metadata_only_feasibility_audit",
        "source_h7_fresh_photos": int(len(df)),
        "source_species": int(len(species)),
        "matrix": matrix,
        "automatic_selection": {
            "minimum_species_required": minimum_species,
            "selected_rule": selected_rule,
            "selected_species": int(len(selected)),
            "h9_metadata_feasible": bool(selected_rule is not None),
        },
        "prospective_h9_design_if_feasible": contract["prospective_h9_if_feasible"],
        "outcome_firewall": {
            "h7_colour_read": False,
            "h6_h6b_effects_read": False,
            "h9_api_queries_opened": False,
            "h9_pixels_opened": False,
        },
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "h7_fresh_metadata_sha256": sha256_file(H7),
            "h7_manifest_sha256": sha256_file(H7_MANIFEST),
            "h8_design_result_sha256": sha256_file(H8_RESULT),
            "selected_species_sha256": sha256_file(OUT_SPECIES),
        },
        "claim_ceiling": contract["claim_ceiling"],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
