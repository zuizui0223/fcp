#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h7_balanced_itv_contract_v1.json"
AUDIT = ROOT / "docs/supporting/random_photo_first_h7_metadata_feasibility_result_v1.json"
POOL = ROOT / "data/frozen/random_photo_first_candidate_pool_v1.csv"
OUT = ROOT / "data/frozen/random_photo_first_h7_species_cells_v1.csv"
MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_target_manifest_v1.json"
EARTH_RADIUS_KM = 6371.0088


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def great_circle_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    cosang = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return float(EARTH_RADIUS_KM * math.acos(max(-1.0, min(1.0, cosang))))


def select_maximin(cell_ids: list[int], centers: dict[int, tuple[float, float]], k: int) -> list[int]:
    ids = sorted({int(x) for x in cell_ids})
    if len(ids) < k:
        raise ValueError("not enough eligible cells")
    if k < 1:
        raise ValueError("k must be positive")
    distance: dict[tuple[int, int], float] = {}
    for a in ids:
        for b in ids:
            if (a, b) not in distance:
                la, loa = centers[a]
                lb, lob = centers[b]
                d = great_circle_km(la, loa, lb, lob)
                distance[(a, b)] = distance[(b, a)] = d

    mean_distance = {a: float(np.mean([distance[(a, b)] for b in ids])) for a in ids}
    first = sorted(ids, key=lambda a: (-mean_distance[a], a))[0]
    selected = [first]
    while len(selected) < k:
        remaining = [a for a in ids if a not in selected]
        score = {a: min(distance[(a, b)] for b in selected) for a in remaining}
        chosen = sorted(remaining, key=lambda a: (-score[a], a))[0]
        selected.append(chosen)
    return selected


def main() -> None:
    contract = json.loads(CONTRACT.read_text())
    audit = json.loads(AUDIT.read_text())
    if contract["status"] != "prospective_h7_frozen_after_metadata_feasibility_before_any_fresh_h7_query_or_pixel":
        raise RuntimeError("H7 contract is not frozen")
    if audit["status"] != "complete_metadata_only_feasibility_audit":
        raise RuntimeError("H7 metadata audit is incomplete")
    if audit["proposed_target_check"]["qualifying_species"] != 52:
        raise RuntimeError("unexpected metadata feasibility result")
    if sha256_file(POOL) != contract["species_discovery"]["required_source_sha256"]:
        raise RuntimeError("discovery pool lineage mismatch")

    df = pd.read_csv(POOL, usecols=["species", "inat_taxon_id", "cell_id", "observation_id", "observer_id"])
    counts = (
        df.groupby(["species", "inat_taxon_id", "cell_id"], observed=True)
        .agg(old_metadata_records=("observation_id", "size"), old_unique_observers=("observer_id", "nunique"))
        .reset_index()
    )
    eligible = counts[counts["old_metadata_records"] >= 2].copy()
    n_cells = eligible.groupby(["species", "inat_taxon_id"], observed=True)["cell_id"].nunique()
    qualifying_keys = n_cells[n_cells >= 6].index
    if len(qualifying_keys) != int(contract["species_discovery"]["expected_species"]):
        raise RuntimeError(f"expected 52 qualifying species, found {len(qualifying_keys)}")

    grid = EqualAreaGrid(
        n_lon=int(contract["target_cells"]["grid"]["n_lon"]),
        n_sinlat=int(contract["target_cells"]["grid"]["n_sinlat"]),
    )
    ids, lats, lons = equal_area_cell_centers(grid)
    centers = {int(i): (float(la), float(lo)) for i, la, lo in zip(ids, lats, lons)}

    rows: list[dict[str, object]] = []
    k = int(contract["target_cells"]["cells_per_species"])
    for species, taxon_id in sorted(qualifying_keys, key=lambda x: (int(x[1]), str(x[0]))):
        sub = eligible[(eligible["species"] == species) & (eligible["inat_taxon_id"] == taxon_id)].copy()
        chosen = select_maximin(sub["cell_id"].astype(int).tolist(), centers, k)
        for order, cell_id in enumerate(chosen, start=1):
            row = sub[sub["cell_id"].astype(int) == int(cell_id)].iloc[0]
            lat, lon = centers[int(cell_id)]
            min_prev = None
            if order > 1:
                min_prev = min(
                    great_circle_km(lat, lon, centers[int(prev)][0], centers[int(prev)][1])
                    for prev in chosen[: order - 1]
                )
            rows.append({
                "species": str(species),
                "inat_taxon_id": int(taxon_id),
                "target_cell_order": int(order),
                "cell_id": int(cell_id),
                "old_metadata_records": int(row["old_metadata_records"]),
                "old_unique_observers": int(row["old_unique_observers"]),
                "cell_center_latitude": lat,
                "cell_center_longitude": lon,
                "min_distance_to_previous_selected_km": None if min_prev is None else float(min_prev),
            })

    out = pd.DataFrame(rows).sort_values(["inat_taxon_id", "target_cell_order"], kind="mergesort")
    if out["inat_taxon_id"].nunique() != 52 or len(out) != 312:
        raise RuntimeError("target freeze did not produce 52 species x 6 cells")
    if not (out.groupby("inat_taxon_id")["cell_id"].nunique() == 6).all():
        raise RuntimeError("duplicate target cells within species")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    manifest = {
        "protocol": contract["protocol"],
        "status": "h7_species_and_target_cells_frozen_before_fresh_query",
        "species": 52,
        "target_cells_per_species": 6,
        "species_cell_targets": 312,
        "selection_rule": contract["target_cells"]["selection"],
        "source_pool_sha256": sha256_file(POOL),
        "target_table_sha256": sha256_file(OUT),
        "outcome_firewall": {
            "colour_columns_read": False,
            "measured_table_read": False,
            "h6_or_h6b_effects_read": False,
            "fresh_h7_api_queries_opened": False,
            "fresh_h7_pixels_opened": False,
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
