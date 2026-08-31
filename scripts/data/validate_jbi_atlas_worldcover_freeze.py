#!/usr/bin/env python3
"""Validate the committed WorldCover freeze and its sampling sensitivity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import (
    LAND_COVER_CODES,
    composition_boundary,
    rook_adjacency_without_repair,
    spearman_rank_correlation,
    validate_environment_contract,
)


DATA_ROOT = Path("data/atlas/environment")
EXPECTED_INVENTORY_SHA256 = "bc6bc85c8c5c92a376036449d2cbaafc1b1febe06ffa94f702461424602d36fe"
EXPECTED_OBJECTS = 2651
EXPECTED_TOTAL_SIZE_BYTES = 124027923380
EXPECTED_HASHES = {
    "worldcover_composition_primary_100km.csv": "31c2df0215bd3f0adec1d921e40514bf5faca3fe1943de29d83eab94dee27154",
    "worldcover_composition_primary_250km.csv": "e873b93e9694b44728e4901a97f7c33ffe147c2dec4eec4a220b991dfa60333f",
    "worldcover_composition_primary_500km.csv": "405e15577ffa7ea34eb251c58536d9e0506cc35f8965b404ba7acb407d39b57c",
    "worldcover_composition_sensitivity_100km.csv": "ca14ea4eaa583c33cc388593b15caed5084eb126095614b822ec10b542d6b47f",
    "worldcover_composition_sensitivity_250km.csv": "044e0ea56aa89d9de88060c2742d4021ef202e3dafff4edf719201a6c03092c2",
    "worldcover_composition_sensitivity_500km.csv": "a0c84d8e3322a0b06097bd9731d7aba6940f68c8943bbe8ac7647fbe4ca61d59",
}
EXPECTED_ROWS = {
    ("primary", 100): 17327,
    ("sensitivity", 100): 16404,
    ("primary", 250): 3187,
    ("sensitivity", 250): 2902,
    ("primary", 500): 891,
    ("sensitivity", 500): 771,
}
MINIMUM_SAMPLES = {100: 25, 250: 100, 500: 400}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json"),
    )
    parser.add_argument(
        "--atlas-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument(
        "--inventory", type=Path, default=DATA_ROOT / "worldcover_2021_v200_map_inventory.csv"
    )
    parser.add_argument(
        "--inventory-manifest",
        type=Path,
        default=DATA_ROOT / "worldcover_2021_v200_inventory_manifest.json",
    )
    parser.add_argument(
        "--composition-manifest",
        type=Path,
        default=DATA_ROOT / "worldcover_composition_manifest.json",
    )
    parser.add_argument("--composition-dir", type=Path, default=DATA_ROOT)
    parser.add_argument(
        "--result",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_worldcover_qualification_result_v1.json"),
    )
    parser.add_argument("--write-result", action="store_true")
    return parser.parse_args()


def validate_rows(path: Path, *, sampling: str, scale: int, n_cells: int) -> list[dict[str, str]]:
    if sha256(path) != EXPECTED_HASHES[path.name]:
        raise RuntimeError(f"WorldCover hash changed: {path.name}")
    rows = read_csv(path)
    if len(rows) != EXPECTED_ROWS[(sampling, scale)]:
        raise RuntimeError(f"WorldCover denominator changed: {path.name}")
    expected_fields = {
        "scale_km", "cell_id", "latitude", "longitude", "valid_samples",
        *(f"worldcover_{code}" for code in LAND_COVER_CODES),
    }
    if not rows or set(rows[0]) != expected_fields:
        raise RuntimeError(f"WorldCover fields changed: {path.name}")
    ids = [int(row["cell_id"]) for row in rows]
    if ids != sorted(ids) or len(ids) != len(set(ids)):
        raise RuntimeError(f"WorldCover cells are duplicated or unsorted: {path.name}")
    for row in rows:
        if int(row["scale_km"]) != scale or not 0 <= int(row["cell_id"]) < n_cells:
            raise RuntimeError(f"WorldCover grid identity changed: {path.name}")
        if int(row["valid_samples"]) < MINIMUM_SAMPLES[scale]:
            raise RuntimeError(f"WorldCover sample minimum failed: {path.name}")
        coordinates = (float(row["latitude"]), float(row["longitude"]))
        composition = [float(row[f"worldcover_{code}"]) for code in LAND_COVER_CODES]
        if not all(math.isfinite(value) for value in (*coordinates, *composition)):
            raise RuntimeError(f"WorldCover contains non-finite values: {path.name}")
        if any(value < 0 or value > 1 for value in composition):
            raise RuntimeError(f"WorldCover composition is invalid: {path.name}")
        if not math.isclose(sum(composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"WorldCover composition does not sum to one: {path.name}")
    return rows


def boundary_by_cell(rows: list[dict[str, str]], *, n_lon: int, n_sinlat: int) -> dict[int, float]:
    cell_ids = np.asarray([int(row["cell_id"]) for row in rows], dtype=int)
    adjacency = rook_adjacency_without_repair(cell_ids, n_lon=n_lon, n_sinlat=n_sinlat)
    composition = np.asarray(
        [[float(row[f"worldcover_{code}"]) for code in LAND_COVER_CODES] for row in rows],
        dtype=float,
    )
    boundary = composition_boundary(composition, adjacency)
    return {
        int(cell): float(value)
        for cell, value in zip(cell_ids, boundary, strict=True)
        if math.isfinite(float(value))
    }


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    environment = json.loads(args.environment_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    inventory_manifest = json.loads(args.inventory_manifest.read_text(encoding="utf-8"))
    composition_manifest = json.loads(args.composition_manifest.read_text(encoding="utf-8"))
    validate_environment_contract(environment)
    if (
        sha256(args.inventory) != EXPECTED_INVENTORY_SHA256
        or inventory_manifest.get("inventory_sha256") != EXPECTED_INVENTORY_SHA256
        or inventory_manifest.get("objects") != EXPECTED_OBJECTS
        or inventory_manifest.get("total_size_bytes") != EXPECTED_TOTAL_SIZE_BYTES
        or inventory_manifest.get("scaleout_colour_opened") is not False
    ):
        raise RuntimeError("WorldCover inventory freeze changed")
    if (
        composition_manifest.get("status") != "pass_worldcover_composition_freeze"
        or composition_manifest.get("formal_complete_inventory") is not True
        or composition_manifest.get("tiles_processed") != EXPECTED_OBJECTS
        or composition_manifest.get("tile_failures") != 0
        or composition_manifest.get("scaleout_colour_opened") is not False
        or composition_manifest.get("output_sha256") != EXPECTED_HASHES
    ):
        raise RuntimeError("WorldCover formal composition manifest changed")
    dimensions = {
        int(row["scale_km"]): (int(row["n_lon"]), int(row["n_sinlat"]))
        for row in atlas["geometry_only_scale_selection"]["candidates"]
    }
    threshold = float(environment["worldcover_sampling"]["minimum_boundary_spearman"])
    scale_results = []
    all_pass = True
    for scale in environment["grid"]["scales_km"]:
        n_lon, n_sinlat = dimensions[scale]
        row_sets = {}
        boundaries = {}
        for sampling in ("primary", "sensitivity"):
            name = f"worldcover_composition_{sampling}_{scale}km.csv"
            rows = validate_rows(
                args.composition_dir / name,
                sampling=sampling,
                scale=scale,
                n_cells=n_lon * n_sinlat,
            )
            row_sets[sampling] = rows
            boundaries[sampling] = boundary_by_cell(rows, n_lon=n_lon, n_sinlat=n_sinlat)
        common = sorted(set(boundaries["primary"]) & set(boundaries["sensitivity"]))
        correlation = spearman_rank_correlation(
            np.asarray([boundaries["primary"][cell] for cell in common]),
            np.asarray([boundaries["sensitivity"][cell] for cell in common]),
        )
        passed = correlation >= threshold
        all_pass = all_pass and passed
        scale_results.append(
            {
                "scale_km": scale,
                "primary_cells": len(row_sets["primary"]),
                "sensitivity_cells": len(row_sets["sensitivity"]),
                "common_finite_boundary_cells": len(common),
                "boundary_spearman": round(correlation, 12),
                "minimum_required_spearman": threshold,
                "passed": passed,
            }
        )
    return {
        "protocol": "jbi-atlas-worldcover-sampling-qualification-v1",
        "status": "pass_land_cover_family_evaluable" if all_pass else "land_cover_family_not_evaluable",
        "inventory_objects": EXPECTED_OBJECTS,
        "inventory_total_size_bytes": EXPECTED_TOTAL_SIZE_BYTES,
        "inventory_sha256": EXPECTED_INVENTORY_SHA256,
        "source_resolution_m": 10,
        "primary_systematic_sample_approx_km": 5,
        "sensitivity_systematic_sample_approx_km": 10,
        "scale_results": scale_results,
        "scaleout_colour_opened": False,
        "environment_colour_join_performed": False,
        "claim_ceiling": "Sampling stability of pre-colour WorldCover boundary intensity only; not exhaustive 10 m aggregation and no flower-colour concordance.",
    }


def main() -> None:
    args = parse_args()
    result = evaluate(args)
    if args.write_result:
        args.result.parent.mkdir(parents=True, exist_ok=True)
        args.result.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        frozen = json.loads(args.result.read_text(encoding="utf-8"))
        if result != frozen:
            raise RuntimeError("committed WorldCover qualification result changed")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
