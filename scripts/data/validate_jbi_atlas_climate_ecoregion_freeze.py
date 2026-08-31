#!/usr/bin/env python3
"""Validate committed pre-colour WorldClim and RESOLVE equal-area grids."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path("data/atlas/environment")
MANIFEST = ROOT / "climate_ecoregion_source_manifest.json"
EXPECTED_ROWS = {100: 15833, 250: 3076, 500: 956}
EXPECTED_GLOBAL = {100: 51200, 250: 8192, 500: 2048}
EXPECTED_FIELDS = {
    "scale_km",
    "cell_id",
    "latitude",
    "longitude",
    "bio1",
    "bio4",
    "bio12",
    "bio15",
    "realm",
    "biome",
    "ecoregion",
    "dominant_eco_id",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "pass_climate_ecoregion_source_freeze":
        raise RuntimeError("climate-ecoregion source freeze has not passed")
    if manifest.get("scaleout_colour_opened") is not False:
        raise RuntimeError("source freeze opened scale-out colour")
    if manifest.get("available_primary_families") != ["macroclimate", "ecoregion"]:
        raise RuntimeError("source family set changed")
    results = {int(row["scale_km"]): row for row in manifest["scale_results"]}
    verified = []
    for scale, expected_rows in EXPECTED_ROWS.items():
        path = ROOT / f"climate_ecoregion_grid_{scale}km.csv"
        if sha256(path) != manifest["output_sha256"][path.name]:
            raise RuntimeError(f"{scale}-km source grid hash mismatch")
        rows = read_csv(path)
        if len(rows) != expected_rows or len({row["cell_id"] for row in rows}) != len(rows):
            raise RuntimeError(f"{scale}-km source grid denominator changed")
        if set(rows[0]) != EXPECTED_FIELDS:
            raise RuntimeError(f"{scale}-km source grid fields changed")
        if int(results[scale]["global_cells"]) != EXPECTED_GLOBAL[scale]:
            raise RuntimeError(f"{scale}-km global grid dimensions changed")
        for row in rows:
            if int(row["scale_km"]) != scale:
                raise RuntimeError("source row assigned to wrong scale")
            if not 0 <= int(row["cell_id"]) < EXPECTED_GLOBAL[scale]:
                raise RuntimeError("source cell ID is outside the grid")
            numeric = [float(row[key]) for key in ("latitude", "longitude", "bio1", "bio4", "bio12", "bio15")]
            if not all(math.isfinite(value) for value in numeric):
                raise RuntimeError("source grid contains non-finite climate")
            if not all(row[key].strip() for key in ("realm", "biome", "ecoregion")):
                raise RuntimeError("source grid contains an empty ecoregion label")
        if len({row["realm"] for row in rows}) != 8 or len({row["biome"] for row in rows}) != 14:
            raise RuntimeError("global realm or biome coverage changed")
        verified.append(
            {
                "scale_km": scale,
                "terrestrial_cells": len(rows),
                "sha256": sha256(path),
                "realms": 8,
                "biomes": 14,
            }
        )
    print(
        json.dumps(
            {
                "status": "pass_committed_climate_ecoregion_freeze",
                "scaleout_colour_opened": False,
                "source_pixels": manifest["source_pixels_with_climate_and_ecoregion"],
                "verified_grids": verified,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
