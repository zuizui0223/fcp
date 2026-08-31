#!/usr/bin/env python3
"""Validate committed pre-colour environmental boundary surfaces."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path("data/atlas/environment")
MANIFEST = ROOT / "environmental_boundary_freeze_manifest.json"
SOURCE_MANIFEST = ROOT / "environmental_source_freeze_manifest.json"
EXPECTED_HASHES = {
    "environmental_boundary_cells_100km.csv": "3b9418a5e11d628facc702adfd026bce2221971231aa129012879d35a9e069b1",
    "environmental_boundary_cells_250km.csv": "0d6715c18003f3945537d4a16ec5bff5f974325a369a551e1916073e707bf50b",
    "environmental_boundary_cells_500km.csv": "25bdddc98e9db12101729fb152be214778c5ad3b989c191217bcb24b65d6c6cc",
}
EXPECTED_ROWS = {100: 17526, 250: 3250, 500: 961}
EXPECTED_FINITE = {
    100: {"macroclimate": 15742, "land_cover": 17308, "ecoregion": 15742},
    250: {"macroclimate": 3036, "land_cover": 3165, "ecoregion": 3036},
    500: {"macroclimate": 940, "land_cover": 877, "ecoregion": 940},
}
BOUNDARY_FIELDS = (
    "macroclimate_boundary",
    "land_cover_boundary",
    "ecoregion_boundary",
    "realm_sensitivity_boundary",
    "biome_sensitivity_boundary",
)
EXPECTED_FIELDS = {"scale_km", "cell_id", "latitude", "longitude", *BOUNDARY_FIELDS}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    source = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
    if (
        manifest.get("status") != "pass_precolour_environmental_boundary_freeze"
        or manifest.get("scaleout_colour_opened") is not False
        or manifest.get("environment_colour_join_performed") is not False
        or source.get("status") != "pass_environmental_source_freeze"
        or source.get("available_primary_families")
        != ["macroclimate", "land_cover", "ecoregion"]
        or source.get("not_evaluable_primary_families") != ["terrain"]
        or manifest.get("source_manifest_sha256") != sha256(SOURCE_MANIFEST)
        or manifest.get("output_sha256") != EXPECTED_HASHES
    ):
        raise RuntimeError("environmental boundary freeze identity changed")
    verified = []
    for scale in (100, 250, 500):
        path = ROOT / f"environmental_boundary_cells_{scale}km.csv"
        if sha256(path) != EXPECTED_HASHES[path.name]:
            raise RuntimeError(f"{scale}-km environmental boundary hash changed")
        rows = read_csv(path)
        ids = [int(row["cell_id"]) for row in rows]
        if (
            len(rows) != EXPECTED_ROWS[scale]
            or ids != sorted(ids)
            or len(ids) != len(set(ids))
            or not rows
            or set(rows[0]) != EXPECTED_FIELDS
        ):
            raise RuntimeError(f"{scale}-km environmental boundary structure changed")
        counts = {family: 0 for family in EXPECTED_FINITE[scale]}
        for row in rows:
            if int(row["scale_km"]) != scale:
                raise RuntimeError("environmental row assigned to wrong scale")
            coordinates = (float(row["latitude"]), float(row["longitude"]))
            if not all(math.isfinite(value) for value in coordinates):
                raise RuntimeError("environmental grid has non-finite coordinates")
            finite_any = False
            for field in BOUNDARY_FIELDS:
                if row[field] == "":
                    continue
                value = float(row[field])
                if not math.isfinite(value) or value < 0:
                    raise RuntimeError("environmental boundary is invalid")
                if field != "macroclimate_boundary" and value > 1:
                    raise RuntimeError("categorical/composition boundary exceeds one")
                finite_any = True
                family = field.removesuffix("_boundary")
                if family in counts:
                    counts[family] += 1
            if not finite_any:
                raise RuntimeError("union grid retained a cell with no finite family boundary")
        if counts != EXPECTED_FINITE[scale]:
            raise RuntimeError(f"{scale}-km family denominator changed")
        verified.append(
            {
                "scale_km": scale,
                "union_cells": len(rows),
                "finite_primary_cells": counts,
                "sha256": sha256(path),
            }
        )
    print(
        json.dumps(
            {
                "status": "pass_committed_precolour_environmental_boundary_freeze",
                "available_primary_families": source["available_primary_families"],
                "not_evaluable_primary_families": source["not_evaluable_primary_families"],
                "scaleout_colour_opened": False,
                "environment_colour_join_performed": False,
                "verified_grids": verified,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
