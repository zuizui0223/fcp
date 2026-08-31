#!/usr/bin/env python3
"""Freeze global environmental boundary surfaces before atlas colour joins."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import (
    environmental_boundary_surfaces,
    rook_adjacency_without_repair,
    validate_environment_contract,
)
from fcp_pipeline.atlas_expansion import validate_expansion_contract


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_argument(value: str) -> tuple[int, Path]:
    try:
        scale, path = value.split("=", 1)
        return int(scale), Path(path)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("grid must be SCALE_KM=PATH") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json"),
    )
    parser.add_argument(
        "--expansion-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"),
    )
    parser.add_argument(
        "--atlas-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--grid", action="append", type=grid_argument, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    environment = json.loads(args.environment_contract.read_text(encoding="utf-8"))
    expansion = json.loads(args.expansion_contract.read_text(encoding="utf-8"))
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    source_manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    validate_environment_contract(environment)
    validate_expansion_contract(expansion)
    if source_manifest.get("status") != "pass_environmental_source_freeze":
        raise RuntimeError("environmental source manifest has not passed")
    grids = dict(args.grid)
    expected_scales = environment["grid"]["scales_km"]
    if sorted(grids) != expected_scales:
        raise RuntimeError("all and only the frozen 100/250/500-km grids are required")
    dimensions = {
        int(row["scale_km"]): (int(row["n_lon"]), int(row["n_sinlat"]))
        for row in atlas["geometry_only_scale_selection"]["candidates"]
    }

    output_files: dict[str, str] = {}
    scale_results: list[dict[str, Any]] = []
    for scale in expected_scales:
        rows = read_csv(grids[scale])
        cell_ids = [int(row["cell_id"]) for row in rows]
        if len(cell_ids) != len(set(cell_ids)):
            raise RuntimeError(f"{scale}-km environment grid has duplicate cells")
        n_lon, n_sinlat = dimensions[scale]
        adjacency = rook_adjacency_without_repair(
            cell_ids, n_lon=n_lon, n_sinlat=n_sinlat
        )
        available = tuple(source_manifest["available_primary_families"])
        surfaces = environmental_boundary_surfaces(rows, adjacency, families=available)
        output_rows: list[dict[str, Any]] = []
        families = tuple(available) + (
            ("realm_sensitivity", "biome_sensitivity")
            if "ecoregion" in available
            else ()
        )
        for index, row in enumerate(rows):
            output = {
                "scale_km": scale,
                "cell_id": cell_ids[index],
                "latitude": float(row["latitude"]),
                "longitude": float(row["longitude"]),
            }
            for family in families:
                value = float(surfaces[family][index])
                output[f"{family}_boundary"] = value if math.isfinite(value) else ""
            output_rows.append(output)
        path = args.output_dir / f"environmental_boundary_cells_{scale}km.csv"
        write_csv(path, output_rows)
        output_files[path.name] = sha256(path)
        scale_results.append(
            {
                "scale_km": scale,
                "terrestrial_cells": len(rows),
                "finite_cells_by_family": {
                    family: sum(
                        math.isfinite(float(value)) for value in surfaces[family]
                    )
                    for family in families
                },
                "continuous_scaling": surfaces["scaling"],
            }
        )

    manifest = {
        "status": "pass_precolour_environmental_boundary_freeze",
        "protocol": environment["protocol"],
        "scaleout_colour_opened": False,
        "environment_colour_join_performed": False,
        "source_manifest_sha256": sha256(args.source_manifest),
        "input_grid_sha256": {
            f"{scale}km": sha256(path) for scale, path in sorted(grids.items())
        },
        "output_sha256": output_files,
        "scale_results": scale_results,
        "claim_ceiling": "Independent environmental boundary surfaces only; no flower-colour concordance has been tested.",
    }
    manifest_path = args.output_dir / "environmental_boundary_freeze_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
