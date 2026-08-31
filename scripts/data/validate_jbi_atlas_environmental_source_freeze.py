#!/usr/bin/env python3
"""Build or validate the joined pre-colour environmental source freeze."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import validate_environment_contract


DATA_ROOT = Path("data/atlas/environment")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json"),
    )
    parser.add_argument(
        "--climate-manifest",
        type=Path,
        default=DATA_ROOT / "climate_ecoregion_source_manifest.json",
    )
    parser.add_argument(
        "--worldcover-manifest",
        type=Path,
        default=DATA_ROOT / "worldcover_composition_manifest.json",
    )
    parser.add_argument(
        "--worldcover-result",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_worldcover_qualification_result_v1.json"),
    )
    parser.add_argument(
        "--terrain-result",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_terrain_access_result_v1.json"),
    )
    parser.add_argument(
        "--manifest", type=Path, default=DATA_ROOT / "environmental_source_freeze_manifest.json"
    )
    parser.add_argument("--write", action="store_true")
    return parser.parse_args()


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    contract = load(args.environment_contract)
    climate = load(args.climate_manifest)
    worldcover = load(args.worldcover_manifest)
    worldcover_result = load(args.worldcover_result)
    terrain = load(args.terrain_result)
    validate_environment_contract(contract)
    if climate.get("status") != "pass_climate_ecoregion_source_freeze":
        raise RuntimeError("climate and ecoregion source freeze did not pass")
    if worldcover.get("status") != "pass_worldcover_composition_freeze":
        raise RuntimeError("WorldCover composition freeze did not pass")
    if worldcover_result.get("status") != "pass_land_cover_family_evaluable":
        raise RuntimeError("WorldCover sampling sensitivity did not pass")
    if terrain.get("status") != "terrain_not_evaluable_before_colour_join":
        raise RuntimeError("terrain access result changed")
    parents = {
        str(args.climate_manifest).replace("\\", "/"): canonical_sha256(args.climate_manifest),
        str(args.worldcover_manifest).replace("\\", "/"): canonical_sha256(args.worldcover_manifest),
        str(args.worldcover_result).replace("\\", "/"): canonical_sha256(args.worldcover_result),
        str(args.terrain_result).replace("\\", "/"): canonical_sha256(args.terrain_result),
    }
    inputs = {}
    for scale in contract["grid"]["scales_km"]:
        for family_file in (
            DATA_ROOT / f"climate_ecoregion_grid_{scale}km.csv",
            DATA_ROOT / f"worldcover_composition_primary_{scale}km.csv",
        ):
            inputs[str(family_file).replace("\\", "/")] = sha256(family_file)
    return {
        "protocol": contract["protocol"],
        "status": "pass_environmental_source_freeze",
        "available_primary_families": ["macroclimate", "land_cover", "ecoregion"],
        "not_evaluable_primary_families": ["terrain"],
        "terrain_decision": terrain["decision"],
        "parent_sha256_lf_canonical_v1": parents,
        "input_grid_sha256": inputs,
        "scales_km": contract["grid"]["scales_km"],
        "coverage_against_final_atlas_opportunity_pending": True,
        "scaleout_colour_opened": False,
        "environment_colour_join_performed": False,
        "claim_ceiling": "Three independent environmental families are frozen before colour; terrain is not evaluable and no flower-colour concordance has been tested.",
    }


def main() -> None:
    args = parse_args()
    result = evaluate(args)
    if args.write:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    elif load(args.manifest) != result:
        raise RuntimeError("committed environmental source freeze changed")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
