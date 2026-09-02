#!/usr/bin/env python3
"""Evaluate final pre-colour environmental coverage after UUID dated-source v4.

This runner reuses the already-frozen environmental boundary surfaces and terminal
geometry.  It does not use image pixels or colour and cannot authorize acquisition
by itself.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import evaluate_environmental_coverage_gate
from scripts.data.run_jbi_atlas_precolour_environmental_coverage import (
    file_sha256,
    load_frozen_boundaries,
    read_csv,
    selected_geometry,
)


DATED_PROTOCOL = "jbi-atlas-dated-source-uuid-bucket-amendment-v4"
DATED_PASS = "pass_dated_source_uuid_bucket_scaleout_freeze"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def validate_v4_dated_source(dated: Mapping[str, Any]) -> None:
    if (
        dated.get("protocol") != DATED_PROTOCOL
        or dated.get("status") != DATED_PASS
        or dated.get("candidate_image_pixels_opened") is not False
        or dated.get("continuous_colour_used") is not False
        or dated.get("selected_species") != 200
        or dated.get("selected_photo_assets") != 60000
        or dated.get("frozen_observations") != 60000
        or dated.get("replacement_permitted") is not False
        or dated.get("image_acquisition_authorized") is not False
    ):
        raise ValueError("UUID dated-source v4 did not pass unchanged")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-feasibility", type=Path, required=True)
    parser.add_argument("--species-panels", type=Path, required=True)
    parser.add_argument("--dated-source-reconciliation", type=Path, required=True)
    parser.add_argument(
        "--environment-dir", type=Path, default=ROOT / "data/atlas/environment"
    )
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json",
    )
    parser.add_argument(
        "--boundary-manifest",
        type=Path,
        default=ROOT / "data/atlas/environment/environmental_boundary_freeze_manifest.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    feasibility = load_json(args.metadata_feasibility)
    panels = read_csv(args.species_panels)
    dated = load_json(args.dated_source_reconciliation)
    validate_v4_dated_source(dated)
    contract = load_json(args.environment_contract)
    boundary_manifest = load_json(args.boundary_manifest)
    rows_by_scale, boundary_hashes = load_frozen_boundaries(
        args.environment_dir, boundary_manifest
    )
    taxa, geometry = selected_geometry(feasibility, panels)
    gate = evaluate_environmental_coverage_gate(
        geometry,
        taxa,
        rows_by_scale,
        contract,
    )
    result = {
        **gate,
        "status": gate["status"],
        "coverage_gate_status": gate["status"],
        "source_stage": "final-dated-source",
        "dated_source_protocol": DATED_PROTOCOL,
        "final_dated_source_required": False,
        "image_acquisition_authorized": False,
        "selected_taxon_ids": taxa,
        "claim_ceiling": (
            "Final pre-colour opportunity-cell coverage only; no image pixels, flower colour, "
            "environment-colour concordance, or biological result."
        ),
        "parents": {
            "metadata_feasibility_sha256": file_sha256(args.metadata_feasibility),
            "species_panels_sha256": file_sha256(args.species_panels),
            "dated_source_reconciliation_sha256": file_sha256(
                args.dated_source_reconciliation
            ),
            "environment_contract_sha256": file_sha256(args.environment_contract),
            "boundary_manifest_sha256": file_sha256(args.boundary_manifest),
            "boundary_files_sha256": boundary_hashes,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if gate["status"] == "pass_precolour_environmental_coverage" else 2


if __name__ == "__main__":
    raise SystemExit(main())
