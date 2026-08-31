#!/usr/bin/env python3
"""Run the metadata-only feasibility audit for all eight atlas panels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_scaleout import (
    live_api_scaleout_feasibility,
    validate_geometry_admission_amendment,
)
from scripts.data.build_jbi_image_first_atlas_metadata import InaturalistMetadataAdapter


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument(
        "--expansion-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidate-species", type=int, default=500)
    parser.add_argument("--maximum-candidates-per-species", type=int, default=1000)
    parser.add_argument(
        "--geometry-amendment",
        type=Path,
        default=Path(
            "docs/supporting/jbi_atlas_scaleout_geometry_admission_amendment_v1.json"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    atlas = json.loads(args.atlas_contract.read_text(encoding="utf-8"))
    expansion = json.loads(args.expansion_contract.read_text(encoding="utf-8"))
    geometry_amendment = json.loads(args.geometry_amendment.read_text(encoding="utf-8"))
    validate_geometry_admission_amendment(geometry_amendment)
    adapter = InaturalistMetadataAdapter(
        base_url=str(atlas["metadata_source"]["base_url"]),
        pause_seconds=float(atlas["metadata_source"]["request_pause_seconds"]),
    )
    frozen = live_api_scaleout_feasibility(
        atlas,
        expansion,
        geometry_amendment,
        adapter,
        candidate_species_pool_size=args.candidate_species,
        maximum_candidates_per_species=args.maximum_candidates_per_species,
    )
    audit_path = args.output_dir / "scaleout_metadata_feasibility.json"
    write_json(audit_path, frozen.audit)
    files: dict[str, str] = {audit_path.name: file_sha256(audit_path)}
    if frozen.panels:
        panels_path = args.output_dir / "scaleout_species_panels.csv"
        observations_path = args.output_dir / "scaleout_observation_manifest.csv"
        write_csv(panels_path, frozen.panels)
        write_csv(observations_path, frozen.observations)
        files[panels_path.name] = file_sha256(panels_path)
        files[observations_path.name] = file_sha256(observations_path)
    manifest = {
        "protocol": expansion["protocol"],
        "status": frozen.audit["status"],
        "source_role": "live API feasibility; final dated export remains required",
        "candidate_image_pixels_opened": False,
        "geometry_admission_contract": {
            "path": args.geometry_amendment.as_posix(),
            "sha256": file_sha256(args.geometry_amendment),
        },
        "files": files,
    }
    write_json(args.output_dir / "scaleout_metadata_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    if frozen.audit["status"] != "pass_live_api_scaleout_feasibility":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
