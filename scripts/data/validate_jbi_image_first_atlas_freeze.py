#!/usr/bin/env python3
"""Validate the atlas contract or the complete committed metadata/geometry freeze."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.image_first_atlas import validate_atlas_contract


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument("--contract-only", action="store_true")
    parser.add_argument(
        "--cohort",
        type=Path,
        default=Path("data/atlas/jbi_image_first_atlas_cohort_v1.csv"),
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path("data/atlas/jbi_image_first_atlas_observation_manifest_v1.csv"),
    )
    parser.add_argument(
        "--feasibility",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_metadata_feasibility_v1.json"),
    )
    parser.add_argument(
        "--geometry",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_geometry_freeze_v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_freeze_manifest_v1.json"),
    )
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_atlas_contract(contract)
    if args.contract_only:
        print("Atlas contract is frozen before metadata queries and image pixels.")
        return 0

    required = [args.cohort, args.observations, args.feasibility, args.geometry, args.manifest]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit(f"partial atlas freeze; missing: {missing}")

    cohort = read_csv(args.cohort)
    observations = read_csv(args.observations)
    feasibility = json.loads(args.feasibility.read_text(encoding="utf-8"))
    geometry = json.loads(args.geometry.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if len(cohort) != 50:
        raise SystemExit(f"atlas cohort must contain 50 species, found {len(cohort)}")
    if len({row["species"] for row in cohort}) != 50:
        raise SystemExit("atlas cohort species are not unique")

    counts = Counter(row["species"] for row in observations)
    if set(counts) != {row["species"] for row in cohort}:
        raise SystemExit("observation and cohort species differ")
    if any(count not in {300, 400, 500} for count in counts.values()):
        raise SystemExit(f"unexpected per-species sample sizes: {sorted(set(counts.values()))}")
    if any(float(row["positional_accuracy_m"]) > 5000 for row in observations):
        raise SystemExit("observation exceeds the frozen 5-km accuracy ceiling")
    if len({row["photo_id"] for row in observations}) != len(observations):
        raise SystemExit("photo IDs are not globally unique")
    if len({row["observation_id"] for row in observations}) != len(observations):
        raise SystemExit("observation IDs are not globally unique")

    if feasibility.get("status") != "pass_50_species_metadata_only":
        raise SystemExit("metadata feasibility did not pass all 50 species")
    if geometry.get("status") != "geometry_scale_frozen":
        raise SystemExit("geometry-only scale selection did not pass")
    if geometry.get("selected_primary_scale_km") not in {100, 250, 500}:
        raise SystemExit("geometry primary scale is outside the frozen candidates")
    if manifest.get("status") != "metadata_and_geometry_frozen_before_image_pixels":
        raise SystemExit("freeze manifest is not at the pre-image gate")
    for key in (
        "candidate_image_pixels_opened",
        "flower_roi_used",
        "continuous_colour_used",
        "literature_classification_used_for_admission",
    ):
        if manifest.get(key) is not False:
            raise SystemExit(f"freeze manifest outcome firewall is open for {key}")

    expected_files = {
        args.contract: manifest["files"].get(str(args.contract).replace("\\", "/")),
        args.cohort: manifest["files"].get(str(args.cohort).replace("\\", "/")),
        args.observations: manifest["files"].get(str(args.observations).replace("\\", "/")),
        args.feasibility: manifest["files"].get(str(args.feasibility).replace("\\", "/")),
        args.geometry: manifest["files"].get(str(args.geometry).replace("\\", "/")),
    }
    for path, expected in expected_files.items():
        if expected != sha256(path):
            raise SystemExit(f"freeze hash mismatch: {path}")

    print(
        f"Validated image-first atlas freeze: 50 species, {len(observations)} observations, "
        f"primary scale {geometry['selected_primary_scale_km']} km."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
