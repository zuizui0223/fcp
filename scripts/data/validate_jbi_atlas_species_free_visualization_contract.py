#!/usr/bin/env python3
"""Validate the pre-pixel species-free atlas display contract and parents."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/jbi_atlas_species_free_visualization_contract_v1.json"


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def validate_visualization_contract(value: dict[str, object]) -> None:
    if (
        value.get("protocol") != "jbi-atlas-species-free-visualization-v1"
        or value.get("status")
        != "prospectively_frozen_before_any_scaleout_candidate_pixel_or_colour"
        or any(flag is not False for flag in value.get("outcome_firewall", {}).values())
    ):
        raise RuntimeError("species-free visualization was not frozen before outcomes")
    map_rule = value.get("map", {})
    photo = value.get("photo_bar", {})
    if (
        map_rule.get("projection") != "Mollweide"
        or map_rule.get("species_labels_visible") is not False
        or photo.get("count") != 48
        or photo.get("selection_fields")
        != ["longitude", "latitude", "measurement_id"]
        or any(
            photo.get(field) is not False
            for field in (
                "colour_used_for_selection",
                "species_or_cohort_used_for_selection",
                "environmental_or_inference_outcome_used_for_selection",
            )
        )
        or photo.get("species_labels_visible") is not False
        or "without replacement" not in str(photo.get("failure_rule", ""))
    ):
        raise RuntimeError("species-free visualization rule changed")
    for parent in value.get("parents", {}).values():
        path = ROOT / parent["path"]
        if canonical_sha256(path) != parent["sha256_lf_canonical_v1"]:
            raise RuntimeError(f"visualization parent changed: {path}")


def main() -> None:
    value = json.loads(CONTRACT.read_text(encoding="utf-8"))
    validate_visualization_contract(value)
    print(
        json.dumps(
            {
                "status": "pass_species_free_visualization_contract",
                "photo_bar_count": value["photo_bar"]["count"],
                "candidate_image_pixels_opened": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
