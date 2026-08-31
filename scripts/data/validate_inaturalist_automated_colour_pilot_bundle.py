#!/usr/bin/env python3
"""Fail-closed validation for the publishable automated-colour pilot bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SUPPORTING = ROOT / "docs" / "supporting"
PROTOCOL = "fcp-inaturalist-automated-colour-state-v2"
ADMITTED = "automated_colour_state_admitted"
EXPECTED = {
    "Erythranthe lewisii": {
        "admitted": 101,
        "rho": 0.007735753532729214,
        "q": 0.59685,
        "contrast": 0.04356050348766461,
        "contrast_q": 0.5105999999999999,
    },
    "Hesperis matronalis": {
        "admitted": 97,
        "rho": -0.0582161511391038,
        "q": 0.9636999999999999,
        "contrast": -0.10022023377470603,
        "contrast_q": 0.9664,
    },
    "Orchis mascula": {
        "admitted": 108,
        "rho": 0.01995986362441413,
        "q": 0.59685,
        "contrast": -0.08759430241419748,
        "contrast_q": 0.9664,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def close(first: float, second: float) -> bool:
    return math.isclose(first, second, rel_tol=0.0, abs_tol=1e-14)


def validate() -> dict[str, Any]:
    plan_path = SUPPORTING / "jbi_inaturalist_automated_colour_plan_manifest_v2.json"
    development_path = SUPPORTING / "jbi_inaturalist_automated_colour_development_gate_v2.json"
    cache_path = SUPPORTING / "jbi_inaturalist_automated_colour_locked_cache_audit_v2.json"
    extraction_path = SUPPORTING / "jbi_inaturalist_automated_colour_locked_extraction_manifest_v2.json"
    locked_path = SUPPORTING / "jbi_inaturalist_automated_colour_locked_gate_v2.json"
    table_path = SUPPORTING / "jbi_inaturalist_automated_colour_locked_analysis_input_v2.csv"
    table_manifest_path = SUPPORTING / "jbi_inaturalist_automated_colour_locked_analysis_input_manifest_v2.json"
    spatial_manifest_path = SUPPORTING / "jbi_inaturalist_automated_colour_spatial_manifest_v2.json"
    result_path = SUPPORTING / "jbi_inaturalist_automated_colour_spatial_v2" / "species_spatial_mark_results.csv"
    variogram_path = SUPPORTING / "jbi_inaturalist_automated_colour_spatial_v2" / "descriptive_equal_pair_variogram.csv"
    figure_manifest_path = SUPPORTING / "jbi_inaturalist_automated_colour_pilot_figure_manifest_v2.json"

    plan = load_json(plan_path)
    development = load_json(development_path)
    cache = load_json(cache_path)
    extraction = load_json(extraction_path)
    locked = load_json(locked_path)
    table_manifest = load_json(table_manifest_path)
    spatial = load_json(spatial_manifest_path)
    figure_manifest = load_json(figure_manifest_path)

    require(plan["protocol"] == PROTOCOL, "plan protocol mismatch")
    require(development["protocol"] == PROTOCOL, "development protocol mismatch")
    require(development["species_passed"] == 3 and development["species_not_evaluable"] == 3, "development decision mismatch")
    require(cache["status"] == "complete_validated_cache", "cache audit is not complete")
    require(cache["valid_cache_records"] == cache["expected_photos"] == 717, "cache count mismatch")
    require(not cache["missing_cache_records"] and not cache["partial_files"] and not cache["unexpected_cache_records"], "cache integrity failure")
    require(extraction["selected_encounters"] == 360 and extraction["selected_photos"] == 717, "locked extraction count mismatch")
    require(locked["species_passing_locked_gate"] == 3 and not locked["coordinate_withheld_species"], "coordinate firewall decision mismatch")
    require(locked["locked_spatial_input_sha256"] == table_manifest["source_private_sha256"], "private-to-public provenance mismatch")

    table_rows = read_csv(table_path)
    require(len(table_rows) == 360, "public table must contain 360 rows")
    require(sha256(table_path) == table_manifest["public_table_sha256"], "public table hash mismatch")
    admitted_rows = [row for row in table_rows if row["encounter_status"] == ADMITTED]
    require(len(admitted_rows) == table_manifest["admitted_rows"] == 306, "public admitted count mismatch")
    require(len({row["encounter_blind_id"] for row in table_rows}) == 360, "duplicate blind encounter")
    require(all(re.fullmatch(r"observer_\d{4}_[0-9a-f]{8}", row["observer_id"]) for row in admitted_rows), "observer pseudonym format mismatch")
    require(all(row["latitude"] and row["longitude"] for row in admitted_rows), "admitted coordinate missing")
    rejected_rows = [row for row in table_rows if row["encounter_status"] != ADMITTED]
    require(all(not row[field] for row in rejected_rows for field in ("latitude", "longitude", "observer_id")), "non-admitted privacy fields were not closed")

    require(spatial["status"] == "complete_locked_species_conditioned_random_mark_test", "spatial run incomplete")
    require(spatial["permutations"] == 9999 and spatial["fdr_q"] == 0.05, "spatial contract mismatch")
    require(not spatial["universality_claim_allowed"] and not spatial["mechanism_claim_allowed"], "claim ceiling opened")
    require(spatial["source_sha256"]["development_gate"] == sha256(development_path), "development source hash mismatch")
    require(spatial["source_sha256"]["locked_table"] == sha256(table_path), "public table source hash mismatch")
    runner_path = ROOT / "scripts" / "data" / "run_inaturalist_automated_colour_spatial_marks.py"
    require(spatial["source_sha256"]["runner"] == sha256(runner_path), "runner source hash mismatch")
    require(spatial["output_sha256"][result_path.name] == sha256(result_path), "result table hash mismatch")
    require(spatial["output_sha256"][variogram_path.name] == sha256(variogram_path), "variogram hash mismatch")

    result_rows = read_csv(result_path)
    require({row["canonical_name"] for row in result_rows} == set(EXPECTED), "result species mismatch")
    for row in result_rows:
        expected = EXPECTED[row["canonical_name"]]
        require(int(row["admitted_encounters"]) == expected["admitted"], "admitted result mismatch")
        require(close(float(row["primary_rho"]), expected["rho"]), "primary rho changed")
        require(close(float(row["primary_bh_q"]), expected["q"]), "primary q changed")
        require(close(float(row["flower_minus_background_rho"]), expected["contrast"]), "contrast changed")
        require(close(float(row["contrast_bh_q"]), expected["contrast_q"]), "contrast q changed")
        require(row["spatial_claim_status"] == "spatial_organization_not_detected", "negative decision changed")

    figure_script = ROOT / "scripts" / "analysis" / "make_jbi_inaturalist_automated_colour_pilot_figures.py"
    require(figure_manifest["display_species_labels_removed"], "figure display is not species-free")
    require(figure_manifest["inference_species_conditioned"], "figure manifest lost species conditioning")
    require(figure_manifest["source_sha256"]["analysis_table"] == sha256(table_path), "figure table hash mismatch")
    require(figure_manifest["source_sha256"]["development_gate"] == sha256(development_path), "figure development hash mismatch")
    require(figure_manifest["source_sha256"]["spatial_results"] == sha256(result_path), "figure result hash mismatch")
    require(figure_manifest["source_sha256"]["figure_script"] == sha256(figure_script), "figure script hash mismatch")
    for name, expected_hash in figure_manifest["output_sha256"].items():
        require(sha256(ROOT / "docs" / "figures" / name) == expected_hash, f"figure hash mismatch: {name}")
    require(len(figure_manifest["displayed_photo_audit"]) == 24, "photo bar must contain 24 audited crops")
    require(all(row["photo_license"] in {"cc0", "cc-by"} for row in figure_manifest["displayed_photo_audit"]), "non-reusable displayed photo")

    manuscript = (ROOT / "docs" / "JBI_INATURALIST_AUTOMATED_COLOUR_PILOT_MANUSCRIPT.md").read_text(encoding="utf-8")
    status = (ROOT / "docs" / "JBI_IMAGE_FIRST_ATLAS_STATUS.md").read_text(encoding="utf-8")
    require("did not detect" in manuscript and "do not prove spatial randomness" in manuscript, "manuscript claim ceiling missing")
    require("20,200 candidate images remain unopened" in manuscript, "atlas image firewall missing from manuscript")
    require("bulk atlas image opening is stopped" in status, "status does not record STOP")

    return {
        "status": "pass_publishable_negative_validation_bundle",
        "protocol": PROTOCOL,
        "development_species_passed": 3,
        "locked_species_tested": 3,
        "spatial_organization_supported": 0,
        "public_analysis_rows": 360,
        "displayed_photo_crops": 24,
        "atlas_candidate_pixels_opened": False,
        "claim_ceiling": "Three-species model-consensus flower-candidate non-detection; no randomness, tissue, morph, mechanism or universality claim.",
    }


def main() -> None:
    print(json.dumps(validate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
