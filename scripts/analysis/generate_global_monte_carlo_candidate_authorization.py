#!/usr/bin/env python3
"""Generate the exact metadata-only candidate acquisition authorization from a frozen capacity result.

This helper is intentionally fail-closed. It never chooses a target or species set itself;
it only copies the already-frozen values from the completed capacity manifest after
verifying the selected-species frame and Git blob lineage.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
CAPACITY = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"
SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv"
OUT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_authorization_v1.json"


def git_blob(path: Path) -> str:
    rel = path.relative_to(ROOT).as_posix()
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{rel}"], text=True).strip()


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing authorization: {OUT}")
    if not CAPACITY.exists() or not SELECTED.exists():
        raise RuntimeError("frozen capacity result is incomplete")

    capacity = json.loads(CAPACITY.read_text(encoding="utf-8"))
    if capacity.get("status") != "complete_metadata_only_capacity_scan_target_selected_v2":
        raise RuntimeError("capacity result did not pass the frozen request-error and target-selection gates")
    if capacity.get("request_coverage_ok") is not True:
        raise RuntimeError("capacity request coverage did not pass")
    if capacity.get("candidate_image_pixels_opened") is not False or capacity.get("flower_colour_used") is not False:
        raise RuntimeError("capacity stage opened forbidden outcomes")
    if capacity.get("actual_image_acquisition_authorized") is not False:
        raise RuntimeError("capacity manifest unexpectedly self-authorized acquisition")

    target = int(capacity.get("selected_raw_photo_target"))
    if target not in (100, 80, 60):
        raise RuntimeError("capacity target is outside the prospectively frozen set")
    selected_n = int(capacity.get("selected_species") or 0)
    if selected_n < int(capacity.get("minimum_metadata_eligible_species") or 300):
        raise RuntimeError("capacity selected fewer than the frozen minimum species")

    selected = pd.read_csv(SELECTED)
    if len(selected) != selected_n or selected["inat_taxon_id"].nunique() != selected_n:
        raise RuntimeError("selected-species frame does not match capacity manifest")
    if "selected_raw_photo_target" not in selected.columns:
        raise RuntimeError("selected-species frame lacks selected_raw_photo_target")
    if not (selected["selected_raw_photo_target"].astype(int) == target).all():
        raise RuntimeError("selected-species frame target differs from capacity manifest")

    authorization = {
        "status": "authorize_exactly_one_metadata_only_global_candidate_acquisition",
        "branch": BRANCH,
        "capacity_manifest_blob_sha": git_blob(CAPACITY),
        "selected_raw_photo_target": target,
        "selected_species": selected_n,
        "candidate_image_pixels_may_open": False,
        "flower_colour_may_open": False,
        "target_relaxation_allowed": False,
        "additional_pages_after_result_allowed": False,
        "rerun_for_favourable_species_set_allowed": False,
    }
    OUT.write_text(json.dumps(authorization, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(authorization, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
