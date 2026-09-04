#!/usr/bin/env python3
"""Generate the one allowed global location-blind measurement authorization."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"
CANDIDATES = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
OUT = ROOT / "docs/supporting/global_monte_carlo_measurement_authorization_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob(path: Path) -> str:
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{path.relative_to(ROOT).as_posix()}"], text=True).strip()


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing authorization: {OUT}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete_global_candidate_acquisition_premeasurement_gate_passed":
        raise RuntimeError("candidate acquisition did not pass the frozen premeasurement gate")
    if manifest.get("premeasurement_gate", {}).get("pass") is not True:
        raise RuntimeError("candidate premeasurement gate is not true")
    if manifest.get("candidate_image_pixels_opened") is not False or manifest.get("flower_colour_used") is not False:
        raise RuntimeError("candidate stage already opened forbidden outcomes")
    if manifest.get("measurement_authorized") is not False:
        raise RuntimeError("candidate manifest unexpectedly self-authorized measurement")
    if manifest.get("target_relaxed") is not False or manifest.get("additional_pages_after_result") is not False:
        raise RuntimeError("candidate sampling rules drifted")
    if sha256_file(CANDIDATES) != manifest.get("lineage", {}).get("candidate_photos_sha256"):
        raise RuntimeError("candidate photo pool SHA differs from manifest")

    target = int(manifest["capacity_selected_raw_photo_target"])
    species = int(manifest["full_target_species"])
    expected_rows = species * target
    frame = pd.read_csv(CANDIDATES, usecols=["inat_taxon_id", "photo_id"])
    if len(frame) != expected_rows or int(frame["inat_taxon_id"].nunique()) != species or int(frame["photo_id"].nunique()) != expected_rows:
        raise RuntimeError("candidate photo denominator differs from frozen manifest")

    auth = {
        "status": "authorize_exactly_one_global_location_blind_measurement_run",
        "branch": BRANCH,
        "candidate_manifest_blob_sha": git_blob(MANIFEST),
        "candidate_rows": expected_rows,
        "full_target_species": species,
        "selected_raw_photo_target": target,
        "measurement_rule_change_allowed": False,
        "replacement_after_failure_allowed": False,
        "external_overlay_may_open": False,
    }
    OUT.write_text(json.dumps(auth, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(auth, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
