#!/usr/bin/env python3
"""Generate the one allowed global location-blind measurement authorization."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_measurement_budget import select_measurement_rows

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"
CANDIDATES = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
BUDGET = ROOT / "docs/supporting/global_monte_carlo_measurement_budget_amendment_v1.json"
OUT = ROOT / "docs/supporting/global_monte_carlo_measurement_authorization_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob(path: Path) -> str:
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{path.relative_to(ROOT).as_posix()}"], text=True).strip()


def taxon_digest(values: pd.Series) -> str:
    text = "\n".join(str(int(x)) for x in sorted(set(values.astype(int)))) + "\n"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing authorization: {OUT}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    budget = json.loads(BUDGET.read_text(encoding="utf-8"))
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
    if budget.get("status") != "frozen_after_metadata_capacity_scale_observed_before_candidate_acquisition_outcome_and_before_any_global_candidate_pixels":
        raise RuntimeError("measurement budget amendment is not frozen pre-pixel")
    if sha256_file(CANDIDATES) != manifest.get("lineage", {}).get("candidate_photos_sha256"):
        raise RuntimeError("candidate photo pool SHA differs from manifest")

    target = int(manifest["capacity_selected_raw_photo_target"])
    full_species = int(manifest["full_target_species"])
    candidate_rows = full_species * target
    frame = pd.read_csv(CANDIDATES)
    if len(frame) != candidate_rows or int(frame["inat_taxon_id"].nunique()) != full_species or int(frame["photo_id"].nunique()) != candidate_rows:
        raise RuntimeError("candidate photo denominator differs from frozen manifest")

    b = budget["measurement_species_budget"]
    maximum_species = int(b["maximum_species"])
    seed = int(b["selection_seed"])
    measured = select_measurement_rows(
        frame,
        target_photos_per_species=target,
        maximum_species=maximum_species,
        seed=seed,
    )
    measurement_species = int(measured["inat_taxon_id"].nunique())
    measurement_rows = int(len(measured))
    if measurement_species != min(full_species, maximum_species):
        raise RuntimeError("measurement species budget did not produce the frozen denominator")
    if measurement_rows != measurement_species * target:
        raise RuntimeError("measurement row denominator drifted")

    auth = {
        "status": "authorize_exactly_one_global_location_blind_measurement_run",
        "branch": BRANCH,
        "candidate_manifest_blob_sha": git_blob(MANIFEST),
        "measurement_budget_amendment_blob_sha": git_blob(BUDGET),
        "candidate_rows": candidate_rows,
        "full_target_species": full_species,
        "measurement_species_budget": maximum_species,
        "measurement_species_selection_seed": seed,
        "measurement_species": measurement_species,
        "measurement_rows": measurement_rows,
        "measurement_taxon_id_sha256": taxon_digest(measured["inat_taxon_id"]),
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
