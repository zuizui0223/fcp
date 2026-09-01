#!/usr/bin/env python3
"""Validate the committed 200-species pre-image scale-out evidence bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = "jbi-atlas-preimage-scaleout-evidence-v1"
RUNNER_PATH = "scripts/data/reconcile_jbi_atlas_dated_source_m2m_v2.py"
EXPECTED_FILES = {
    "scaleout_metadata_manifest.json",
    "scaleout_metadata_feasibility.json",
    "scaleout_observation_manifest.csv",
    "scaleout_species_panels.csv",
    "dated_source_m2m_manifest.json",
    "dated_source_m2m_reconciliation.json",
    "dated_source_m2m_observation_manifest.csv",
    "precolour_environmental_coverage_final.json",
    "dated_source_m2m_executable.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _false(value: object) -> bool:
    return value is False or str(value).strip().casefold() == "false"


def validate_preimage_evidence(evidence_dir: Path) -> dict[str, Any]:
    evidence_dir = evidence_dir.resolve()
    manifest = json.loads(
        (evidence_dir / "preimage_scaleout_evidence_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        manifest.get("protocol") != PROTOCOL
        or manifest.get("status") != "pass_frozen_preimage_scaleout_evidence"
        or manifest.get("candidate_image_pixels_opened") is not False
        or manifest.get("metadata_github_actions_run")
        != "https://github.com/zuizui0223/fcp/actions/runs/33405153936"
        or set(manifest.get("files", {})) != EXPECTED_FILES
    ):
        raise RuntimeError("pre-image scale-out evidence manifest changed")
    for name, expected in manifest["files"].items():
        path = evidence_dir / name
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f"pre-image scale-out artifact changed: {name}")

    metadata_manifest = json.loads(
        (evidence_dir / "scaleout_metadata_manifest.json").read_text(encoding="utf-8")
    )
    feasibility_path = evidence_dir / "scaleout_metadata_feasibility.json"
    observations_path = evidence_dir / "scaleout_observation_manifest.csv"
    panels_path = evidence_dir / "scaleout_species_panels.csv"
    if (
        metadata_manifest.get("status") != "pass_live_api_scaleout_feasibility"
        or metadata_manifest.get("candidate_image_pixels_opened") is not False
        or metadata_manifest.get("files", {}).get(feasibility_path.name)
        != sha256(feasibility_path)
        or metadata_manifest.get("files", {}).get(observations_path.name)
        != sha256(observations_path)
        or metadata_manifest.get("files", {}).get(panels_path.name)
        != sha256(panels_path)
    ):
        raise RuntimeError("metadata feasibility evidence changed")
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (
        feasibility.get("status") != "pass_live_api_scaleout_feasibility"
        or feasibility.get("frozen_species") != 200
        or feasibility.get("frozen_observations") != 60000
        or feasibility.get("cohorts") != 8
        or not _false(feasibility.get("candidate_image_pixels_opened"))
        or not _false(feasibility.get("continuous_colour_used"))
    ):
        raise RuntimeError("metadata feasibility denominator changed")

    panels = read_csv(panels_path)
    observations = read_csv(observations_path)
    cohort_counts = Counter(row["cohort_id"] for row in panels)
    if (
        len(panels) != 200
        or len({row["taxon_id"] for row in panels}) != 200
        or len({row["genus"] for row in panels}) != 200
        or cohort_counts != Counter({f"C{index:02d}": 25 for index in range(1, 9)})
        or len(observations) != 60000
        or len({row["observation_id"] for row in observations}) != 60000
        or len({row["photo_id"] for row in observations}) != 60000
        or any(not _false(row["candidate_image_pixels_opened"]) for row in observations)
    ):
        raise RuntimeError("metadata panel or observation denominator changed")

    source_manifest = json.loads(
        (evidence_dir / "dated_source_m2m_manifest.json").read_text(encoding="utf-8")
    )
    reconciliation_path = evidence_dir / "dated_source_m2m_reconciliation.json"
    dated_path = evidence_dir / "dated_source_m2m_observation_manifest.csv"
    if (
        source_manifest.get("status") != "pass_dated_source_m2m_scaleout_freeze"
        or source_manifest.get("candidate_image_pixels_opened") is not False
        or source_manifest.get("files", {}).get(reconciliation_path.name)
        != sha256(reconciliation_path)
        or source_manifest.get("files", {}).get(dated_path.name) != sha256(dated_path)
    ):
        raise RuntimeError("dated-source manifest changed")
    reconciliation = json.loads(reconciliation_path.read_text(encoding="utf-8"))
    if (
        reconciliation.get("status") != "pass_dated_source_m2m_scaleout_freeze"
        or reconciliation.get("selected_species") != 200
        or reconciliation.get("selected_photo_assets") != 60000
        or reconciliation.get("frozen_observations") != 60000
        or reconciliation.get("association_resolution_failure_count") != 0
        or reconciliation.get("asset_field_conflict_count") != 0
        or any(reconciliation.get("missing_counts", {}).values())
        or reconciliation.get("candidate_image_pixels_opened") is not False
        or reconciliation.get("continuous_colour_used") is not False
        or reconciliation.get("replacement_permitted") is not False
        or reconciliation.get("image_acquisition_authorized") is not False
        or reconciliation.get("parents", {}).get("metadata_feasibility_sha256")
        != sha256(feasibility_path)
        or reconciliation.get("parents", {}).get("species_panels_sha256")
        != sha256(panels_path)
        or reconciliation.get("parents", {}).get("selected_observations_sha256")
        != sha256(observations_path)
    ):
        raise RuntimeError("dated-source reconciliation changed")
    dated = read_csv(dated_path)
    if (
        len(dated) != 60000
        or len({row["observation_id"] for row in dated}) != 60000
        or len({row["photo_id"] for row in dated}) != 60000
        or any(not _false(row["candidate_image_pixels_opened"]) for row in dated)
    ):
        raise RuntimeError("dated-source observation denominator changed")

    coverage_path = evidence_dir / "precolour_environmental_coverage_final.json"
    coverage = json.loads(coverage_path.read_text(encoding="utf-8"))
    if (
        coverage.get("status") != "pass_precolour_environmental_coverage"
        or coverage.get("coverage_gate_status")
        != "pass_precolour_environmental_coverage"
        or coverage.get("source_stage") != "final-dated-source"
        or coverage.get("final_dated_source_required") is not False
        or coverage.get("scaleout_colour_opened") is not False
        or coverage.get("image_acquisition_authorized") is not False
        or coverage.get("parents", {}).get("metadata_feasibility_sha256")
        != sha256(feasibility_path)
        or coverage.get("parents", {}).get("species_panels_sha256")
        != sha256(panels_path)
        or coverage.get("parents", {}).get("dated_source_reconciliation_sha256")
        != sha256(reconciliation_path)
    ):
        raise RuntimeError("final pre-colour environmental coverage changed")

    execution = manifest.get("source_execution", {})
    runner = evidence_dir / "dated_source_m2m_executable.py"
    try:
        committed = subprocess.run(
            ["git", "show", f"{execution['git_commit']}:{RUNNER_PATH}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (KeyError, subprocess.CalledProcessError) as exc:
        raise RuntimeError("dated-source execution commit is unavailable") from exc
    if (
        hashlib.sha256(committed).hexdigest() != execution.get("runner_sha256_exact")
        or sha256(runner) != execution.get("runner_sha256_exact")
    ):
        raise RuntimeError("dated-source execution program changed")
    return {
        "status": "pass_committed_preimage_scaleout_evidence",
        "species": 200,
        "observations": 60000,
        "candidate_image_pixels_opened": False,
        "next_gate": "build the measurement firewall",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate_preimage_evidence(args.evidence_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
