#!/usr/bin/env python3
"""Freeze passing metadata, dated-source and environmental pre-image evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
RUNNER_PATH = "scripts/data/reconcile_jbi_atlas_dated_source_m2m_v2.py"

from scripts.data.validate_jbi_atlas_preimage_scaleout_evidence import (
    PROTOCOL,
    sha256,
    validate_preimage_evidence,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--dated-source-dir", type=Path, required=True)
    parser.add_argument("--final-coverage-result", type=Path, required=True)
    parser.add_argument("--source-execution-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise RuntimeError("refusing to replace frozen pre-image scale-out evidence")
    output.mkdir(parents=True)
    sources = {
        "scaleout_metadata_manifest.json": args.metadata_dir
        / "scaleout_metadata_manifest.json",
        "scaleout_metadata_feasibility.json": args.metadata_dir
        / "scaleout_metadata_feasibility.json",
        "scaleout_observation_manifest.csv": args.metadata_dir
        / "scaleout_observation_manifest.csv",
        "scaleout_species_panels.csv": args.metadata_dir / "scaleout_species_panels.csv",
        "dated_source_m2m_manifest.json": args.dated_source_dir
        / "dated_source_m2m_manifest.json",
        "dated_source_m2m_reconciliation.json": args.dated_source_dir
        / "dated_source_m2m_reconciliation.json",
        "dated_source_m2m_observation_manifest.csv": args.dated_source_dir
        / "dated_source_m2m_observation_manifest.csv",
        "precolour_environmental_coverage_final.json": args.final_coverage_result,
    }
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        raise RuntimeError(f"pre-image evidence inputs are missing: {missing}")
    for name, source in sources.items():
        shutil.copy2(source, output / name)
    runner_payload = subprocess.run(
        ["git", "show", f"{args.source_execution_commit}:{RUNNER_PATH}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    runner = output / "dated_source_m2m_executable.py"
    runner.write_bytes(runner_payload)
    files = {
        path.name: sha256(path)
        for path in sorted(output.iterdir(), key=lambda item: item.name)
        if path.is_file()
    }
    manifest = {
        "protocol": PROTOCOL,
        "status": "pass_frozen_preimage_scaleout_evidence",
        "metadata_github_actions_run": (
            "https://github.com/zuizui0223/fcp/actions/runs/33405153936"
        ),
        "source_execution": {
            "git_commit": args.source_execution_commit,
            "runner_path": RUNNER_PATH,
            "runner_sha256_exact": sha256(runner),
        },
        "candidate_image_pixels_opened": False,
        "files": files,
        "claim_ceiling": (
            "Metadata, exact dated-source identity and environmental opportunity coverage "
            "only; no candidate image pixel, flower colour or ecological outcome."
        ),
    }
    (output / "preimage_scaleout_evidence_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(validate_preimage_evidence(output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
