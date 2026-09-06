#!/usr/bin/env python3
"""Generate the one allowed global RGFCA G1/G3 inference authorization."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
OUT = ROOT / "docs/supporting/global_monte_carlo_inference_authorization_v1.json"


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
    measurement = json.loads(MANIFEST.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    gate = contract["input_gate"]
    if measurement.get("status") != gate["measurement_status_required"]:
        raise RuntimeError("measurement status does not permit frozen inference")
    post = measurement.get("postmeasurement_gate", {})
    if post.get("pass") is not True or post.get("decision") != gate["measurement_postgate_decision_required"]:
        raise RuntimeError("measurement postgate did not permit frozen G1/G3 inference")
    if measurement.get("g1_g3_inference_run") is not False or measurement.get("external_overlay_opened") is not False:
        raise RuntimeError("measurement result indicates inference/overlay already opened")
    if sha256_file(MEASURED) != measurement.get("lineage", {}).get("measured_table_sha256"):
        raise RuntimeError("measured table SHA differs from measurement manifest")

    minimum = int(gate["minimum_classifiable_photos_per_species"])
    frame = pd.read_csv(MEASURED, usecols=["species", "global_classifiable"])
    flag = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    counts = frame.loc[flag].groupby("species", observed=True).size()
    evaluable = int((counts >= minimum).sum())
    if evaluable != int(post.get("evaluable_species")):
        raise RuntimeError("measured table and measurement manifest disagree on evaluable species")
    if evaluable < int(gate["minimum_inferential_species"]):
        raise RuntimeError("measurement pool is below the frozen inferential species minimum")

    auth = {
        "status": "authorize_exactly_one_global_rgfca_g1_g3_primary_inference",
        "branch": BRANCH,
        "measurement_manifest_blob_sha": git_blob(MANIFEST),
        "evaluable_species": evaluable,
        "g1_g3_rule_change_allowed": False,
        "null_rule_change_allowed": False,
        "species_or_photo_replacement_allowed": False,
        "external_overlay_may_open_before_g1": False,
    }
    OUT.write_text(json.dumps(auth, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(auth, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
