#!/usr/bin/env python3
"""Generate the one allowed transport-only HTTP429 capacity recovery authorization."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_capacity_recovery import frozen_recovery_rows

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"
AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_capacity_429_recovery_contract_v1.json"
OUT = ROOT / "docs/supporting/global_monte_carlo_capacity_429_recovery_authorization_v1.json"


def git_blob(path: Path) -> str:
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{path.relative_to(ROOT).as_posix()}"], text=True).strip()


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing authorization: {OUT}")
    failed = json.loads(MANIFEST.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    trigger = contract["trigger"]
    if failed.get("status") != trigger["required_parent_status"]:
        raise RuntimeError("capacity v2 is not formally not-evaluable from request failure")
    if float(failed.get("request_error_fraction") or 0.0) <= float(trigger["required_parent_request_error_fraction_above"]):
        raise RuntimeError("capacity v2 did not exceed the frozen request-error ceiling")
    if failed.get("candidate_image_pixels_opened") is not False or failed.get("flower_colour_used") is not False:
        raise RuntimeError("capacity v2 opened forbidden outcomes")
    audit = pd.read_csv(AUDIT)
    if len(audit) != int(failed.get("discovered_species_scanned") or -1):
        raise RuntimeError("capacity v2 audit denominator differs from manifest")
    retry = frozen_recovery_rows(audit)
    if len(retry) == 0:
        raise RuntimeError("capacity v2 contains no exact HTTP429 rows eligible for recovery")

    auth = {
        "status": "authorize_exactly_one_transport_429_capacity_recovery",
        "branch": BRANCH,
        "failed_capacity_manifest_blob_sha": git_blob(MANIFEST),
        "retry_rows": int(len(retry)),
        "candidate_image_pixels_may_open": False,
        "flower_colour_may_open": False,
        "successful_rows_may_be_requeried": False,
        "non_429_rows_may_be_requeried": False,
        "second_recovery_allowed": False,
    }
    OUT.write_text(json.dumps(auth, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(auth, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
