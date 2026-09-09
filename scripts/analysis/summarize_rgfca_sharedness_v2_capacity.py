#!/usr/bin/env python3
"""Summarize frozen metadata capacity for a high-photo-count RGFCA sharedness v2 design."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_capacity_extension_contract_v1.json"
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json"
AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def summarize() -> dict[str, object]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    if contract["status"] != "frozen_before_sharedness_v2_capacity_summary_and_before_any_new_capacity_query":
        raise RuntimeError("sharedness v2 capacity contract drifted")
    if parent["status"] != contract["parent_capacity"]["required_status"]:
        raise RuntimeError("parent capacity status drifted")
    if parent["candidate_image_pixels_opened"] is not False or parent["flower_colour_used"] is not False:
        raise RuntimeError("parent capacity unexpectedly opened biological outcomes")
    frame = pd.read_csv(AUDIT)
    required = {"inat_taxon_id", "raw_results", "after_observer_cap", "request_error"}
    if not required.issubset(frame.columns):
        raise RuntimeError("capacity audit lacks required columns")
    if len(frame) != int(contract["parent_capacity"]["species_census"]):
        raise RuntimeError("species census differs from frozen contract")
    if frame["inat_taxon_id"].nunique() != len(frame):
        raise RuntimeError("duplicate species/taxon rows in capacity audit")
    if frame["request_error"].fillna("").astype(str).str.len().gt(0).any():
        raise RuntimeError("parent v3 audit still contains request errors")
    raw = pd.to_numeric(frame["raw_results"], errors="raise").astype(int)
    capped = pd.to_numeric(frame["after_observer_cap"], errors="raise").astype(int)
    ceiling = int(contract["parent_capacity"]["old_query_per_page_ceiling"])
    if int(raw.max()) != ceiling or (raw > ceiling).any():
        raise RuntimeError("parent request ceiling assumption is false")
    thresholds = [int(x) for x in contract["fixed_existing-census_thresholds_after_observer_cap"]]
    censored = raw.eq(ceiling)
    counts = {str(t): int(capped.ge(t).sum()) for t in thresholds}
    censored_counts = {str(t): int((censored & capped.ge(t)).sum()) for t in thresholds}
    result = {
        "protocol": contract["protocol"],
        "status": "complete_metadata_only_existing_census_capacity_summary",
        "species_census": int(len(frame)),
        "request_errors": 0,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "old_query_per_page_ceiling": ceiling,
        "right_censored_species_raw_results_equal_ceiling": int(censored.sum()),
        "maximum_observed_raw_results": int(raw.max()),
        "maximum_observed_after_observer_cap": int(capped.max()),
        "after_observer_cap_threshold_counts": counts,
        "right_censored_species_after_observer_cap_threshold_counts": censored_counts,
        "availability_above_200_raw_results_identifiable_from_parent_census": False,
        "can_infer_300_or_400_photo_capacity_from_parent_census": False,
        "next_required_step": "prospective_metadata_only_extension_on_right_censored_species_before_any_300_or_400_photo_design_claim",
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "parent_manifest_sha256": sha256_file(PARENT),
            "parent_species_audit_sha256": sha256_file(AUDIT),
        },
        "claim_ceiling": "metadata capacity only; not classifiable yield, sharedness power, flower-specific ecology, or an acquisition authorization",
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
