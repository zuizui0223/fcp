#!/usr/bin/env python3
"""Merge every frozen HTTP429 recovery row back into the immutable v2 audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_capacity_recovery import merge_transport_recovery, nonempty_error_mask

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_capacity_429_recovery_contract_v1.json"
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
FAILED_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"
FAILED_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv"
OUT_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"
OUT_SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv"
OUT_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    for path in (OUT_AUDIT, OUT_SELECTED, OUT_MANIFEST):
        if path.exists():
            raise RuntimeError(f"refusing to overwrite frozen recovery output: {path}")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    failed = json.loads(FAILED_MANIFEST.read_text(encoding="utf-8"))
    if failed.get("status") != contract["trigger"]["required_parent_status"]:
        raise RuntimeError("v2 capacity result does not authorize transport recovery")
    original = pd.read_csv(FAILED_AUDIT)
    if len(original) != int(failed.get("discovered_species_scanned") or -1):
        raise RuntimeError("v2 capacity audit denominator differs from its manifest")

    shard_count = int(contract["execution"]["deterministic_shards"])
    pieces: list[pd.DataFrame] = []
    manifests: list[dict[str, object]] = []
    for shard in range(shard_count):
        csvs = list(args.input_root.rglob(f"capacity_429_recovery_shard_{shard:02d}.csv"))
        js = list(args.input_root.rglob(f"capacity_429_recovery_shard_{shard:02d}.json"))
        if len(csvs) != 1 or len(js) != 1:
            raise RuntimeError(f"recovery shard {shard} requires exactly one CSV and manifest")
        m = json.loads(js[0].read_text(encoding="utf-8"))
        if m.get("status") != "complete_transport_429_recovery_shard":
            raise RuntimeError(f"recovery shard {shard} status invalid")
        if int(m.get("shard_index", -1)) != shard or int(m.get("shard_count", -1)) != shard_count:
            raise RuntimeError(f"recovery shard {shard} identity drift")
        if m.get("candidate_image_pixels_opened") is not False or m.get("flower_colour_used") is not False:
            raise RuntimeError(f"recovery shard {shard} opened forbidden outcomes")
        if sha256_file(csvs[0]) != m["lineage"]["recovery_csv_sha256"]:
            raise RuntimeError(f"recovery shard {shard} CSV SHA mismatch")
        frame = pd.read_csv(csvs[0])
        if len(frame) != int(m["shard_rows"]):
            raise RuntimeError(f"recovery shard {shard} row count mismatch")
        pieces.append(frame)
        manifests.append(m)

    recovered = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
    if recovered["global_row_index"].duplicated().any() if len(recovered) else False:
        raise RuntimeError("duplicate global_row_index across 429 recovery shards")
    merged = merge_transport_recovery(original, recovered)
    if len(merged) != len(original):
        raise RuntimeError("transport recovery changed the global species denominator")
    if merged["global_row_index"].astype(int).tolist() != original["global_row_index"].astype(int).tolist():
        raise RuntimeError("transport recovery changed deterministic global row ordering")
    if merged["inat_taxon_id"].astype(int).tolist() != original["inat_taxon_id"].astype(int).tolist():
        raise RuntimeError("transport recovery changed the taxon denominator")

    targets = [int(x) for x in parent["target_rule"]["candidate_raw_photos_per_species"]]
    target_counts = {str(t): int(merged[f"eligible_raw_{t}"].astype(bool).sum()) for t in targets}
    errors = int(nonempty_error_mask(merged).sum())
    error_fraction = float(errors / len(merged)) if len(merged) else 1.0
    ceiling = float(contract["final_gate"]["maximum_remaining_error_fraction"])
    minimum_species = int(contract["target_rule"]["minimum_metadata_eligible_species"])
    selected_target = next((t for t in targets if target_counts[str(t)] >= minimum_species), None) if error_fraction <= ceiling else None
    if error_fraction > ceiling:
        status = contract["final_gate"]["if_remaining_error_fraction_above_ceiling"]
    elif selected_target is None:
        status = parent["target_rule"]["if_none_pass"]
    else:
        status = contract["final_gate"]["if_pass_and_target_exists"]

    if selected_target is None:
        selected = merged.iloc[0:0][["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"]].copy()
    else:
        selected = merged.loc[
            merged[f"eligible_raw_{selected_target}"].astype(bool),
            ["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"],
        ].copy().sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
        selected["selected_raw_photo_target"] = int(selected_target)

    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_AUDIT, index=False, lineterminator="\n")
    selected.to_csv(OUT_SELECTED, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": str(status),
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "original_v2_status": failed["status"],
        "original_v2_discovered_species_scanned": int(failed["discovered_species_scanned"]),
        "discovered_species_scanned": int(len(merged)),
        "original_v2_request_errors": int(failed["request_errors"]),
        "recovered_429_rows": int(len(recovered)),
        "remaining_request_errors": errors,
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": ceiling,
        "target_counts": target_counts,
        "selected_raw_photo_target": None if selected_target is None else int(selected_target),
        "selected_species": int(len(selected)),
        "minimum_metadata_eligible_species": minimum_species,
        "actual_image_acquisition_authorized": False,
        "second_recovery_permitted": False,
        "biological_rules_changed": False,
        "lineage": {
            "recovery_contract_sha256": sha256_file(CONTRACT),
            "failed_capacity_manifest_sha256": sha256_file(FAILED_MANIFEST),
            "failed_capacity_audit_sha256": sha256_file(FAILED_AUDIT),
            "recovery_csv_sha256": [m["lineage"]["recovery_csv_sha256"] for m in manifests],
            "audit_sha256": sha256_file(OUT_AUDIT),
            "selected_sha256": sha256_file(OUT_SELECTED)
        },
        "files": {
            "species_audit": str(OUT_AUDIT.relative_to(ROOT)),
            "selected_species_frame": str(OUT_SELECTED.relative_to(ROOT))
        }
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
