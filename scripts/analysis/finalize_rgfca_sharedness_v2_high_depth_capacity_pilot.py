#!/usr/bin/env python3
"""Strictly combine all shards of the frozen high-depth RGFCA capacity pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_high_depth_capacity_pilot_contract_v1.json"
PARENT_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    args = p.parse_args()
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    shards = int(c["execution"]["deterministic_shards"])
    frames = []
    shard_hashes = []
    for shard in range(shards):
        matches = list(args.input_root.rglob(f"rgfca_sharedness_v2_capacity_pilot_shard_{shard:02d}.csv"))
        manifests = list(args.input_root.rglob(f"rgfca_sharedness_v2_capacity_pilot_shard_{shard:02d}.json"))
        if len(matches) != 1 or len(manifests) != 1:
            raise RuntimeError(f"shard {shard} census is incomplete or duplicated")
        m = json.loads(manifests[0].read_text(encoding="utf-8"))
        if m["status"] != "complete_metadata_only_high_depth_capacity_pilot_shard" or int(m["shard_index"]) != shard:
            raise RuntimeError(f"bad shard manifest {shard}")
        if m["image_pixels_opened"] is not False or m["flower_colour_used"] is not False:
            raise RuntimeError("shard opened prohibited outcomes")
        if m["lineage"]["contract_sha256"] != sha256_file(CONTRACT):
            raise RuntimeError("shard contract lineage mismatch")
        if m["lineage"]["parent_audit_sha256"] != sha256_file(PARENT_AUDIT):
            raise RuntimeError("shard parent lineage mismatch")
        if m["lineage"]["csv_sha256"] != sha256_file(matches[0]):
            raise RuntimeError("shard CSV hash mismatch")
        frames.append(pd.read_csv(matches[0]))
        shard_hashes.append(sha256_file(matches[0]))
    frame = pd.concat(frames, ignore_index=True).sort_values("pilot_rank", kind="mergesort").reset_index(drop=True)
    n = int(c["pilot_species_selection"]["n"])
    if len(frame) != n or frame["pilot_rank"].tolist() != list(range(n)) or frame["inat_taxon_id"].nunique() != n:
        raise RuntimeError("fixed 500-species pilot census is incomplete")
    error_mask = frame["terminal_request_error"].fillna("").astype(str).str.len().gt(0)
    error_n = int(error_mask.sum())
    error_fraction = error_n / n
    targets = [int(x) for x in c["fixed_targets_after_observer_cap"]]
    counts = {str(t): int(frame[f"confirmed_at_least_{t}"].sum()) for t in targets}
    below = {str(t): int(frame[f"confirmed_below_{t}"].sum()) for t in targets}
    indeterminate = {str(t): int(frame[f"indeterminate_{t}"].sum()) for t in targets}
    primary_target = 300
    threshold = 150
    evaluable = error_fraction <= 0.05
    primary_pass = bool(evaluable and counts[str(primary_target)] >= threshold)
    result = {
        "protocol": c["protocol"],
        "status": "complete_metadata_only_high_depth_capacity_pilot" if evaluable else "not_evaluable_high_depth_capacity_pilot_due_request_failure",
        "pilot_species": n,
        "terminal_request_error_species": error_n,
        "terminal_request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": 0.05,
        "confirmed_at_least_after_observer_cap": counts,
        "confirmed_below_after_api_exhaustion": below,
        "indeterminate_due_error_or_page_cap": indeterminate,
        "primary_feasibility_target_after_observer_cap": primary_target,
        "primary_required_species": threshold,
        "primary_feasibility_pass": primary_pass,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "ecological_claim_changed": False,
        "reserve_flower_specific_gate_reclassified": False,
        "next_step_if_pass": "sharedness_v2_synthetic_power_design_with_no_observed_colour_opening",
        "next_step_if_fail_or_indeterminate": "do_not_relax_target_post_outcome; either prefreeze a larger metadata pilot or conclude current high-depth design is not established",
        "claim_ceiling": c["claim_ceiling"],
        "lineage": {"contract_sha256": sha256_file(CONTRACT), "parent_audit_sha256": sha256_file(PARENT_AUDIT), "shard_csv_sha256": shard_hashes},
    }
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_csv, index=False, lineterminator="\n")
    result["combined_csv_sha256"] = sha256_file(args.output_csv)
    args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
