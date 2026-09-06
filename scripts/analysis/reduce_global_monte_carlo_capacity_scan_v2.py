#!/usr/bin/env python3
"""Strictly reassemble the eight frozen metadata-only capacity-scan shards."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_v2_execution_contract.json"
V2_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_v2_manifest_v1.json"
SPECIES_FRAME = ROOT / "data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv"
OUT_AUDIT = ROOT / "data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv"
OUT_SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv"
OUT_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    a = args()
    for out in (OUT_AUDIT, OUT_SELECTED, OUT_MANIFEST):
        if out.exists():
            raise RuntimeError(f"refusing to overwrite frozen output: {out}")

    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    v2 = json.loads(V2_MANIFEST.read_text(encoding="utf-8"))
    if execution.get("biological_rules_changed") is not False:
        raise RuntimeError("execution contract changed biological rules")
    if v2.get("status") != "complete_cache_resistant_metadata_only_global_species_discovery_v2":
        raise RuntimeError("V2 discovery is not complete")
    if v2.get("candidate_image_pixels_opened") is not False or v2.get("flower_colour_used") is not False:
        raise RuntimeError("V2 discovery opened forbidden outcomes")

    shard_count = int(execution["execution"]["deterministic_shards"])
    parts: list[pd.DataFrame] = []
    shard_manifests: list[dict[str, object]] = []
    for shard_index in range(shard_count):
        csv_matches = list(a.input_root.rglob(f"capacity_scan_shard_{shard_index:02d}.csv"))
        json_matches = list(a.input_root.rglob(f"capacity_scan_shard_{shard_index:02d}.json"))
        if len(csv_matches) != 1 or len(json_matches) != 1:
            raise RuntimeError(f"shard {shard_index} requires exactly one CSV and one manifest")
        m = json.loads(json_matches[0].read_text(encoding="utf-8"))
        if m.get("status") != "complete_metadata_only_capacity_scan_shard":
            raise RuntimeError(f"shard {shard_index} status invalid")
        if int(m.get("shard_index", -1)) != shard_index or int(m.get("shard_count", -1)) != shard_count:
            raise RuntimeError(f"shard {shard_index} identity drift")
        if m.get("candidate_image_pixels_opened") is not False or m.get("flower_colour_used") is not False:
            raise RuntimeError(f"shard {shard_index} opened forbidden outcomes")
        if sha256_file(csv_matches[0]) != m["lineage"]["audit_sha256"]:
            raise RuntimeError(f"shard {shard_index} audit SHA mismatch")
        frame = pd.read_csv(csv_matches[0])
        if len(frame) != int(m["shard_species_count"]):
            raise RuntimeError(f"shard {shard_index} row-count mismatch")
        if frame["shard_index"].nunique() != 1 or int(frame["shard_index"].iloc[0]) != shard_index:
            raise RuntimeError(f"shard {shard_index} CSV identity mismatch")
        parts.append(frame)
        shard_manifests.append(m)

    audit = pd.concat(parts, ignore_index=True)
    if audit["inat_taxon_id"].duplicated().any():
        raise RuntimeError("duplicate taxon across capacity shards")
    if audit["global_row_index"].duplicated().any():
        raise RuntimeError("duplicate global row index across capacity shards")

    expected = pd.read_csv(SPECIES_FRAME)[["species", "inat_taxon_id"]].drop_duplicates(["inat_taxon_id"], keep="first")
    expected["inat_taxon_id"] = expected["inat_taxon_id"].astype(int)
    expected["species"] = expected["species"].astype(str)
    expected = expected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    expected["global_row_index"] = range(len(expected))
    got_ids = set(audit["inat_taxon_id"].astype(int))
    expected_ids = set(expected["inat_taxon_id"].astype(int))
    if got_ids != expected_ids:
        raise RuntimeError(f"capacity taxon set mismatch: missing={len(expected_ids-got_ids)} extra={len(got_ids-expected_ids)}")
    if len(audit) != len(expected):
        raise RuntimeError("capacity audit row count does not equal expected species count")

    audit = audit.sort_values("global_row_index", kind="mergesort").reset_index(drop=True)
    if audit["global_row_index"].astype(int).tolist() != expected["global_row_index"].astype(int).tolist():
        raise RuntimeError("capacity global row-index coverage drift")
    if audit["inat_taxon_id"].astype(int).tolist() != expected["inat_taxon_id"].astype(int).tolist():
        raise RuntimeError("capacity deterministic ordering drift")

    targets = [int(x) for x in parent["target_rule"]["candidate_raw_photos_per_species"]]
    target_counts = {str(t): int(audit[f"eligible_raw_{t}"].astype(bool).sum()) for t in targets}
    minimum_species = int(parent["target_rule"]["minimum_metadata_eligible_species"])
    errors = int(audit["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    error_fraction = float(errors / len(audit)) if len(audit) else 1.0
    ceiling = float(parent["request_failure_policy"]["maximum_error_fraction_for_capacity_decision"])
    request_coverage_ok = error_fraction <= ceiling
    selected_target = next((t for t in targets if target_counts[str(t)] >= minimum_species), None) if request_coverage_ok else None
    if not request_coverage_ok:
        status = str(parent["request_failure_policy"]["if_exceeded"])
    elif selected_target is None:
        status = str(parent["target_rule"]["if_none_pass"])
    else:
        status = "complete_metadata_only_capacity_scan_target_selected_v2"

    if selected_target is None:
        selected = audit.iloc[0:0][["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"]].copy()
    else:
        selected = audit.loc[
            audit[f"eligible_raw_{selected_target}"].astype(bool),
            ["species", "inat_taxon_id", "after_observer_cap", "maximum_span_km"],
        ].copy()
        selected = selected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
        selected["selected_raw_photo_target"] = int(selected_target)

    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUT_AUDIT, index=False, lineterminator="\n")
    selected.to_csv(OUT_SELECTED, index=False, lineterminator="\n")
    manifest = {
        "protocol": execution["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "shards_required": shard_count,
        "shards_reassembled": len(parts),
        "discovered_species_scanned": int(len(audit)),
        "request_attempts": int(len(audit)),
        "request_errors": errors,
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": ceiling,
        "request_coverage_ok": bool(request_coverage_ok),
        "target_counts": target_counts,
        "minimum_metadata_eligible_species": minimum_species,
        "selected_raw_photo_target": (None if selected_target is None else int(selected_target)),
        "selected_species": int(len(selected)),
        "actual_image_acquisition_authorized": False,
        "biological_rules_changed_from_parent": False,
        "lineage": {
            "parent_contract_sha256": sha256_file(PARENT),
            "execution_contract_sha256": sha256_file(EXECUTION),
            "v2_manifest_sha256": sha256_file(V2_MANIFEST),
            "combined_species_sha256": sha256_file(SPECIES_FRAME),
            "shard_manifest_audit_sha256": [m["lineage"]["audit_sha256"] for m in shard_manifests],
            "audit_sha256": sha256_file(OUT_AUDIT),
            "selected_sha256": sha256_file(OUT_SELECTED),
        },
        "files": {
            "species_audit": str(OUT_AUDIT.relative_to(ROOT)),
            "selected_species_frame": str(OUT_SELECTED.relative_to(ROOT)),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
