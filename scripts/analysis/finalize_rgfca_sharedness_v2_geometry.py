#!/usr/bin/env python3
"""Strictly combine and freeze the complete 150 x 300 RGFCA v2 metadata geometry."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_geometry_freeze_contract_v1.json"
CAPACITY_RESULT = ROOT / "docs/supporting/rgfca_sharedness_v2_high_depth_capacity_pilot_result_v1.json"
CAPACITY_ROWS = ROOT / "data/frozen/rgfca_sharedness_v2_high_depth_capacity_pilot_v1.csv"
SELECTION_RULE = ROOT / "docs/supporting/rgfca_sharedness_v2_post_capacity_selection_rule_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-photos", type=Path, required=True)
    p.add_argument("--output-audit", type=Path, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    args = p.parse_args()
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    shard_count = int(c["deterministic_shards"])
    photo_frames, audit_frames, shard_hashes = [], [], []
    for shard in range(shard_count):
        photos = list(args.input_root.rglob(f"rgfca_sharedness_v2_geometry_photos_shard_{shard:02d}.csv"))
        audits = list(args.input_root.rglob(f"rgfca_sharedness_v2_geometry_audit_shard_{shard:02d}.csv"))
        manifests = list(args.input_root.rglob(f"rgfca_sharedness_v2_geometry_shard_{shard:02d}.json"))
        if len(photos) != 1 or len(audits) != 1 or len(manifests) != 1:
            raise RuntimeError(f"shard {shard} census missing or duplicated")
        m = json.loads(manifests[0].read_text(encoding="utf-8"))
        if m["protocol"] != c["protocol"] or m["status"] != "complete_metadata_only_geometry_shard" or int(m["shard_index"]) != shard:
            raise RuntimeError(f"bad geometry shard manifest {shard}")
        if m["image_pixels_opened"] is not False or m["flower_colour_opened"] is not False:
            raise RuntimeError("geometry shard opened prohibited outcomes")
        expected = {
            "contract_sha256": sha256_file(CONTRACT),
            "capacity_result_sha256": sha256_file(CAPACITY_RESULT),
            "capacity_rows_sha256": sha256_file(CAPACITY_ROWS),
            "selection_rule_sha256": sha256_file(SELECTION_RULE),
            "photos_csv_sha256": sha256_file(photos[0]),
            "audit_csv_sha256": sha256_file(audits[0]),
        }
        for key, value in expected.items():
            if m["lineage"][key] != value:
                raise RuntimeError(f"shard {shard} lineage mismatch for {key}")
        photo_frames.append(pd.read_csv(photos[0]))
        audit_frames.append(pd.read_csv(audits[0]))
        shard_hashes.append({"shard": shard, "photos": expected["photos_csv_sha256"], "audit": expected["audit_csv_sha256"]})

    photos = pd.concat(photo_frames, ignore_index=True).sort_values(["species_order", "photo_order"], kind="mergesort").reset_index(drop=True)
    audit = pd.concat(audit_frames, ignore_index=True).sort_values("species_order", kind="mergesort").reset_index(drop=True)
    all_species = int(c["activation"]["selected_species"])
    target = int(c["fixed_photo_frame"]["photos_per_species"])
    errors = audit["terminal_request_error"].fillna("").astype(str).str.len().gt(0)
    full = bool(
        len(audit) == all_species
        and audit["species_order"].tolist() == list(range(all_species))
        and audit["inat_taxon_id"].nunique() == all_species
        and int(errors.sum()) == 0
        and audit["full_fixed_300"].astype(bool).all()
        and len(photos) == all_species * target
        and photos["photo_id"].nunique() == len(photos)
        and photos["observation_id"].nunique() == len(photos)
        and photos.groupby("inat_taxon_id").size().eq(target).all()
        and (photos["sample_role"] == "training").sum() == 75 * target
        and (photos["sample_role"] == "evaluation").sum() == 75 * target
    )
    status = "complete_frozen_metadata_geometry_150x300" if full else "geometry_freeze_failed_incomplete_fixed_frame"
    args.output_photos.parent.mkdir(parents=True, exist_ok=True)
    args.output_audit.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    photos.to_csv(args.output_photos, index=False, lineterminator="\n")
    audit.to_csv(args.output_audit, index=False, lineterminator="\n")
    result = {
        "protocol": c["protocol"],
        "status": status,
        "complete_fixed_frame": full,
        "species": int(len(audit)),
        "training_species": int(audit.loc[audit["sample_role"] == "training", "inat_taxon_id"].nunique()),
        "evaluation_species": int(audit.loc[audit["sample_role"] == "evaluation", "inat_taxon_id"].nunique()),
        "target_photos_per_species": target,
        "retained_photos": int(len(photos)),
        "terminal_request_error_species": int(errors.sum()),
        "species_below_target": int((~audit["full_fixed_300"].astype(bool)).sum()),
        "image_pixels_opened": False,
        "flower_colour_opened": False,
        "synthetic_sharedness_qualification_permitted": full,
        "image_acquisition_permitted": False,
        "ecological_claim_changed": False,
        "reserve_flower_specific_gate_reclassified": False,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "capacity_result_sha256": sha256_file(CAPACITY_RESULT),
            "selection_rule_sha256": sha256_file(SELECTION_RULE),
            "combined_photos_sha256": sha256_file(args.output_photos),
            "combined_audit_sha256": sha256_file(args.output_audit),
            "shards": shard_hashes
        },
        "claim_ceiling": c["claim_ceiling"]
    }
    args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
