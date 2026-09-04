#!/usr/bin/env python3
"""Strictly reduce all frozen global candidate-acquisition shards before pixels open."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_capacity_handoff import resolve_capacity_handoff

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_contract_v1.json"
AUTHORIZATION = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_authorization_v1.json"
OUT_AUDIT = ROOT / "data/frozen/global_monte_carlo_candidate_species_audit_v1.csv"
OUT_CANDIDATES = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
OUT_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (OUT_AUDIT, OUT_CANDIDATES, OUT_MANIFEST):
        if path.exists():
            raise RuntimeError(f"refusing to overwrite frozen output: {path}")

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    handoff = resolve_capacity_handoff(
        ROOT,
        AUTHORIZATION,
        minimum_species=int(contract["premeasurement_gate"]["minimum_full_target_species"]),
    )
    capacity = handoff.manifest
    target = int(handoff.target)
    expected = handoff.selected[["species", "inat_taxon_id"]].drop_duplicates("inat_taxon_id", keep="first").copy()
    expected["inat_taxon_id"] = expected["inat_taxon_id"].astype(int)
    expected["species"] = expected["species"].astype(str)
    expected = expected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    expected["global_row_index"] = range(len(expected))
    if len(expected) != int(capacity.get("selected_species") or 0):
        raise RuntimeError("capacity selected-species lineage mismatch")

    shard_count = int(contract["execution"]["deterministic_shards"])
    audit_parts: list[pd.DataFrame] = []
    photo_parts: list[pd.DataFrame] = []
    shard_manifests: list[dict[str, object]] = []
    for shard_index in range(shard_count):
        audit_matches = list(args.input_root.rglob(f"candidate_acquisition_shard_{shard_index:02d}_species.csv"))
        photo_matches = list(args.input_root.rglob(f"candidate_acquisition_shard_{shard_index:02d}_photos.csv"))
        manifest_matches = list(args.input_root.rglob(f"candidate_acquisition_shard_{shard_index:02d}.json"))
        if len(audit_matches) != 1 or len(photo_matches) != 1 or len(manifest_matches) != 1:
            raise RuntimeError(f"candidate shard {shard_index} requires exactly one audit, photo file and manifest")
        manifest = json.loads(manifest_matches[0].read_text(encoding="utf-8"))
        if manifest.get("status") != "complete_metadata_only_candidate_acquisition_shard":
            raise RuntimeError(f"candidate shard {shard_index} status invalid")
        if int(manifest.get("shard_index", -1)) != shard_index or int(manifest.get("shard_count", -1)) != shard_count:
            raise RuntimeError(f"candidate shard {shard_index} identity drift")
        if manifest.get("capacity_source") != handoff.source:
            raise RuntimeError(f"candidate shard {shard_index} capacity source drift")
        if int(manifest.get("selected_raw_photo_target", -1)) != target:
            raise RuntimeError(f"candidate shard {shard_index} target drift")
        if manifest.get("candidate_image_pixels_opened") is not False or manifest.get("flower_colour_used") is not False:
            raise RuntimeError(f"candidate shard {shard_index} opened forbidden outcomes")
        if manifest["lineage"].get("capacity_manifest_sha256") != sha256_file(handoff.manifest_path):
            raise RuntimeError(f"candidate shard {shard_index} capacity manifest lineage drift")
        if manifest["lineage"].get("selected_species_sha256") != sha256_file(handoff.selected_path):
            raise RuntimeError(f"candidate shard {shard_index} selected-species lineage drift")
        if sha256_file(audit_matches[0]) != manifest["lineage"]["species_audit_sha256"]:
            raise RuntimeError(f"candidate shard {shard_index} audit SHA mismatch")
        if sha256_file(photo_matches[0]) != manifest["lineage"]["candidate_photos_sha256"]:
            raise RuntimeError(f"candidate shard {shard_index} photo SHA mismatch")
        audit = pd.read_csv(audit_matches[0])
        photos = pd.read_csv(photo_matches[0])
        if len(audit) != int(manifest["shard_species"]):
            raise RuntimeError(f"candidate shard {shard_index} species count mismatch")
        if len(photos) != int(manifest["candidate_rows"]):
            raise RuntimeError(f"candidate shard {shard_index} candidate-row mismatch")
        audit_parts.append(audit)
        photo_parts.append(photos)
        shard_manifests.append(manifest)

    audit = pd.concat(audit_parts, ignore_index=True)
    photos = pd.concat(photo_parts, ignore_index=True)
    if audit["inat_taxon_id"].duplicated().any() or audit["global_row_index"].duplicated().any():
        raise RuntimeError("duplicate species/global index across candidate shards")
    got_ids = set(audit["inat_taxon_id"].astype(int))
    expected_ids = set(expected["inat_taxon_id"].astype(int))
    if got_ids != expected_ids or len(audit) != len(expected):
        raise RuntimeError("candidate acquisition did not cover the exact capacity-selected taxon set")
    audit = audit.sort_values("global_row_index", kind="mergesort").reset_index(drop=True)
    if audit["inat_taxon_id"].astype(int).tolist() != expected["inat_taxon_id"].astype(int).tolist():
        raise RuntimeError("candidate acquisition deterministic species ordering drifted")

    if len(photos):
        if photos["observation_id"].duplicated().any() or photos["photo_id"].duplicated().any():
            raise RuntimeError("duplicate observation/photo IDs across final candidate pool")
        photos = photos.sort_values(["global_row_index", "row_hash"], kind="mergesort").reset_index(drop=True)

    full_species = int(audit["full_target"].astype(bool).sum())
    expected_candidate_rows = full_species * target
    if len(photos) != expected_candidate_rows:
        raise RuntimeError(f"candidate row total {len(photos)} != full_species*target {expected_candidate_rows}")
    per_species = photos.groupby("inat_taxon_id", observed=True).size() if len(photos) else pd.Series(dtype=int)
    if len(per_species) and not (per_species.astype(int) == target).all():
        raise RuntimeError("not every retained candidate species contributes exactly the selected target")
    retained_ids = set(per_species.index.astype(int).tolist()) if len(per_species) else set()
    expected_retained = set(audit.loc[audit["full_target"].astype(bool), "inat_taxon_id"].astype(int).tolist())
    if retained_ids != expected_retained:
        raise RuntimeError("candidate photo species do not match full-target audit species")

    attempts = int(sum(int(m["request_attempts"]) for m in shard_manifests))
    errors = int(sum(int(m["request_errors"]) for m in shard_manifests))
    error_fraction = float(errors / attempts) if attempts else 1.0
    error_ceiling = float(contract["request_failure_policy"]["maximum_error_fraction_for_acquisition_gate"])
    minimum_species = int(contract["premeasurement_gate"]["minimum_full_target_species"])
    error_ok = error_fraction <= error_ceiling
    species_ok = full_species >= minimum_species
    gate_pass = bool(error_ok and species_ok)
    if not error_ok:
        status = str(contract["request_failure_policy"]["if_exceeded"])
        decision = status
    elif not species_ok:
        status = str(contract["premeasurement_gate"]["if_fewer_than_300"])
        decision = status
    else:
        status = "complete_global_candidate_acquisition_premeasurement_gate_passed"
        decision = "authorize_location_blind_candidate_pixel_measurement_only_under_separate_exact_authorization"

    OUT_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(OUT_AUDIT, index=False, lineterminator="\n")
    photos.to_csv(OUT_CANDIDATES, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "capacity_source": handoff.source,
        "capacity_manifest_path": str(handoff.manifest_path.relative_to(ROOT)),
        "capacity_selected_species_path": str(handoff.selected_path.relative_to(ROOT)),
        "shards_required": shard_count,
        "shards_reassembled": shard_count,
        "capacity_selected_raw_photo_target": target,
        "capacity_selected_species": int(len(expected)),
        "request_attempts": attempts,
        "request_errors": errors,
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": error_ceiling,
        "full_target_species": full_species,
        "minimum_full_target_species": minimum_species,
        "candidate_rows": int(len(photos)),
        "premeasurement_gate": {
            "request_coverage_pass": bool(error_ok),
            "species_count_pass": bool(species_ok),
            "pass": gate_pass,
            "decision": decision,
        },
        "measurement_authorized": False,
        "target_relaxed": False,
        "additional_pages_after_result": False,
        "lineage": {
            "candidate_contract_sha256": sha256_file(CONTRACT),
            "capacity_manifest_sha256": sha256_file(handoff.manifest_path),
            "capacity_selected_species_sha256": sha256_file(handoff.selected_path),
            "shard_species_audit_sha256": [m["lineage"]["species_audit_sha256"] for m in shard_manifests],
            "shard_candidate_photos_sha256": [m["lineage"]["candidate_photos_sha256"] for m in shard_manifests],
            "species_audit_sha256": sha256_file(OUT_AUDIT),
            "candidate_photos_sha256": sha256_file(OUT_CANDIDATES),
        },
        "files": {
            "species_audit": str(OUT_AUDIT.relative_to(ROOT)),
            "candidate_photo_pool": str(OUT_CANDIDATES.relative_to(ROOT)),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
