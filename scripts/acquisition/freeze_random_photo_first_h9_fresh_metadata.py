#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h9_individual_distance_contract_v1.json"
FEASIBILITY = ROOT / "docs/supporting/random_photo_first_h9_metadata_feasibility_result_v1.json"
SPECIES = ROOT / "data/frozen/random_photo_first_h9_metadata_species_v1.csv"
EXCLUSION = ROOT / "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv"
EXCLUSION_MANIFEST = ROOT / "docs/supporting/random_photo_first_h9_exclusion_manifest_v1.json"
OUT = ROOT / "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv"
AUDIT = ROOT / "data/frozen/random_photo_first_h9_fresh_metadata_species_audit_v1.csv"
MANIFEST = ROOT / "docs/supporting/random_photo_first_h9_fresh_metadata_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if OUT.exists() or AUDIT.exists() or MANIFEST.exists():
        raise RuntimeError("H9 durable fresh metadata output already exists; one-shot query rerun is forbidden")
    contract = json.loads(CONTRACT.read_text())
    feasibility = json.loads(FEASIBILITY.read_text())
    exclusion_manifest = json.loads(EXCLUSION_MANIFEST.read_text())
    if contract["status"] != "prospective_h9_frozen_after_metadata_feasibility_before_any_h9_query_or_pixel":
        raise RuntimeError("H9 contract is not in prequery state")
    if feasibility["automatic_selection"]["h9_metadata_feasible"] is not True:
        raise RuntimeError("H9 metadata feasibility did not pass")
    species = pd.read_csv(SPECIES)
    if len(species) != int(contract["design_lineage"]["expected_species"]):
        raise RuntimeError("H9 species-frame size mismatch")
    exclusion = pd.read_csv(EXCLUSION)
    if sha256_file(EXCLUSION) != exclusion_manifest["h9_exclusion_sha256"]:
        raise RuntimeError("H9 exclusion ledger lineage mismatch")
    excluded_obs = set(exclusion["observation_id"].astype(int).tolist())
    excluded_photo = set(exclusion["photo_id"].astype(int).tolist())

    acq = contract["fresh_metadata_acquisition"]
    client = InaturalistObservationClient(
        request_interval_seconds=1.05,
        timeout_seconds=45.0,
        max_retries=int(acq["request_retries"]),
        user_agent="fcp-random-photo-first-h9/1.0 (github.com/zuizui0223/fcp)",
    )
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=species[["species", "inat_taxon_id"]],
        exclusion_observation_ids=excluded_obs,
        exclusion_photo_ids=excluded_photo,
        per_page=int(acq["per_page"]),
        observer_cap_n=int(acq["observer_cap"]),
        fixed_raw_photos=int(acq["fixed_raw_photos_per_species"]),
        maximum_positional_accuracy_m=int(acq["maximum_positional_accuracy_m"]),
        allowed_photo_licenses=tuple(acq["allowed_photo_licenses"]),
    )
    if int(frozen.manifest["query_attempts"]) != int(acq["queries"]):
        raise RuntimeError("H9 did not execute exactly the frozen query count")
    if len(frozen.observations):
        if set(frozen.observations["observation_id"].astype(int)) & excluded_obs:
            raise RuntimeError("prior observation entered H9")
        if set(frozen.observations["photo_id"].astype(int)) & excluded_photo:
            raise RuntimeError("prior photo entered H9")
        if frozen.observations["observation_id"].duplicated().any() or frozen.observations["photo_id"].duplicated().any():
            raise RuntimeError("H9 retained IDs are not unique")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    frozen.observations.to_csv(OUT, index=False, lineterminator="\n")
    frozen.species_audit.to_csv(AUDIT, index=False, lineterminator="\n")
    gate = contract["premeasurement_gate"]
    full_species = int(frozen.manifest["full_fixed_n_species"])
    required = int(gate["minimum_full_species"])
    passed = full_species >= required
    manifest = {
        "protocol": contract["protocol"],
        "status": "h9_fresh_metadata_frozen_before_pixels",
        "query_attempts": int(frozen.manifest["query_attempts"]),
        "expected_query_attempts": int(acq["queries"]),
        "query_retries": int(acq["request_retries"]),
        "request_error_species": int(frozen.manifest["request_errors"]),
        "selected_species": int(frozen.manifest["selected_species"]),
        "retained_fresh_photos": int(frozen.manifest["retained_fresh_photos"]),
        "retained_unique_observation_ids": int(frozen.observations["observation_id"].nunique()) if len(frozen.observations) else 0,
        "retained_unique_photo_ids": int(frozen.observations["photo_id"].nunique()) if len(frozen.observations) else 0,
        "premeasurement_gate": {
            "fixed_raw_photos_per_species": int(acq["fixed_raw_photos_per_species"]),
            "full_species": full_species,
            "required_species": required,
            "pass": bool(passed),
            "decision": "authorize_h9_location_blind_measurement" if passed else str(gate["if_fail"]),
        },
        "freshness_firewall": {
            "prior_exclusion_observation_ids": len(excluded_obs),
            "prior_exclusion_photo_ids": len(excluded_photo),
            "prior_observations_retained": 0,
            "prior_photos_retained": 0,
        },
        "candidate_image_pixels_opened": False,
        "colour_used_for_selection": False,
        "query_opened_once": True,
        "query_completed_without_replacement": True,
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "feasibility_result_sha256": sha256_file(FEASIBILITY),
            "selected_species_sha256": sha256_file(SPECIES),
            "exclusion_table_sha256": sha256_file(EXCLUSION),
        },
        "files": {
            "observations": {"path": str(OUT.relative_to(ROOT)), "sha256": sha256_file(OUT)},
            "species_audit": {"path": str(AUDIT.relative_to(ROOT)), "sha256": sha256_file(AUDIT)},
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
