#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd

from fcp_pipeline.h7_acquisition import freeze_h7_fresh_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient
from fcp_pipeline.shared_transition_surface import EqualAreaGrid

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h7_balanced_itv_contract_v1.json"
AMENDMENT = ROOT / "docs/supporting/random_photo_first_h7_freshness_source_amendment_v1.json"
TARGETS = ROOT / "data/frozen/random_photo_first_h7_species_cells_v1.csv"
TARGET_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_target_manifest_v1.json"
EXCLUSIONS = ROOT / "data/frozen/random_photo_first_h7_exclusion_ledger_v1.csv"
EXCLUSION_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_exclusion_manifest_v1.json"

OUT_OBS = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_v1.csv"
OUT_AUDIT = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_target_audit_v1.csv"
OUT_SUPPORT = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_species_support_v1.csv"
OUT_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_fresh_metadata_manifest_v1.json"
OUTPUTS = (OUT_OBS, OUT_AUDIT, OUT_SUPPORT, OUT_MANIFEST)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv_atomic(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def write_json_atomic(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(tmp, path)


def main() -> None:
    # One-shot guard must execute before the network client is constructed.
    existing = [p.as_posix() for p in OUTPUTS if p.exists()]
    if existing:
        raise RuntimeError(f"H7 fresh metadata freeze is one-shot; durable outputs already exist: {existing}")

    contract = json.loads(CONTRACT.read_text())
    amendment = json.loads(AMENDMENT.read_text())
    target_manifest = json.loads(TARGET_MANIFEST.read_text())
    exclusion_manifest = json.loads(EXCLUSION_MANIFEST.read_text())
    if contract["status"] != "prospective_h7_frozen_after_metadata_feasibility_before_any_fresh_h7_query_or_pixel":
        raise RuntimeError("H7 prospective contract is not frozen")
    if amendment["status"] != "technical_lineage_correction_frozen_before_any_fresh_h7_query_or_pixel":
        raise RuntimeError("H7 provenance amendment is not frozen")
    if target_manifest["status"] != "h7_species_and_target_cells_frozen_before_fresh_query":
        raise RuntimeError("H7 target table is not frozen")
    if exclusion_manifest["status"] != "h7_prior_photo_exclusion_ledger_frozen_before_fresh_query":
        raise RuntimeError("H7 exclusion ledger is not frozen")
    if sha256_file(TARGETS) != target_manifest["target_table_sha256"]:
        raise RuntimeError("H7 target table hash mismatch")
    if sha256_file(EXCLUSIONS) != exclusion_manifest["exclusion_table_sha256"]:
        raise RuntimeError("H7 exclusion table hash mismatch")
    if target_manifest["outcome_firewall"]["fresh_h7_api_queries_opened"] is not False:
        raise RuntimeError("target manifest says H7 query already opened")
    if exclusion_manifest["fresh_h7_api_queries_opened"] is not False:
        raise RuntimeError("exclusion manifest says H7 query already opened")

    q = contract["fresh_metadata_acquisition"]
    gate = contract["premeasurement_gate"]
    grid_spec = contract["target_cells"]["grid"]
    targets = pd.read_csv(TARGETS)
    exclusions = pd.read_csv(EXCLUSIONS)
    grid = EqualAreaGrid(n_lon=int(grid_spec["n_lon"]), n_sinlat=int(grid_spec["n_sinlat"]))

    # max_retries=0 is binding: each frozen target gets one random page attempt only.
    client = InaturalistObservationClient(
        user_agent="fcp-h7-balanced-itv/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=1.05,
        timeout_seconds=45.0,
        max_retries=0,
    )
    result = freeze_h7_fresh_metadata(
        client=client,
        targets=targets,
        exclusions=exclusions,
        grid=grid,
        per_page=int(q["per_page"]),
        observer_cap=int(q["observer_cap_per_species_cell"]),
        retained_cap=int(q["retained_photo_cap_per_species_cell"]),
        selection_seed=int(q["retention_seed"]),
        expected_species=int(contract["species_discovery"]["expected_species"]),
        expected_targets=int(q["expected_queries"]),
        required_species_for_gate=int(gate["required_species"]),
        required_full_cells_per_species=5,
        maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
        flowering_term_id=int(q["flowering_annotation"]["term_id"]),
        flowering_term_value_id=int(q["flowering_annotation"]["term_value_id"]),
        allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
    )

    write_csv_atomic(result.observations, OUT_OBS)
    write_csv_atomic(result.target_audit, OUT_AUDIT)
    write_csv_atomic(result.species_support, OUT_SUPPORT)
    manifest = dict(result.manifest)
    manifest.update(
        {
            "status": "h7_fresh_metadata_frozen_before_pixels",
            "query_opened_once": True,
            "query_completed_without_replacement": True,
            "lineage": {
                "contract_path": CONTRACT.relative_to(ROOT).as_posix(),
                "freshness_amendment_path": AMENDMENT.relative_to(ROOT).as_posix(),
                "target_table_sha256": sha256_file(TARGETS),
                "exclusion_table_sha256": sha256_file(EXCLUSIONS),
                "pr21_terminal_manifest_sha256": exclusion_manifest["pr21_inner_manifest_sha256"],
            },
        }
    )
    # Hash CSV outputs before writing the final manifest.
    manifest["files"] = {
        "observations": {"path": OUT_OBS.relative_to(ROOT).as_posix(), "sha256": sha256_file(OUT_OBS)},
        "target_audit": {"path": OUT_AUDIT.relative_to(ROOT).as_posix(), "sha256": sha256_file(OUT_AUDIT)},
        "species_support": {"path": OUT_SUPPORT.relative_to(ROOT).as_posix(), "sha256": sha256_file(OUT_SUPPORT)},
    }
    write_json_atomic(manifest, OUT_MANIFEST)
    print(json.dumps(manifest, indent=2))
    if not manifest["premeasurement_gate"]["pass"]:
        print("H7 premeasurement gate did not pass; candidate image pixels MUST remain unopened.")


if __name__ == "__main__":
    main()
