#!/usr/bin/env python3
"""Generate the exact metadata-only candidate acquisition authorization from one frozen capacity result.

The handoff is fail-closed. A normal v2 capacity success is used directly. If a v3
transport-recovery manifest exists, it must itself be the successful frozen recovery
result and becomes the only admissible source. The helper never chooses a target from
biological outcomes. It applies only the prospectively fixed, outcome-blind candidate
species budget before fresh candidate acquisition.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

from fcp_pipeline.global_measurement_budget import select_hashed_taxa

ROOT = Path(__file__).resolve().parents[2]
BRANCH = "analysis/global-monte-carlo-barrier-atlas"
V2_CAPACITY = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json"
V2_SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv"
V3_CAPACITY = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json"
V3_SELECTED = ROOT / "data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv"
BUDGET = ROOT / "docs/supporting/global_monte_carlo_candidate_species_budget_amendment_v1.json"
OUT = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_authorization_v1.json"


def git_blob(path: Path) -> str:
    rel = path.relative_to(ROOT).as_posix()
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{rel}"], text=True).strip()


def taxon_digest(values: pd.Series) -> str:
    text = "\n".join(str(int(x)) for x in sorted(set(values.astype(int)))) + "\n"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def choose_capacity_source() -> tuple[Path, Path, str]:
    if V3_CAPACITY.exists():
        if not V3_SELECTED.exists():
            raise RuntimeError("v3 capacity manifest exists without its selected-species frame")
        manifest = json.loads(V3_CAPACITY.read_text(encoding="utf-8"))
        if manifest.get("status") != "complete_metadata_only_capacity_scan_target_selected_after_transport_recovery_v3":
            raise RuntimeError("v3 recovery exists but did not yield an admissible capacity result")
        if manifest.get("original_v2_status") != "not_evaluable_capacity_scan_due_request_failure":
            raise RuntimeError("v3 recovery lineage does not originate from the frozen v2 request-failure state")
        if manifest.get("second_recovery_permitted") is not False or manifest.get("biological_rules_changed") is not False:
            raise RuntimeError("v3 recovery contract drifted")
        return V3_CAPACITY, V3_SELECTED, "v3_transport_recovery"

    if not V2_CAPACITY.exists() or not V2_SELECTED.exists():
        raise RuntimeError("no complete frozen capacity source is available")
    manifest = json.loads(V2_CAPACITY.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete_metadata_only_capacity_scan_target_selected_v2":
        raise RuntimeError("v2 capacity result did not pass the frozen request-error and target-selection gates")
    if manifest.get("request_coverage_ok") is not True:
        raise RuntimeError("v2 capacity request coverage did not pass")
    return V2_CAPACITY, V2_SELECTED, "v2_primary"


def main() -> int:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing authorization: {OUT}")

    capacity_path, selected_path, source = choose_capacity_source()
    capacity = json.loads(capacity_path.read_text(encoding="utf-8"))
    budget = json.loads(BUDGET.read_text(encoding="utf-8"))
    if capacity.get("candidate_image_pixels_opened") is not False or capacity.get("flower_colour_used") is not False:
        raise RuntimeError("capacity stage opened forbidden outcomes")
    if capacity.get("actual_image_acquisition_authorized") is not False:
        raise RuntimeError("capacity manifest unexpectedly self-authorized acquisition")
    if budget.get("status") != "frozen_after_metadata_capacity_scale_observed_before_candidate_acquisition_outcome_and_before_any_global_candidate_pixels":
        raise RuntimeError("candidate species budget amendment is not frozen pre-outcome")

    target = int(capacity.get("selected_raw_photo_target"))
    if target not in (100, 80, 60):
        raise RuntimeError("capacity target is outside the prospectively frozen set")
    capacity_selected_n = int(capacity.get("selected_species") or 0)
    if capacity_selected_n < int(capacity.get("minimum_metadata_eligible_species") or 300):
        raise RuntimeError("capacity selected fewer than the frozen minimum species")

    selected = pd.read_csv(selected_path)
    if len(selected) != capacity_selected_n or selected["inat_taxon_id"].nunique() != capacity_selected_n:
        raise RuntimeError("selected-species frame does not match capacity manifest")
    if "selected_raw_photo_target" not in selected.columns:
        raise RuntimeError("selected-species frame lacks selected_raw_photo_target")
    if not (selected["selected_raw_photo_target"].astype(int) == target).all():
        raise RuntimeError("selected-species frame target differs from capacity manifest")

    b = budget["candidate_species_budget"]
    maximum_species = int(b["maximum_species"])
    seed = int(b["selection_seed"])
    taxa = set(select_hashed_taxa(selected["inat_taxon_id"], maximum_species=maximum_species, seed=seed))
    bounded = selected.loc[selected["inat_taxon_id"].astype(int).isin(taxa)].copy()
    bounded_n = int(len(bounded))
    if bounded_n != min(capacity_selected_n, maximum_species):
        raise RuntimeError("bounded candidate species selection produced the wrong denominator")

    authorization = {
        "status": "authorize_exactly_one_metadata_only_global_candidate_acquisition",
        "branch": BRANCH,
        "capacity_source": source,
        "capacity_manifest_path": capacity_path.relative_to(ROOT).as_posix(),
        "selected_species_path": selected_path.relative_to(ROOT).as_posix(),
        "capacity_manifest_blob_sha": git_blob(capacity_path),
        "selected_species_blob_sha": git_blob(selected_path),
        "candidate_species_budget_amendment_path": BUDGET.relative_to(ROOT).as_posix(),
        "candidate_species_budget_amendment_blob_sha": git_blob(BUDGET),
        "selected_raw_photo_target": target,
        "capacity_selected_species": capacity_selected_n,
        "candidate_species_budget": maximum_species,
        "candidate_species_selection_seed": seed,
        "selected_species": bounded_n,
        "candidate_taxon_id_sha256": taxon_digest(bounded["inat_taxon_id"]),
        "candidate_image_pixels_may_open": False,
        "flower_colour_may_open": False,
        "target_relaxation_allowed": False,
        "additional_pages_after_result_allowed": False,
        "rerun_for_favourable_species_set_allowed": False,
    }
    OUT.write_text(json.dumps(authorization, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(authorization, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
