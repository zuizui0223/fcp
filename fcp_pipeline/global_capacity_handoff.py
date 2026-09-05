"""Resolve the one authorized capacity source for bounded global candidate acquisition.

A normal v2 capacity result and a successful v3 transport-recovery result are both
admissible, but only when the authorization names the exact source paths and Git
blob SHAs. The full capacity-selected set remains immutable provenance; only a
prospectively fixed hash-ranked subset is handed to the fresh candidate stage.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess

import pandas as pd

from .global_measurement_budget import select_hashed_taxa

BRANCH = "analysis/global-monte-carlo-barrier-atlas"
BUDGET_PATH = "docs/supporting/global_monte_carlo_candidate_species_budget_amendment_v1.json"
ALLOWED_SOURCES = {
    "v2_primary": {
        "manifest": "docs/supporting/global_monte_carlo_capacity_scan_manifest_v2.json",
        "selected": "data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv",
        "status": "complete_metadata_only_capacity_scan_target_selected_v2",
    },
    "v3_transport_recovery": {
        "manifest": "docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json",
        "selected": "data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv",
        "status": "complete_metadata_only_capacity_scan_target_selected_after_transport_recovery_v3",
    },
}


@dataclass(frozen=True)
class CapacityHandoff:
    source: str
    manifest_path: Path
    selected_path: Path
    manifest: dict[str, object]
    selected: pd.DataFrame
    target: int
    capacity_selected_species: int
    selected_species: int


def git_blob_sha(root: Path, path: Path) -> str:
    rel = path.relative_to(root).as_posix()
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{rel}"], text=True).strip()


def taxon_digest(values: pd.Series) -> str:
    text = "\n".join(str(int(x)) for x in sorted(set(values.astype(int)))) + "\n"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def resolve_capacity_handoff(
    root: Path,
    authorization_path: Path,
    *,
    minimum_species: int = 300,
) -> CapacityHandoff:
    root = root.resolve()
    authorization_path = authorization_path.resolve()
    auth = json.loads(authorization_path.read_text(encoding="utf-8"))
    if auth.get("status") != "authorize_exactly_one_metadata_only_global_candidate_acquisition":
        raise RuntimeError("candidate acquisition is not explicitly authorized")
    if auth.get("branch") != BRANCH:
        raise RuntimeError("candidate acquisition authorization branch mismatch")
    if auth.get("candidate_image_pixels_may_open") is not False or auth.get("flower_colour_may_open") is not False:
        raise RuntimeError("candidate acquisition authorization unexpectedly opens pixels/colour")
    if auth.get("target_relaxation_allowed") is not False or auth.get("additional_pages_after_result_allowed") is not False:
        raise RuntimeError("candidate authorization permits result-dependent sampling changes")
    if auth.get("rerun_for_favourable_species_set_allowed") is not False:
        raise RuntimeError("candidate authorization permits favourable reruns")

    source = str(auth.get("capacity_source") or "")
    spec = ALLOWED_SOURCES.get(source)
    if spec is None:
        raise RuntimeError("candidate authorization names an unsupported capacity source")
    if auth.get("capacity_manifest_path") != spec["manifest"] or auth.get("selected_species_path") != spec["selected"]:
        raise RuntimeError("candidate authorization capacity paths do not match the named source")

    manifest_path = root / str(spec["manifest"])
    selected_path = root / str(spec["selected"])
    budget_path = root / BUDGET_PATH
    if not manifest_path.exists() or not selected_path.exists() or not budget_path.exists():
        raise RuntimeError("authorized capacity/budget source files are missing")
    if auth.get("capacity_manifest_blob_sha") != git_blob_sha(root, manifest_path):
        raise RuntimeError("authorized capacity manifest blob SHA mismatch")
    if auth.get("selected_species_blob_sha") != git_blob_sha(root, selected_path):
        raise RuntimeError("authorized selected-species blob SHA mismatch")
    if auth.get("candidate_species_budget_amendment_path") != BUDGET_PATH:
        raise RuntimeError("candidate authorization budget path drifted")
    if auth.get("candidate_species_budget_amendment_blob_sha") != git_blob_sha(root, budget_path):
        raise RuntimeError("candidate species budget blob SHA mismatch")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    budget = json.loads(budget_path.read_text(encoding="utf-8"))
    if budget.get("status") != "frozen_after_metadata_capacity_scale_observed_before_candidate_acquisition_outcome_and_before_any_global_candidate_pixels":
        raise RuntimeError("candidate species budget is not frozen pre-outcome")
    if manifest.get("status") != spec["status"]:
        raise RuntimeError("authorized capacity source is not a successful frozen capacity result")
    if manifest.get("candidate_image_pixels_opened") is not False or manifest.get("flower_colour_used") is not False:
        raise RuntimeError("authorized capacity result opened forbidden outcomes")
    if manifest.get("actual_image_acquisition_authorized") is not False:
        raise RuntimeError("capacity result unexpectedly self-authorized candidate acquisition")
    if source == "v3_transport_recovery":
        if manifest.get("original_v2_status") != "not_evaluable_capacity_scan_due_request_failure":
            raise RuntimeError("v3 capacity source lacks the required failed-v2 lineage")
        if manifest.get("second_recovery_permitted") is not False or manifest.get("biological_rules_changed") is not False:
            raise RuntimeError("v3 capacity recovery rules drifted")

    target = int(manifest.get("selected_raw_photo_target"))
    if target not in (100, 80, 60):
        raise RuntimeError("authorized capacity target is outside the frozen set")
    capacity_selected_n = int(manifest.get("selected_species") or 0)
    if capacity_selected_n < int(minimum_species):
        raise RuntimeError("authorized capacity source has fewer than the frozen minimum species")
    if int(auth.get("selected_raw_photo_target")) != target or int(auth.get("capacity_selected_species")) != capacity_selected_n:
        raise RuntimeError("candidate authorization capacity denominator differs from capacity result")

    full_selected = pd.read_csv(selected_path)
    if not {"species", "inat_taxon_id", "selected_raw_photo_target"}.issubset(full_selected.columns):
        raise RuntimeError("selected-species frame lacks required handoff columns")
    if len(full_selected) != capacity_selected_n or full_selected["inat_taxon_id"].nunique() != capacity_selected_n:
        raise RuntimeError("selected-species frame denominator differs from capacity result")
    if not (full_selected["selected_raw_photo_target"].astype(int) == target).all():
        raise RuntimeError("selected-species frame target differs from capacity result")

    b = budget["candidate_species_budget"]
    maximum_species = int(b["maximum_species"])
    seed = int(b["selection_seed"])
    taxa = select_hashed_taxa(
        full_selected["inat_taxon_id"].astype(int),
        maximum_species=maximum_species,
        seed=seed,
    )
    selected_set = set(taxa)
    selected = full_selected.loc[full_selected["inat_taxon_id"].astype(int).isin(selected_set)].copy()
    selected = selected.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    selected_n = int(len(selected))
    if selected_n != min(capacity_selected_n, maximum_species):
        raise RuntimeError("bounded candidate species selection produced the wrong denominator")
    if int(auth.get("candidate_species_budget")) != maximum_species or int(auth.get("candidate_species_selection_seed")) != seed:
        raise RuntimeError("candidate authorization species-budget rule drifted")
    if int(auth.get("selected_species")) != selected_n:
        raise RuntimeError("candidate authorization bounded species count drifted")
    if auth.get("candidate_taxon_id_sha256") != taxon_digest(selected["inat_taxon_id"]):
        raise RuntimeError("candidate authorization taxon subset differs from frozen hash selection")

    return CapacityHandoff(
        source=source,
        manifest_path=manifest_path,
        selected_path=selected_path,
        manifest=manifest,
        selected=selected,
        target=target,
        capacity_selected_species=capacity_selected_n,
        selected_species=selected_n,
    )


__all__ = ["ALLOWED_SOURCES", "CapacityHandoff", "resolve_capacity_handoff", "taxon_digest"]
