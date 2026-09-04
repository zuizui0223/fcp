"""Resolve the one authorized capacity source for global candidate acquisition.

A normal v2 capacity result and a successful v3 transport-recovery result are both
admissible, but only when the authorization names the exact source paths and Git
blob SHAs. This module keeps that technical handoff identical across gate, shard,
and reducer code.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess

import pandas as pd

BRANCH = "analysis/global-monte-carlo-barrier-atlas"
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
    selected_species: int


def git_blob_sha(root: Path, path: Path) -> str:
    rel = path.relative_to(root).as_posix()
    return subprocess.check_output(["git", "rev-parse", f"HEAD:{rel}"], text=True).strip()


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
    if not manifest_path.exists() or not selected_path.exists():
        raise RuntimeError("authorized capacity source files are missing")
    if auth.get("capacity_manifest_blob_sha") != git_blob_sha(root, manifest_path):
        raise RuntimeError("authorized capacity manifest blob SHA mismatch")
    if auth.get("selected_species_blob_sha") != git_blob_sha(root, selected_path):
        raise RuntimeError("authorized selected-species blob SHA mismatch")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
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
    selected_n = int(manifest.get("selected_species") or 0)
    if selected_n < int(minimum_species):
        raise RuntimeError("authorized capacity source has fewer than the frozen minimum species")
    if int(auth.get("selected_raw_photo_target")) != target or int(auth.get("selected_species")) != selected_n:
        raise RuntimeError("candidate authorization denominator differs from capacity result")

    selected = pd.read_csv(selected_path)
    if not {"species", "inat_taxon_id", "selected_raw_photo_target"}.issubset(selected.columns):
        raise RuntimeError("selected-species frame lacks required handoff columns")
    if len(selected) != selected_n or selected["inat_taxon_id"].nunique() != selected_n:
        raise RuntimeError("selected-species frame denominator differs from capacity result")
    if not (selected["selected_raw_photo_target"].astype(int) == target).all():
        raise RuntimeError("selected-species frame target differs from capacity result")

    return CapacityHandoff(
        source=source,
        manifest_path=manifest_path,
        selected_path=selected_path,
        manifest=manifest,
        selected=selected,
        target=target,
        selected_species=selected_n,
    )


__all__ = ["ALLOWED_SOURCES", "CapacityHandoff", "resolve_capacity_handoff"]
