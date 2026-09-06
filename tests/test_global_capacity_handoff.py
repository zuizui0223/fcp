from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import fcp_pipeline.global_capacity_handoff as handoff


def _write(root: Path, rel: str, text: str) -> Path:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _selected(n: int, target: int) -> pd.DataFrame:
    return pd.DataFrame({
        "species": [f"sp_{i}" for i in range(n)],
        "inat_taxon_id": range(1000, 1000 + n),
        "after_observer_cap": [target] * n,
        "maximum_span_km": [1000.0] * n,
        "selected_raw_photo_target": [target] * n,
    })


def _budget(root: Path, maximum: int = 1000, seed: int = 20260916) -> Path:
    return _write(
        root,
        handoff.BUDGET_PATH,
        json.dumps({
            "status": "frozen_after_metadata_capacity_scale_observed_before_candidate_acquisition_outcome_and_before_any_global_candidate_pixels",
            "candidate_species_budget": {"maximum_species": maximum, "selection_seed": seed},
        }),
    )


def _auth(
    source: str,
    manifest: str,
    selected: str,
    target: int,
    capacity_n: int,
    bounded_n: int,
    *,
    maximum: int = 1000,
    seed: int = 20260916,
    digest: str = "DIGEST",
) -> dict[str, object]:
    return {
        "status": "authorize_exactly_one_metadata_only_global_candidate_acquisition",
        "branch": handoff.BRANCH,
        "capacity_source": source,
        "capacity_manifest_path": manifest,
        "selected_species_path": selected,
        "capacity_manifest_blob_sha": "blob",
        "selected_species_blob_sha": "blob",
        "candidate_species_budget_amendment_path": handoff.BUDGET_PATH,
        "candidate_species_budget_amendment_blob_sha": "blob",
        "selected_raw_photo_target": target,
        "capacity_selected_species": capacity_n,
        "candidate_species_budget": maximum,
        "candidate_species_selection_seed": seed,
        "selected_species": bounded_n,
        "candidate_taxon_id_sha256": digest,
        "candidate_image_pixels_may_open": False,
        "flower_colour_may_open": False,
        "target_relaxation_allowed": False,
        "additional_pages_after_result_allowed": False,
        "rerun_for_favourable_species_set_allowed": False,
    }


def _resolved_auth(root: Path, source: str, manifest: dict[str, object], n: int, target: int, *, maximum: int = 1000):
    spec = handoff.ALLOWED_SOURCES[source]
    _write(root, spec["manifest"], json.dumps(manifest))
    selected_path = root / spec["selected"]
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    frame = _selected(n, target)
    frame.to_csv(selected_path, index=False)
    _budget(root, maximum=maximum)
    taxa = handoff.select_hashed_taxa(frame["inat_taxon_id"], maximum_species=maximum, seed=20260916)
    bounded = frame.loc[frame["inat_taxon_id"].isin(set(taxa))]
    digest = handoff.taxon_digest(bounded["inat_taxon_id"])
    auth = _auth(source, spec["manifest"], spec["selected"], target, n, len(bounded), maximum=maximum, digest=digest)
    return _write(root, "auth.json", json.dumps(auth))


def test_resolve_normal_v2_handoff(tmp_path, monkeypatch):
    spec = handoff.ALLOWED_SOURCES["v2_primary"]
    manifest = {
        "status": spec["status"],
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "actual_image_acquisition_authorized": False,
        "selected_raw_photo_target": 100,
        "selected_species": 300,
        "minimum_metadata_eligible_species": 300,
    }
    auth_path = _resolved_auth(tmp_path, "v2_primary", manifest, 300, 100)
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    result = handoff.resolve_capacity_handoff(tmp_path, auth_path)
    assert result.source == "v2_primary"
    assert result.target == 100
    assert result.capacity_selected_species == 300
    assert result.selected_species == 300


def test_resolve_transport_recovered_v3_handoff(tmp_path, monkeypatch):
    spec = handoff.ALLOWED_SOURCES["v3_transport_recovery"]
    manifest = {
        "status": spec["status"],
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "actual_image_acquisition_authorized": False,
        "selected_raw_photo_target": 80,
        "selected_species": 350,
        "minimum_metadata_eligible_species": 300,
        "original_v2_status": "not_evaluable_capacity_scan_due_request_failure",
        "second_recovery_permitted": False,
        "biological_rules_changed": False,
    }
    auth_path = _resolved_auth(tmp_path, "v3_transport_recovery", manifest, 350, 80)
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    result = handoff.resolve_capacity_handoff(tmp_path, auth_path)
    assert result.source == "v3_transport_recovery"
    assert result.target == 80
    assert result.capacity_selected_species == 350
    assert result.selected_species == 350


def test_handoff_caps_large_capacity_at_1000_before_candidate_query(tmp_path, monkeypatch):
    spec = handoff.ALLOWED_SOURCES["v2_primary"]
    manifest = {
        "status": spec["status"],
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "actual_image_acquisition_authorized": False,
        "selected_raw_photo_target": 100,
        "selected_species": 2400,
        "minimum_metadata_eligible_species": 300,
    }
    auth_path = _resolved_auth(tmp_path, "v2_primary", manifest, 2400, 100, maximum=1000)
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    result = handoff.resolve_capacity_handoff(tmp_path, auth_path)
    assert result.capacity_selected_species == 2400
    assert result.selected_species == 1000
    assert len(result.selected) == 1000


def test_handoff_rejects_source_path_mismatch(tmp_path, monkeypatch):
    spec = handoff.ALLOWED_SOURCES["v2_primary"]
    manifest = {
        "status": spec["status"],
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "actual_image_acquisition_authorized": False,
        "selected_raw_photo_target": 100,
        "selected_species": 300,
    }
    auth_path = _resolved_auth(tmp_path, "v2_primary", manifest, 300, 100)
    bad = json.loads(auth_path.read_text())
    bad["capacity_manifest_path"] = "docs/supporting/wrong.json"
    auth_path.write_text(json.dumps(bad))
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    with pytest.raises(RuntimeError, match="paths"):
        handoff.resolve_capacity_handoff(tmp_path, auth_path)


def test_v3_handoff_rejects_missing_failed_v2_lineage(tmp_path, monkeypatch):
    spec = handoff.ALLOWED_SOURCES["v3_transport_recovery"]
    manifest = {
        "status": spec["status"],
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "actual_image_acquisition_authorized": False,
        "selected_raw_photo_target": 100,
        "selected_species": 300,
        "original_v2_status": "something_else",
        "second_recovery_permitted": False,
        "biological_rules_changed": False,
    }
    auth_path = _resolved_auth(tmp_path, "v3_transport_recovery", manifest, 300, 100)
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    with pytest.raises(RuntimeError, match="failed-v2"):
        handoff.resolve_capacity_handoff(tmp_path, auth_path)
