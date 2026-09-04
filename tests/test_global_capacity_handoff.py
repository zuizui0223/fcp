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


def _auth(source: str, manifest: str, selected: str, target: int, n: int) -> dict[str, object]:
    return {
        "status": "authorize_exactly_one_metadata_only_global_candidate_acquisition",
        "branch": handoff.BRANCH,
        "capacity_source": source,
        "capacity_manifest_path": manifest,
        "selected_species_path": selected,
        "capacity_manifest_blob_sha": "blob",
        "selected_species_blob_sha": "blob",
        "selected_raw_photo_target": target,
        "selected_species": n,
        "candidate_image_pixels_may_open": False,
        "flower_colour_may_open": False,
        "target_relaxation_allowed": False,
        "additional_pages_after_result_allowed": False,
        "rerun_for_favourable_species_set_allowed": False,
    }


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
    _write(tmp_path, spec["manifest"], json.dumps(manifest))
    selected_path = tmp_path / spec["selected"]
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    _selected(300, 100).to_csv(selected_path, index=False)
    auth_path = _write(tmp_path, "auth.json", json.dumps(_auth("v2_primary", spec["manifest"], spec["selected"], 100, 300)))
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    result = handoff.resolve_capacity_handoff(tmp_path, auth_path)
    assert result.source == "v2_primary"
    assert result.target == 100
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
    _write(tmp_path, spec["manifest"], json.dumps(manifest))
    selected_path = tmp_path / spec["selected"]
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    _selected(350, 80).to_csv(selected_path, index=False)
    auth_path = _write(tmp_path, "auth.json", json.dumps(_auth("v3_transport_recovery", spec["manifest"], spec["selected"], 80, 350)))
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    result = handoff.resolve_capacity_handoff(tmp_path, auth_path)
    assert result.source == "v3_transport_recovery"
    assert result.target == 80
    assert result.selected_species == 350


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
    _write(tmp_path, spec["manifest"], json.dumps(manifest))
    selected_path = tmp_path / spec["selected"]
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    _selected(300, 100).to_csv(selected_path, index=False)
    bad = _auth("v2_primary", "docs/supporting/wrong.json", spec["selected"], 100, 300)
    auth_path = _write(tmp_path, "auth.json", json.dumps(bad))
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
    _write(tmp_path, spec["manifest"], json.dumps(manifest))
    selected_path = tmp_path / spec["selected"]
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    _selected(300, 100).to_csv(selected_path, index=False)
    auth_path = _write(tmp_path, "auth.json", json.dumps(_auth("v3_transport_recovery", spec["manifest"], spec["selected"], 100, 300)))
    monkeypatch.setattr(handoff, "git_blob_sha", lambda root, path: "blob")
    with pytest.raises(RuntimeError, match="failed-v2"):
        handoff.resolve_capacity_handoff(tmp_path, auth_path)
