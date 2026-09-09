"""Preserve the real closed recovery outcome and exercise its pre-inference stop.

Artificial status rows below test software only, not ecological frequencies.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from scripts.analysis.audit_rgfca_background_recovery_completion import (
    EXECUTION, RESULT_COMMIT, RESULT_PATH, ROOT, git_bytes,
)

FINALIZER = "scripts/analysis/finalize_global_rgfca_background_control.py"
FINALIZER_SHA = "89b6659a2eccd85c85b0b46c44d44504cc51eb431c52eac36825614ab9508e35"


def test_closed_receipt_matches_immutable_result():
    result = json.loads(git_bytes(RESULT_COMMIT, RESULT_PATH))
    audit = json.loads((ROOT / "docs/supporting/rgfca_background_recovery_completion_audit_v1.json").read_text())
    assert result["status"] == "not_evaluable_incomplete_exact_background_recovery"
    assert result["exact_rows"] == audit["exact_rows"] == 21339
    assert result["failed_rows"] == audit["failed_rows"] == 85
    assert audit["rows_verified"] == result["matched_rows_required"] == 21424
    assert audit["partitions_verified"] == 128 and audit["species_in_fixed_frame"] == 369
    assert audit["lineage"] == result["lineage"]
    assert audit["failure_combinations"] == {
        "flower_mask_pixels_exact + flower_palette_exact": 79,
        "background_pixels_exact": 6,
    }
    assert result["primary_statistic_computed"] is False
    assert result["replacement_photos_used"] is False
    assert result["denominator_adapted_after_recovery"] is False
    assert audit["reserve_outcomes_read"] is False
    assert audit["raw_images_reacquired_by_audit"] is False


def real_finalizer(monkeypatch, tmp_path, failed_rows=85, missing=False, duplicate=False):
    raw = git_bytes(EXECUTION, FINALIZER)
    assert hashlib.sha256(raw).hexdigest() == FINALIZER_SHA
    tree = ast.parse(raw)
    # Only the real main and byte-hash helper; no SciPy or image-model import.
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)
             and n.name in {"main", "sha256_file"}]
    source_contract = json.loads(git_bytes(EXECUTION, "docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json"))
    measured = tmp_path / "synthetic_measurement.txt"
    measured.write_bytes(b"synthetic non-colour sentinel")
    source_contract["matched_frame"]["measured_table_sha256"] = hashlib.sha256(measured.read_bytes()).hexdigest()
    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps(source_contract), encoding="utf-8")
    amendment = tmp_path / "amendment.json"
    amendment.write_text(json.dumps({"status": "frozen_after_flower_only_omnibus_outcome_but_before_any_background_colour_recovery"}), encoding="utf-8")
    recovery = tmp_path / "recovery"
    recovery.mkdir()
    for shard in range(128 - int(missing)):
        (recovery / f"recovery_s{shard // 4:02d}_p{shard % 4:02d}.csv").touch()

    def read_statuses(path, **kwargs):
        assert path != measured, "Non-evaluable branch attempted to open the statistical frame"
        name = path.stem
        shard = int(name.split("_s")[1].split("_p")[0]) * 4 + int(name.split("_p")[1])
        indices = list(range(shard, 21424, 128))
        return pd.DataFrame({
            "measurement_id": ["synthetic-duplicate" if duplicate else f"synthetic-{i}" for i in indices],
            "recovery_status": ["original_roi_or_flower_palette_not_exactly_reproduced" if i < failed_rows
                                else "exact_matched_background_recovered" for i in indices],
        })

    namespace = {
        "argparse": argparse, "hashlib": hashlib, "json": json, "Path": Path,
        "pd": SimpleNamespace(read_csv=read_statuses, concat=pd.concat),
        "DEFAULT_CONTRACT": contract, "DEFAULT_AMENDMENT": amendment, "DEFAULT_MEASURED": measured,
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), FINALIZER, "exec"), namespace)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [FINALIZER, "--recovery-dir", str(recovery), "--output-dir", str(out)])
    return namespace["main"], out


@pytest.mark.parametrize("failed_rows", [1, 85])
def test_real_finalizer_stops_before_frame_or_inference(monkeypatch, tmp_path, failed_rows):
    run, out = real_finalizer(monkeypatch, tmp_path, failed_rows=failed_rows)
    assert run() == 0  # A completed workflow is compatible with scientific non-evaluability.
    summary = json.loads((out / "summary.json").read_text())
    assert summary["failed_rows"] == failed_rows
    assert summary["exact_rows"] == 21424 - failed_rows
    assert summary["primary_statistic_computed"] is False
    assert summary["replacement_photos_used"] is False
    assert summary["denominator_adapted_after_recovery"] is False
    assert {p.name for p in out.iterdir()} == {"summary.json", "background_recovery_rows_v2.csv"}


@pytest.mark.parametrize("kwargs,match", [({"missing": True}, "expected 128"), ({"duplicate": True}, "unique matched photos")])
def test_real_finalizer_rejects_missing_or_duplicate_census(monkeypatch, tmp_path, kwargs, match):
    run, out = real_finalizer(monkeypatch, tmp_path, **kwargs)
    with pytest.raises(RuntimeError, match=match):
        run()
    assert not out.exists()


def test_audit_cli_help_does_not_import_models_or_read_reserve():
    completed = subprocess.run([sys.executable, str(ROOT / "scripts/analysis/audit_rgfca_background_recovery_completion.py"), "--help"],
                               cwd=ROOT, capture_output=True, text=True, check=True)
    assert "--partitions-dir" in completed.stdout


def test_manuscript_records_non_evaluability_not_background_null():
    manuscript = (ROOT / "docs/RGFCA_MANUSCRIPT.md").read_text(encoding="utf-8")
    assert "21,339" in manuscript and "85" in manuscript
    assert "not_evaluable_incomplete_exact_background_recovery" in manuscript
    assert "No background-adjusted statistic or p-value was computed" in manuscript
