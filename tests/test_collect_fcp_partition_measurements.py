"""Synthetic-only regression for safe collection of parallel partition outputs."""
from __future__ import annotations

import json

import pytest

from fcp_pipeline.fcp_fresh_study import reassemble_fresh
from scripts.acquisition.collect_fcp_partition_measurements import collect
from test_execute_fcp_full_photo_work_queue import sample, patch_full_method, call
from test_fcp_fresh_study import attestation


def test_partition_workers_are_disjoint_and_reassembled_only_after_all_done(tmp_path, monkeypatch):
    f = sample(tmp_path)
    patch_full_method(monkeypatch, tmp_path)
    root = tmp_path / "partition_outputs"
    p0 = root / "partition_000"
    p1 = root / "partition_001"
    a = call(f, p0, budget=0, partition=0)
    assert a["full_study_original_photo_measurement_complete"] is False
    partial = collect(f, root, tmp_path / "joined", batch_size=2, _test_only_small_source=True)
    assert partial["n_pending_work_units"] > 0
    assert partial["complete_source_photo_terminal_coverage"] is False
    call(f, p1, budget=0, partition=1)
    complete = collect(f, root, tmp_path / "joined", batch_size=2, _test_only_small_source=True)
    assert complete["n_pending_work_units"] == 0
    assert complete["complete_source_photo_terminal_coverage"] is True
    assert complete["W_H2_environment_computed"] is False
    assert complete["source_photo_ids"] == 8
    sealed = tmp_path / "joined" / "sealed_terminal_batches"
    result = reassemble_fresh(f, sealed, attestation(tmp_path), tmp_path / "assembled")
    assert result["status"] == "COMPLETE_FRESH_MEASURED_PHOTO_CENSUS"
    assert result["n_terminal_measurements"] == 8
    again = collect(f, root, tmp_path / "joined", batch_size=2, _test_only_small_source=True)
    assert again["n_completed_work_units"] == complete["n_completed_work_units"]


def test_tampered_collected_bytes_cannot_be_replaced(tmp_path, monkeypatch):
    f = sample(tmp_path)
    patch_full_method(monkeypatch, tmp_path)
    root = tmp_path / "partition_outputs"
    call(f, root / "partition_000", budget=0, partition=0)
    joined = tmp_path / "joined"
    collect(f, root, joined, batch_size=2, _test_only_small_source=True)
    target = next((joined / "sealed_terminal_batches").glob("terminal_*.csv"))
    target.write_text("incorrect historic or corrupted output")
    with pytest.raises(ValueError, match="Previously collected"):
        collect(f, root, joined, batch_size=2, _test_only_small_source=True)


def test_wrong_partition_terminal_is_refused(tmp_path, monkeypatch):
    f = sample(tmp_path)
    patch_full_method(monkeypatch, tmp_path)
    root = tmp_path / "partition_outputs"
    call(f, root / "partition_000", budget=0, partition=0)
    folder = root / "partition_000" / "sealed_terminal_batches"
    (folder / "terminal_p001_b0000000.csv").write_text("wrong partition")
    with pytest.raises(ValueError, match="wrong-partition"):
        collect(f, root, tmp_path / "joined", batch_size=2, _test_only_small_source=True)


def test_cannot_relabel_small_fixture_as_42111_species(tmp_path):
    f = sample(tmp_path)
    with pytest.raises(ValueError, match="42111"):
        collect(f, tmp_path / "empty", tmp_path / "joined")
