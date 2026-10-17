"""No-network fail-closed tests for multi-invocation durable 93M photo-ID source stages."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.acquisition import run_fcp_durable_source_phases as mod


def prepared(tmp_path, monkeypatch):
    root = tmp_path / "durable"
    (root / "inputs").mkdir(parents=True)
    (root / "references").mkdir()
    (root / "snapshots").mkdir()
    (root / "FCP_DURABLE_STORAGE_APPROVED").write_text(
        "I confirm this is persistent FCP source storage, not hosted-runner scratch.")
    frame = root / "inputs/42111_species.csv"
    counts = root / "references/all_42111_species_photo_opportunity_counts.csv"
    audit = root / "references/all_42111_photo_opportunity_audit.json"
    frame.write_text("taxon test frame")
    counts.write_text("photo count test input")
    audit.write_text(json.dumps({
        "status": "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
        "snapshot_date": mod.DATE,
        "source_frame_species": mod.SPECIES,
        "source_frame_sha256": mod.FRAME_SHA,
        "source_archive_sha256": mod.ARCHIVE_SHA,
        "count_ledger_sha256": mod.COUNTS_SHA,
        "all_source_photo_opportunity_ids": mod.SOURCE_PHOTOS,
        "original_source_taxon_IDs_and_names_independently_matched": True,
        "all_source_photo_id_manifests_frozen": False,
        "all_original_image_pixels_measured": False,
    }))
    # A sparse placeholder is NOT a real archive and is never opened by
    # tests: fake source hashes only inside tests, without provider traffic.
    archive = root / "snapshots" / mod.SOURCE_ARCHIVE
    archive.touch()
    with archive.open("r+b") as stream:
        stream.truncate(35_966_549_880)
    repository = tmp_path / "repo"
    script = repository / "scripts/acquisition/execute_fcp_42111_census.py"
    script.parent.mkdir(parents=True)
    script.write_text("synthetic, never launched")
    monkeypatch.setattr(mod, "validate_root", lambda p, runner_temp=None: root)
    monkeypatch.setattr(mod, "sha256", lambda p: {
        frame.name: mod.FRAME_SHA,
        counts.name: mod.COUNTS_SHA,
        audit.name: mod.AUDIT_SHA,
        archive.name: mod.ARCHIVE_SHA,
    }[p.name])
    monkeypatch.setattr(mod.shutil, "disk_usage",
                        lambda p: SimpleNamespace(free=350_000_000_000))
    return root, repository, audit


def test_rejects_ephemeral_and_relative_root(tmp_path):
    with pytest.raises(ValueError, match="absolute"):
        mod.validate_root(Path("relative/source/root"))
    with pytest.raises(ValueError, match="Ephemeral/unsafe"):
        mod.validate_root(tmp_path)


def test_shard_budget_limited_per_invocation_but_not_study():
    for start, stop in [(0, 5), (127, 129), (6, 6), (-1, 1)]:
        with pytest.raises(ValueError, match="Shard range|four"):
            mod.validate_shards(start, stop)
    mod.validate_shards(0, 4)
    mod.validate_shards(124, 128)


def test_source_preflight_checks_exact_independent_real_counts(tmp_path, monkeypatch):
    root, repo, _ = prepared(tmp_path, monkeypatch)
    got = mod.run_phase(root, "preflight", repo=repo)
    assert got["source_photo_opportunities"] == 93_209_351
    assert got["original_photo_pixels_measured"] is False
    assert got["flower_colour_W_estimated"] is False
    assert (root / "status/phase_preflight_000_004.json").is_file()


def test_incompatible_count_audit_blocks_photo_id_source(tmp_path, monkeypatch):
    root, repo, audit = prepared(tmp_path, monkeypatch)
    source = json.loads(audit.read_text())
    source["all_source_photo_opportunity_ids"] -= 1
    audit.write_text(json.dumps(source))
    with pytest.raises(ValueError, match="Source-reference audit"):
        mod.run_phase(root, "preflight", repo=repo)


def test_stage_requires_zero_existing_work_and_exact_photo_id_denominator(tmp_path, monkeypatch):
    root, repo, _ = prepared(tmp_path, monkeypatch)
    calls = []
    receipt = root / "source_20260927/full_source_execution_receipt.json"

    def stage_subprocess(command, *, check, cwd):
        calls.append(command)
        receipt.parent.mkdir(parents=True, exist_ok=True)
        receipt.write_text(json.dumps({
            "status": mod.SOURCE_STAGED,
            "staged_original_source_photo_ids": mod.SOURCE_PHOTOS,
            "complete_metadata_census": False,
        }))
    monkeypatch.setattr(mod.subprocess, "run", stage_subprocess)
    got = mod.run_phase(root, "stage", repo=repo)
    assert got["status"] == "ORIGINAL_SOURCE_STAGE_ONLY_NO_COLOURS"
    assert len(calls) == 1
    assert "--stage-only" in calls[0]
    assert "--release-verified-temporary-tables" in calls[0]
    assert "--archive" in calls[0]
    assert calls[0][1:3] == ["-m", "scripts.acquisition.execute_fcp_42111_census"]
    with pytest.raises(ValueError, match="inspection"):
        mod.run_phase(root, "stage", repo=repo)


def test_emit_requires_staged_source_and_is_four_shard_bounded(tmp_path, monkeypatch):
    root, repo, _ = prepared(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="No completely staged"):
        mod.run_phase(root, "emit", 12, 16, repo=repo)
    receipt = root / "source_20260927/full_source_execution_receipt.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps({"status": mod.SOURCE_STAGED}))
    calls = []
    monkeypatch.setattr(mod.subprocess, "run",
                        lambda cmd, *, check, cwd: calls.append(cmd))
    report = mod.run_phase(root, "emit", 12, 16, repo=repo)
    assert report["phase"] == "emit"
    assert report["start"] == 12 and report["stop"] == 16
    assert "--shard-stop" in calls[0] and "16" in calls[0]
    assert calls[0][1:3] == ["-m", "scripts.acquisition.emit_fcp_staged_source_shards"]
    assert report["original_photo_pixels_measured"] is False


def test_finalization_requires_independent_93m_photo_id_count_admission(tmp_path, monkeypatch):
    root, repo, _ = prepared(tmp_path, monkeypatch)
    receipt = root / "source_20260927/full_source_execution_receipt.json"
    receipt.parent.mkdir(parents=True)
    receipt.write_text(json.dumps({"status": mod.SOURCE_STAGED}))
    calls = []

    def faux_process(cmd, *, check, cwd):
        calls.append(cmd)
        if "scripts.acquisition.verify_fcp_complete_photo_ids_vs_original_counts" in " ".join(cmd):
            admission = root / "source_20260927/independent_original_all_ID_source_admission.json"
            admission.write_text(json.dumps({
                "status": mod.ADMITTED,
                "distinct_original_photo_ids": mod.SOURCE_PHOTOS,
                "full_study_original_pixels_classified": False,
            }))
    monkeypatch.setattr(mod.subprocess, "run", faux_process)
    got = mod.run_phase(root, "finalize", repo=repo)
    assert len(calls) == 2
    assert "--skip-package" in calls[0]
    assert calls[0][1:3] == ["-m", "scripts.acquisition.finalize_fcp_staged_source"]
    assert calls[1][1:3] == ["-m", "scripts.acquisition.verify_fcp_complete_photo_ids_vs_original_counts"]
    assert got["full_source_photo_ID_universe_frozen"] is True
    assert got["original_photo_pixels_measured"] is False

    # False claim of a complete original image measurement cannot pass.
    def false_admission(cmd, *, check, cwd):
        if "scripts.acquisition.verify_fcp_complete_photo_ids_vs_original_counts" in " ".join(cmd):
            admission = root / "source_20260927/independent_original_all_ID_source_admission.json"
            admission.write_text(json.dumps({
                "status": mod.ADMITTED,
                "distinct_original_photo_ids": mod.SOURCE_PHOTOS,
                "full_study_original_pixels_classified": True,
            }))
    monkeypatch.setattr(mod.subprocess, "run", false_admission)
    with pytest.raises(ValueError, match="photo-ID admission"):
        mod.run_phase(root, "finalize", repo=repo)
