from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.analysis import monarda_execution_guard as guard
from scripts.analysis import run_rgfca_monarda_region_agreement as runner


def census():
    return [{"split": "train", "coco_image_id": i, "member_path": f"artificial/{i}.jpg",
             "annotation_count": int(i != 69)} for i in range(110)]


def test_original_contract_is_preserved_and_amendment_is_pre_outcome():
    amendment = json.loads((guard.ROOT / guard.AMENDMENT).read_text(encoding="utf-8"))
    assert guard.blob_id((guard.ROOT / guard.PARENT).read_bytes(), text=True) == "b18286a286cca46c370a98d520b44372345fafaa"
    assert amendment["local_pixels_opened_before_amendment"] is False
    assert amendment["scope_clarification"]["incidental_internal_cielab_computation_permitted"] is True
    assert amendment["scope_clarification"]["continuous_colour_values_saved_or_analyzed"] is False
    assert amendment["scope_clarification"]["interface"].endswith(".measure(image)")
    assert amendment["canonical_output_directory"] == ".artifacts/monarda-region-agreement-v1"


def test_text_identity_is_explicit_not_applied_to_binary():
    assert guard.blob_id(b"a\r\nb\r\n", text=True) == guard.blob_id(b"a\nb\n")
    assert guard.blob_id(b"a\r\nb\r\n") != guard.blob_id(b"a\nb\n")


@pytest.mark.parametrize("key", ["python", "system", "machine", "packages"])
def test_environment_drift_is_rejected(key):
    expected = {"python": "3.10.11", "system": "Windows", "machine": "AMD64",
                "packages": {"Pillow": "12.2.0"}}
    observed = {**expected, "packages": {"pillow": "12.2.0"}}
    guard.check_environment(expected, observed)
    observed[key] = {"pillow": "12.0.0"} if key == "packages" else "changed"
    with pytest.raises(RuntimeError, match="mismatch"):
        guard.check_environment(expected, observed)


def test_missing_authorization_blocks_before_archive_or_model(tmp_path):
    with pytest.raises(RuntimeError, match="pixels remain closed"):
        guard.verify_authorization(tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_committed_authorization_binds_source_and_rejects_tampering(tmp_path, monkeypatch):
    source_bytes = {p: (guard.ROOT / p).read_bytes().replace(b"\r\n", b"\n")
                    for p in guard.QUALIFIED_PATHS}
    for p, raw in source_bytes.items():
        destination = tmp_path / p
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(raw)
    auth = {"status": "authorized_once_after_pre_pixel_qualification",
            "qualification": {"head_sha": "artificial-qualified-commit", "run_id": 1, "conclusion": "success"},
            "qualified_source_blobs": {p: guard.blob_id(b) for p, b in source_bytes.items()}}
    raw_auth = (json.dumps(auth) + "\n").encode()
    (tmp_path / guard.AUTHORIZATION).write_bytes(raw_auth)
    monkeypatch.setattr(guard, "git_bytes", lambda rev, path, root: raw_auth if path == guard.AUTHORIZATION else source_bytes[path])
    monkeypatch.setattr(guard.subprocess, "run", lambda *a, **kw: None)
    expected = json.loads(source_bytes[guard.AMENDMENT])["execution_environment"]
    monkeypatch.setattr(guard, "environment_snapshot", lambda: expected)
    verified, _, receipt = guard.verify_authorization(tmp_path)
    assert verified == auth and receipt["monarda_pixels_decoded"] is False
    target = tmp_path / "scripts/analysis/run_rgfca_monarda_region_agreement.py"
    target.write_bytes(target.read_bytes() + b"# changed after qualification\n")
    with pytest.raises(RuntimeError, match="source identity mismatch"):
        guard.verify_authorization(tmp_path)


def test_ledger_retains_110_pending_including_unknown(tmp_path):
    output = tmp_path / "once"
    ledger = guard.ExecutionLedger(output, census(), {"model_loaded": False})
    rows = json.loads((output / "all_110_checkpoint.json").read_text())
    assert len(rows) == 110
    assert sum(r["annotation_count"] == 0 for r in rows) == 1
    assert all(r["scoring_status"] == "not_started" for r in rows)
    assert "segmentations" not in (output / "all_110_checkpoint.json").read_text()
    ledger.event("artificial_start")
    ledger.fail("artificial_setup_failure")
    terminal = json.loads((output / "execution_terminal.json").read_text())
    assert terminal["all_110_images_accounted_for"] is True
    assert terminal["terminal_scoring_rows"] == 0
    assert terminal["status"] == "not_evaluable_execution_incomplete"


def test_existing_output_rejects_rerun_without_overwrite(tmp_path):
    output = tmp_path / "once"
    guard.ExecutionLedger(output, census(), {})
    before = (output / "execution_start.json").read_bytes()
    with pytest.raises(FileExistsError):
        guard.ExecutionLedger(output, census(), {})
    assert (output / "execution_start.json").read_bytes() == before


def test_duplicate_or_partial_census_cannot_claim_execution(tmp_path):
    for rows in (census()[:-1], census()[:-1] + [census()[0]]):
        with pytest.raises(RuntimeError, match="110 unique"):
            guard.ExecutionLedger(tmp_path / "never", rows, {})
    assert not (tmp_path / "never").exists()


def test_incremental_score_persisted_and_duplicate_rejected(tmp_path):
    ledger = guard.ExecutionLedger(tmp_path / "once", census(), {})
    row = {"split": "train", "coco_image_id": 0, "model_status": "model_runtime_failure:Artificial"}
    ledger.update("score", row)
    with pytest.raises(RuntimeError, match="duplicate image score"):
        ledger.update("score", row)
    ledger.fail("artificial_interruption")
    rows = json.loads((ledger.output / "all_110_checkpoint.json").read_text())
    assert len(rows) == 110 and rows[0]["scoring_status"] == "terminal"
    assert sum(r["scoring_status"] == "not_started" for r in rows) == 109


def test_runtime_setup_precedes_alignment_and_exception_leaves_no_pixels(tmp_path, monkeypatch):
    ledger = guard.ExecutionLedger(tmp_path / "once", census(), {})
    calls = []
    def fail_import(cache):
        calls.append("runtime_setup")
        raise RuntimeError("artificial missing runtime")
    monkeypatch.setattr(runner, "import_runtime", fail_import)
    monkeypatch.setattr(runner, "_alignment_audit", lambda *a: calls.append("pixels_opened"))
    args = SimpleNamespace(output_dir=ledger.output)
    with pytest.raises(RuntimeError, match="artificial missing runtime"):
        runner._execute(args, {}, {}, tmp_path, ledger, census())
    assert calls == ["runtime_setup"]


def test_isolated_cache_rejects_changed_runtime_without_overwrite(tmp_path, monkeypatch):
    raw = b"# artificial source only\n"
    monkeypatch.setattr(guard, "RUNTIME_BLOBS", {"fcp_pipeline/artificial.py": guard.blob_id(raw)})
    monkeypatch.setattr(guard, "git_bytes", lambda *a: raw)
    cache, receipt = guard.materialize_runtime(tmp_path)
    assert receipt["fcp_pipeline/artificial.py"]["sha256"] == guard.sha_bytes(raw)
    target = cache / "fcp_pipeline/artificial.py"
    target.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="existing immutable cache differs"):
        guard.materialize_runtime(tmp_path)
    assert target.read_bytes() == b"changed"


def test_complete_artificial_execution_retains_failure_unknown_and_no_colour(tmp_path, monkeypatch):
    # Entirely artificial 8x8 JPEGs; neither target images nor a real model.
    from io import BytesIO
    import zipfile
    import numpy as np
    from PIL import Image

    data = BytesIO()
    Image.new("RGB", (8, 8), (100, 100, 100)).save(data, format="JPEG")
    rows = census()
    for row in rows:
        row.update(width_declared=8, height_declared=8,
                   segmentations=[[0, 0, 7, 0, 7, 7, 0, 7]] if row["annotation_count"] else [])
        row["split"] = "train" if row["coco_image_id"] < 77 else "valid" if row["coco_image_id"] < 99 else "test"
    archive = tmp_path / "artificial.zip"
    with zipfile.ZipFile(archive, "w") as z:
        for row in rows:
            z.writestr(row["member_path"], data.getvalue())
    empty_json = tmp_path / "artificial.json"
    empty_json.write_text("{}")
    class FakeEstimator:
        def __init__(self, *a, **kw):
            self.calls = 0
        def measure(self, image):
            self.calls += 1
            if self.calls == 1:
                raise RuntimeError("artificial whole-measure failure")
            return {"flower_mask": np.ones((8, 8), dtype=bool), "retained_instances": 1,
                    "automated_colour_state_status": "artificial_status",
                    "unused_internal_colour_value": "must_not_be_serialized"}
    runtime = SimpleNamespace(
        validate_roi_v4_contract=lambda c: None,
        file_sha256=lambda p: "artificial_weight_hash",
        validate_scaleout_authorization=lambda *a, **kw: None,
        FrozenFlowerColourEstimator=FakeEstimator,
    )
    monkeypatch.setattr(runner, "import_runtime", lambda cache: runtime)
    ledger = guard.ExecutionLedger(tmp_path / "once", rows, {"artificial": True})
    args = SimpleNamespace(output_dir=ledger.output, archive=archive,
        roi_contract=empty_json, locked_result=empty_json, detector_weight=empty_json,
        efficient_sam_dir=tmp_path, torch_threads=2, contract=empty_json,
        intake=empty_json, mapping=empty_json)
    contract = json.loads((guard.ROOT / guard.PARENT).read_text(encoding="utf-8"))
    contract["parent_checkpoint"]["archive_sha256"] = runner.sha256_file(archive)
    assert runner._execute(args, contract, {"protocol": "artificial"}, tmp_path, ledger, rows) == 0
    result = json.loads((ledger.output / "monarda_region_agreement_result_v1.json").read_text())
    checkpoint = json.loads((ledger.output / "all_110_checkpoint.json").read_text())
    assert result["model_runtime_failure_images"] == 1
    assert result["aggregate"]["reference_pixels_total"] == 109 * 64
    assert result["aggregate"]["intersection_pixels_total"] == 108 * 64
    assert result["unknown_zero_annotation_image"]["counted_as_verified_negative"] is False
    assert len(checkpoint) == 110 and all(r["scoring_status"] == "terminal" for r in checkpoint)
    assert result["continuous_colour_values_saved_or_analyzed"] is False
    for path in ledger.output.iterdir():
        assert "must_not_be_serialized" not in path.read_text()
