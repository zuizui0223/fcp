import hashlib
import json
from pathlib import Path

import pytest

from scripts.data.freeze_jbi_atlas_roi_v4_training_evidence import (
    post_training_summary,
    validate_source_run,
)


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def make_completed_run(tmp_path: Path) -> tuple[Path, Path, dict[str, object]]:
    run = tmp_path / "run"
    frozen_train = run / "frozen_train"
    weights = frozen_train / "weights"
    frozen_train.mkdir(parents=True)
    weights.mkdir()
    weight = b"last-epoch-weight"
    result_rows = [
        f"{epoch},1.0,{0.8 if epoch == 50 else 0},{0.7 if epoch == 50 else 0},"
        f"{0.75 if epoch == 50 else 0},{0.5 if epoch == 50 else 0}"
        for epoch in range(1, 51)
    ]
    results = (
        "epoch,loss,metrics/precision(B),metrics/recall(B),metrics/mAP50(B),"
        "metrics/mAP50-95(B)\n" + "\n".join(result_rows) + "\n"
    ).encode()
    (run / "jrc_yolo11n_last_v4.pt").write_bytes(weight)
    (weights / "last.pt").write_bytes(weight)
    (frozen_train / "results.csv").write_bytes(results)
    (frozen_train / "args.yaml").write_text(
        "epochs: 50\nval: false\n", encoding="utf-8"
    )
    (tmp_path / "data.yaml").write_text(
        "train: images/train\n"
        "val: images/train  # schema-only; train script fixes val=False\n",
        encoding="utf-8",
    )
    write_json(
        run / "training_result_manifest.json",
        {
            "protocol": "jbi-atlas-roi-estimator-v4",
            "status": "complete_roi_v4_detector_training_not_yet_qualified",
            "epochs": 50,
            "weight_selection": "last epoch only; never best epoch",
            "trained_weight_sha256": digest(weight),
            "training_results_sha256": digest(results),
            "jrc_test_images_decoded_or_scored": False,
            "scaleout_candidate_pixels_opened": False,
            "environment": {"device": "cpu"},
        },
    )
    materialization_path = tmp_path / "training_manifest.json"
    materialization = {
        "protocol": "jbi-atlas-roi-estimator-v4",
        "status": "pass_roi_v4_training_materialization",
        "images": 400,
        "source_annotation_boxes": 6992,
        "evaluable_training_boxes": 6991,
        "source_not_evaluable_boxes": 1,
        "training_validation_enabled": False,
        "jrc_test_directory_read": False,
        "jrc_test_images_decoded_or_scored": False,
        "scaleout_candidate_pixels_opened": False,
    }
    write_json(materialization_path, materialization)
    contract = {"protocol": "jbi-atlas-roi-estimator-v4"}
    return run, materialization_path, contract


def test_completed_training_freezes_only_after_all_firewalls_close(
    tmp_path: Path,
) -> None:
    run, materialization_path, contract = make_completed_run(tmp_path)
    result, materialization, paths = validate_source_run(
        run, materialization_path, contract
    )
    assert result["epochs"] == 50
    assert materialization["images"] == 400
    assert set(paths) == {
        "trained_weight",
        "training_results",
        "training_arguments",
        "training_data",
        "training_result",
        "materialization",
    }


def test_training_summary_records_but_never_uses_final_train_metrics(
    tmp_path: Path,
) -> None:
    run, _, _ = make_completed_run(tmp_path)
    summary = post_training_summary(run / "frozen_train/results.csv")
    assert summary["automatic_post_training_training_set_summary_emitted"] is True
    assert summary["summary_epoch"] == 50
    assert summary["summary_metrics_used_for_weight_selection"] is False
    assert summary["best_pt_frozen_or_used"] is False


def test_training_evidence_rejects_any_test_access(tmp_path: Path) -> None:
    run, materialization_path, contract = make_completed_run(tmp_path)
    manifest = json.loads(materialization_path.read_text(encoding="utf-8"))
    manifest["jrc_test_directory_read"] = True
    write_json(materialization_path, manifest)
    with pytest.raises(RuntimeError, match="materialization changed"):
        validate_source_run(run, materialization_path, contract)
