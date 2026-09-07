"""Reconstruct completed ROI qualification, never remeasure any photograph.

Read only literal immutable Git objects from the 100-image JRC locked test.
No network, model import, reserve result, photograph decode, or new inference.
Artificial masks below replay only the frozen region-pooling code with stubbed
detector/encoder/decoder responses; they do not test botanical accuracy.
"""
import ast
import csv
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import subprocess
from types import SimpleNamespace

import numpy as np
from PIL import Image
import pytest

ROOT = Path(__file__).resolve().parents[1]
QUALIFICATION = "2449c77597cc4c57f18eb0aaa446211143ce3be5"
RUNTIME_SOURCE = "9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
RUNTIME_PATH = "fcp_pipeline/flower_roi_v4_runtime.py"
RUNTIME_SHA256 = "e269930cd9a06a58503277ab0b105f3396aeabe597deebe89103a33611e1e9f9"
DIRECTORY = "data/atlas/qualification/roi_v4_locked_test/"
SOURCES = {
    "contract": ("docs/supporting/jbi_atlas_roi_estimator_contract_v4.json", "e3ffefdcd4aa6719a10b45c57886ee895f10f424a0d0b174a46541d92dd2f7a1"),
    "amendment": ("docs/supporting/jbi_atlas_roi_v4_reference_size_amendment_v1.json", "0c99f99aed483fb036a850a4b5b08eac3f88d7dc41595ec4110e8b9698fc93cd"),
    "manifest": (DIRECTORY + "gate_evidence_manifest.json", "dc467f9ac4e12325c86c5c4e305216b02adc6cc97b9fbebb67fd28236fbb89e8"),
    "result": (DIRECTORY + "jrc_roi_v4_locked_test_result.json", "420ae571e0950bfe252d5f0de34b4e047785c3efb2968d9f14abfc2161c08aea"),
    "rows": (DIRECTORY + "jrc_roi_v4_locked_test_rows.csv", "9f834de4bcdd530581d7713376cc918aa587e80726bba32021fa2dfa1034b248"),
    "runner": (DIRECTORY + "gate_executable.py", "f43048839f5dde8ec1acd56c1163d988a8a7cddb024a289a6a2ba96b76c84074"),
    "rules": ("fcp_pipeline/flower_roi_v4.py", "9c965e4ef5167c3d4ed0034657b3dbcc8399aef0d7feea1df1c8d7656d659b6b"),
}


@lru_cache(None)
def raw_source(key):
    path, expected = SOURCES[key]
    raw = subprocess.check_output(["git", "show", f"{QUALIFICATION}:{path}"], cwd=ROOT)
    assert hashlib.sha256(raw).hexdigest() == expected, path
    return raw


def record(key):
    return json.loads(raw_source(key))


def rows():
    return list(csv.DictReader(io.StringIO(raw_source("rows").decode("utf-8"))))


def total(name):
    return sum(int(row[name]) for row in rows())


def reconstructed_metrics():
    all_rows = rows()
    tp, fp, fn = (total(k) for k in ("true_positive", "false_positive", "false_negative"))
    admitted = sum(row["estimator_admitted"] == "True" for row in all_rows)
    return {
        "images": len(all_rows),
        "admitted_images": admitted,
        "admitted_fraction": admitted / len(all_rows),
        "detector_precision_iou_0_5": tp / (tp + fp),
        "detector_recall_iou_0_5": tp / (tp + fn),
        "medium_reference_object_recall": total("medium_hit_boxes") / total("medium_reference_boxes"),
        "large_reference_object_recall": total("large_hit_boxes") / total("large_reference_boxes"),
        "pooled_mask_pixels_inside_reference_box_union": total("mask_pixels_inside_reference_box_union") / total("mask_pixels"),
        "median_image_mask_pixels_inside_reference_box_union": statistics.median(float(r["image_mask_pixels_inside_reference_box_union"]) for r in all_rows),
    }


def test_exact_completed_evidence_and_runtime_result_identity():
    for key in SOURCES:
        raw_source(key)
    manifest, result = record("manifest"), record("result")
    for key in ("result", "rows", "runner"):
        assert manifest["evidence"][key]["sha256_exact"] == SOURCES[key][1]
    assert result["rows_sha256"] == SOURCES["rows"][1]
    assert result["contract_sha256_lf_canonical_v1"] == SOURCES["contract"][1]
    assert result["reference_size_amendment_sha256_lf_canonical_v1"] == SOURCES["amendment"][1]
    runtime_result = subprocess.check_output(["git", "show", f"{RUNTIME_SOURCE}:{SOURCES['result'][0]}"], cwd=ROOT)
    assert runtime_result == raw_source("result")


def test_full_image_and_box_census_not_admitted_only_denominator():
    all_rows = rows()
    assert len(all_rows) == 100
    assert {int(r["image_id"]) for r in all_rows} == set(range(1, 101))
    for field in ("file_name", "image_sha256"):
        assert len({r[field] for r in all_rows}) == 100
    assert total("reference_boxes") == total("source_annotation_boxes") == 2524
    assert total("source_not_evaluable_boxes") == 0
    assert (total("true_positive"), total("false_positive"), total("false_negative")) == (2008, 741, 516)
    for r in all_rows:
        assert int(r["true_positive"]) + int(r["false_positive"]) == int(r["detector_predictions"])
        assert int(r["true_positive"]) + int(r["false_negative"]) == int(r["reference_boxes"])


def test_all_saved_metrics_reconstructed_without_model_execution():
    assert reconstructed_metrics() == pytest.approx(record("result")["metrics"], abs=1e-14)
    assert reconstructed_metrics()["admitted_fraction"] == .85


def test_all_eight_original_gate_decisions_including_equality():
    metrics = reconstructed_metrics()
    gates = record("contract")["gates"]["locked_test"]
    checks = {name: metrics[name.removeprefix("minimum_")] >= threshold
              for name, threshold in gates.items() if name.startswith("minimum_")}
    assert len(checks) == 8
    assert checks == record("result")["checks"]
    assert all(checks.values())
    assert metrics["large_reference_object_recall"] == gates["minimum_large_reference_object_recall"] == .5


def test_admission_and_all_fifteen_failures_retained():
    cfg = record("contract")["image_measurement"]
    failures = 0
    for row in rows():
        assert row["estimator_admitted"] in {"True", "False"}
        passed = (
            int(row["retained_instances"]) >= 1
            and int(row["mask_pixels"]) >= cfg["minimum_union_flower_pixels_on_original"]
            and int(row["background_pixels"]) >= cfg["minimum_background_pixels_on_original"]
            and float(row["horizontal_flip_mask_iou"]) >= cfg["horizontal_flip_mask_iou_minimum"]
            and float(row["horizontal_flip_colour_delta_e"]) <= cfg["horizontal_flip_colour_delta_e_maximum"]
        )
        assert passed == (row["estimator_admitted"] == "True")
        assert bool(row["failure_reasons"]) == (not passed)
        failures += not passed
    assert failures == 15


def test_reference_size_denominators_and_small_large_sample_ceiling():
    assert [total(f"{size}_reference_boxes") for size in ("small", "medium", "large")] == [2385, 135, 4]
    assert [total(f"{size}_hit_boxes") for size in ("small", "medium", "large")] == [1894, 112, 2]
    assert sum(int(r["large_reference_boxes"]) > 0 for r in rows()) == 3
    for row in rows():
        assert sum(int(row[f"{s}_reference_boxes"]) for s in ("small", "medium", "large")) == int(row["reference_boxes"])
        assert sum(int(row[f"{s}_hit_boxes"]) for s in ("small", "medium", "large")) == int(row["true_positive"])
    amendment = record("amendment")
    assert amendment["outcome_firewall"]["jrc_locked_test_images_decoded_or_scored"] is False
    assert amendment["parent_contract"]["sha256_lf_canonical_v1"] == SOURCES["contract"][1]


def test_mask_metric_is_box_containment_not_petal_intersection_over_union():
    assert total("mask_pixels_inside_reference_box_union") == 5543384
    assert total("mask_pixels") == 6448251
    for row in rows():
        numerator, denominator = (int(row[k]) for k in ("mask_pixels_inside_reference_box_union", "mask_pixels"))
        assert 0 <= numerator <= denominator
        assert float(row["image_mask_pixels_inside_reference_box_union"]) == pytest.approx(numerator / denominator)
    note = (ROOT / "docs/RGFCA_ROI_QUALIFICATION_AUDIT.md").read_text(encoding="utf-8")
    assert "not petal segmentation IoU" in note
    assert "2/4" in note and "three images" in note


def test_historical_failure_and_ecological_firewall_not_reinterpreted():
    contract, result = record("contract"), record("result")
    assert contract["immutable_v3_stop"]["manifest_status"] == "stop_jrc_development_failed"
    assert contract["immutable_v3_stop"]["admitted_images"] == 17
    assert result["status"] == "pass_roi_v4_locked_test"
    assert result["scaleout_candidate_pixels_opened"] is False
    assert "taxon-uniform error" in contract["claim_ceiling"]
    assert "calibrated reflectance" in contract["claim_ceiling"]


@lru_cache(None)
def frozen_runtime_tree():
    raw = subprocess.check_output(["git", "show", f"{RUNTIME_SOURCE}:{RUNTIME_PATH}"], cwd=ROOT)
    assert hashlib.sha256(raw).hexdigest() == RUNTIME_SHA256
    return ast.parse(raw.decode("utf-8"))


def test_measurement_interface_has_no_focal_taxon_input():
    cls = next(n for n in frozen_runtime_tree().body
               if isinstance(n, ast.ClassDef) and n.name == "FrozenFlowerColourEstimator")
    for name in ("measure", "_analyze_orientation"):
        method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name)
        assert [a.arg for a in method.args.args] == ["self", "image"]
        assert not method.args.kwonlyargs and method.args.kwarg is None
    assert record("contract")["detector"]["classes"] == ["flower"]


def replay_pooling(boxes):
    """Run the original pooling/helper functions, never model loading or measure."""
    tree = frozen_runtime_tree()
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef)
               and n.name == "FrozenFlowerColourEstimator")
    orientation = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                       and n.name == "_analyze_orientation")
    helper_names = {"_letterboxed_rgb", "_canvas_mask_to_original", "_background_annulus"}
    rule_names = {"letterbox_geometry", "box_to_canvas", "select_prompt_mask"}
    nodes = ast.parse("from __future__ import annotations").body
    nodes += [n for n in ast.parse(raw_source("rules").decode("utf-8")).body
              if isinstance(n, ast.FunctionDef) and n.name in rule_names]
    nodes += [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in helper_names]
    nodes += [orientation]
    namespace = {"np": np, "Image": Image, "math": math, "CANVAS_SIZE": 1024}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), RUNTIME_PATH, "exec"), namespace)

    def tensor(values):
        stub = SimpleNamespace()
        stub.detach = lambda: stub
        stub.cpu = lambda: stub
        stub.numpy = lambda: np.asarray(values)
        return stub

    predictions = SimpleNamespace(boxes=SimpleNamespace(
        xyxy=tensor(boxes), conf=tensor([.9] * len(boxes))))
    # The segmenter stub supplies positive logits everywhere; the real frozen
    # select_prompt_mask clips each instance to its distinct detector box.
    estimator = SimpleNamespace(
        contract=record("contract"),
        detector=SimpleNamespace(predict=lambda **kwargs: [predictions]),
        encoder=SimpleNamespace(run=lambda *args: [None]),
        decoder=SimpleNamespace(run=lambda *args: (
            np.ones((1, 1, 1, 1024, 1024), dtype=np.float32),
            np.ones((1, 1, 1), dtype=np.float32), None)),
    )
    return namespace["_analyze_orientation"](estimator, Image.new("RGB", (1024, 1024)))


@pytest.mark.parametrize("boxes,pixels", [
    ([[20, 20, 60, 60], [80, 20, 120, 60]], 3200),
    ([[20, 20, 60, 60], [40, 40, 80, 80]], 2800),
])
def test_all_retained_instances_are_pooled_not_one_focal_flower(boxes, pixels):
    result = replay_pooling(boxes)
    expected = np.zeros((1024, 1024), dtype=bool)
    for x0, y0, x1, y1 in boxes:
        expected[y0:y1, x0:x1] = True
    assert result["retained_instances"] == 2
    assert np.array_equal(result["flower_mask"], expected)
    assert result["flower_mask"].sum() == pixels
    assert not np.any(result["flower_mask"] & result["background_mask"])


def test_taxon_mask_limitation_is_explicit_in_manuscript_and_audit():
    manuscript = (ROOT / "docs/RGFCA_MANUSCRIPT.md").read_text(encoding="utf-8")
    audit = (ROOT / "docs/RGFCA_ROI_QUALIFICATION_AUDIT.md").read_text(encoding="utf-8")
    for text in (manuscript, audit):
        assert "observation's taxon label" in text
        assert "focal-species mask" in text
        assert "contamination rate" in text
    assert RUNTIME_SHA256 in audit
