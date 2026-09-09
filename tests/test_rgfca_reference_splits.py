"""Source-list audits need no photographs, masks or model runtime."""
import hashlib
import json
from pathlib import Path

import pytest

from scripts.analysis.audit_rgfca_reference_splits import (
    audit_bundle, parse_names, summarize_splits,
)


@pytest.mark.parametrize("raw", [
    b"", b"\n", b"a.JPG\n\n", b" a.JPG\n", b"a.JPG \n", b"../a.JPG\n",
    b"a/b.JPG\n", b"a\\b.JPG\n", b"a.JPG\x00\n", b"\xef\xbb\xbfa.JPG\n",
    b"a.JPG\vb.JPG\n", b"a.JPG\rb.JPG\n",
    b"not_a_filename\n", b"\xff.JPG\n",
    pytest.param(b"a" * 1_000_001, id="oversized-source"),
])
def test_invalid_or_ambiguous_rows_are_not_silently_repaired(raw):
    with pytest.raises(ValueError):
        parse_names(raw)


@pytest.mark.parametrize("raw", [b"A.JPG\nb.png", b"A.JPG\nb.png\n",
                                 b"A.JPG\r\nb.png\r\n"])
def test_line_endings_are_parsed_without_changing_source_identity(raw):
    assert parse_names(raw) == ["A.JPG", "b.png"]


def test_exact_overlap_and_duplicate_denominators_are_retained():
    train = ["a.JPG", "a.JPG", "b.JPG"]
    valid = ["a.JPG", "c.JPG"]
    out = summarize_splits(train, valid)
    assert train == ["a.JPG", "a.JPG", "b.JPG"]
    assert out["profiles"]["train"]["row_count"] == 3
    assert out["profiles"]["train"]["duplicate_extra_row_count"] == 1
    assert out["profiles"]["train"]["duplicate_filenames"] == {"a.JPG": 2}
    assert out["exact_union_filename_count"] == 3
    assert out["exact_shared_filenames"] == ["a.JPG"]
    assert not out["filename_integrity_gate_passed"]


def test_casefold_collision_is_not_hidden_by_exact_disjointness():
    out = summarize_splits(["A.JPG"], ["a.jpg"])
    assert out["filename_split_disjoint"]
    assert out["casefold_shared_keys"] == ["a.jpg"]
    assert out["casefold_spelling_collisions"] == {"a.jpg": ["A.JPG", "a.jpg"]}
    assert not out["filename_integrity_gate_passed"]


def test_within_split_casefold_collision_fails_gate():
    out = summarize_splits(["A.JPG", "a.jpg"], ["b.JPG"])
    assert not out["filename_integrity_gate_passed"]


def test_disjoint_unique_names_pass_only_the_filename_gate():
    out = summarize_splits(["a.JPG"], ["b.JPG"])
    assert out["filename_integrity_gate_passed"]


@pytest.mark.parametrize("train,valid", [([], ["b.JPG"]), (["a.JPG"], [])])
def test_missing_split_does_not_pass(train, valid):
    with pytest.raises(ValueError):
        summarize_splits(train, valid)


def make_bundle(tmp_path):
    files = []
    for role, raw in [("train", b"a.JPG\n"), ("validation", b"b.JPG\n")]:
        name = f"{role}.txt"
        (tmp_path / name).write_bytes(raw)
        files.append({"role": role, "local_filename": name, "bytes": len(raw),
                      "sha256": hashlib.sha256(raw).hexdigest(),
                      "provider_md5": hashlib.md5(raw, usedforsecurity=False).hexdigest()})
    manifest = {"schema_version": "rgfca-reference-split-source-v1",
                "source_record": "artificial_test_only", "files": files}
    (tmp_path / "source_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return manifest


def test_source_verified_does_not_authorize_benchmark_or_independence(tmp_path):
    make_bundle(tmp_path)
    out = audit_bundle(tmp_path)
    assert out["filename_integrity_gate_passed"]
    for key in ("benchmark_execution_authorized", "image_content_identity_verified",
                "event_observer_training_independence_verified", "image_pixels_decoded",
                "model_run", "coordinates_joined", "ecological_inference_performed"):
        assert out[key] is False


@pytest.mark.parametrize("mutation", ["bytes", "sha256", "provider_md5", "role", "path"])
def test_tampered_sources_or_manifest_fail_closed(tmp_path, mutation):
    manifest = make_bundle(tmp_path)
    first = manifest["files"][0]
    if mutation == "path":
        first["local_filename"] = "../train.txt"
    elif mutation == "role":
        first["role"] = "validation"
    elif mutation == "bytes":
        first["bytes"] += 1
    else:
        first[mutation] = "0" * len(first[mutation])
    (tmp_path / "source_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError):
        audit_bundle(tmp_path)


def test_source_newline_drift_is_detected_before_parsing(tmp_path):
    make_bundle(tmp_path)
    (tmp_path / "train.txt").write_bytes(b"a.JPG\r\n")
    with pytest.raises(ValueError, match="SHA-256"):
        audit_bundle(tmp_path)


def test_retained_usda_lists_reproduce_the_complete_non_disjoint_result():
    bundle = Path(__file__).resolve().parents[1] / "data/validation/usda_flower_split_audit_v1"
    out = audit_bundle(bundle)
    assert out == json.loads((bundle / "result.json").read_bytes())
    assert out["profiles"]["train"]["row_count"] == 100
    assert out["profiles"]["validation"]["row_count"] == 30
    assert out["exact_shared_filenames"] == ["IMG_0339.JPG"]
    assert out["exact_union_filename_count"] == 129
    assert not out["filename_split_disjoint"]
    assert not out["filename_integrity_gate_passed"]


def test_cli_rejects_changed_saved_result_without_writing(tmp_path, monkeypatch):
    from scripts.analysis.audit_rgfca_reference_splits import main
    make_bundle(tmp_path)
    result = audit_bundle(tmp_path)
    result["benchmark_execution_authorized"] = True
    saved = tmp_path / "result.json"
    saved.write_text(json.dumps(result), encoding="utf-8")
    before = saved.read_bytes()
    monkeypatch.setattr("sys.argv", ["audit", str(tmp_path), "--verify", str(saved)])
    with pytest.raises(ValueError, match="recomputed"):
        main()
    assert saved.read_bytes() == before
