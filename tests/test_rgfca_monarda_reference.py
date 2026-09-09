"""Offline receipt reconstruction and artificial-input checks; no image decoder."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import zipfile

import pytest

from scripts.analysis import audit_rgfca_monarda_archive as archive_audit
from scripts.analysis import map_rgfca_monarda_reference_sources as source_audit


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "docs/supporting"
INTAKE_PATH = SUPPORT / "rgfca_monarda_archive_intake_v1.json"
MAPPING_PATH = SUPPORT / "rgfca_monarda_source_mapping_v1.json"
INTAKE_SHA = "1f5acb5af939b83c76b75d6513a0f404019ea88d4a4924d6d1322c8b8cde70f8"
MAPPING_SHA = "fd0c58eb2657b8a7f5c620dca851e195ed341aa6f6781e7c704645b14a308106"


def read(path):
    return json.loads(path.read_bytes())


def test_exact_canonical_receipts_and_archive_identity():
    assert source_audit.sha(INTAKE_PATH) == INTAKE_SHA == source_audit.INTAKE_SHA
    assert source_audit.sha(MAPPING_PATH) == MAPPING_SHA
    intake, mapping = read(INTAKE_PATH), read(MAPPING_PATH)
    assert intake["archive_sha256"] == mapping["archive_sha256"] == archive_audit.EXPECTED_SHA
    assert mapping["intake_sha256"] == INTAKE_SHA
    assert intake["archive_bytes"] == 41055963
    assert len(intake["members"]) == 115
    assert sum(m["bytes"] for m in intake["members"]) == 41434856
    assert intake["archive_crc_all_members_verified"] is True


def test_complete_split_local_pair_census_and_annotation_denominator():
    intake = read(INTAKE_PATH)
    images = intake["image_inventory"]
    members = {m["path"]: m for m in intake["members"]}
    assert len(images) == intake["total_images"] == 110
    assert sum(i["annotation_count"] for i in images) == intake["total_annotations"] == 788
    assert len({(i["split"], i["coco_image_id"]) for i in images}) == 110
    assert len({i["sha256"] for i in images}) == len({i["original_name"] for i in images}) == 110
    assert intake["duplicate_image_byte_groups"] == intake["duplicate_original_name_groups"] == []
    for split, expected in {"train": (77, 537), "valid": (22, 134), "test": (11, 117)}.items():
        rows = [i for i in images if i["split"] == split]
        assert (len(rows), sum(i["annotation_count"] for i in rows)) == expected
        assert 0 in {i["coco_image_id"] for i in rows}  # IDs repeat between splits.
        assert intake["splits"][split]["used_category_ids"] == {"1": expected[1]}
        assert [c["name"] for c in intake["splits"][split]["categories"]] == ["flowers", "flowers"]
    for row in images:
        assert row["bytes"] == members[row["path"]]["bytes"]
        assert row["sha256"] == members[row["path"]]["sha256"]


def test_empty_annotation_and_export_dates_are_retained_not_relabelled():
    intake = read(INTAKE_PATH)
    empty = [r for r in intake["image_inventory"] if r["annotation_count"] == 0]
    assert [(r["split"], r["coco_image_id"], r["original_name"]) for r in empty] == [("train", 69, "22132.jpg")]
    assert all(r["date_captured_equals_export_date_created"] for r in intake["image_inventory"])
    assert intake["checks"] == {"polygon_components": 788, "date_captured_equals_export_date_created": 110, "images_without_annotations": 1}
    assert "not_admitted_for_model_evaluation" in intake["decision"]
    assert intake["read_depth"]["zip_bytes_and_coco_json_read"] is True
    assert all(v is False for k, v in intake["read_depth"].items() if k != "zip_bytes_and_coco_json_read")


def test_source_mapping_and_overlap_totals_reconstruct_for_all_images():
    intake, mapping = read(INTAKE_PATH), read(MAPPING_PATH)
    rows = mapping["image_source_mapping"]
    originals = {i["path"]: i["original_name"] for i in intake["image_inventory"]}
    assert len(rows) == mapping["reference_images_mapped"] == 110
    assert {r["path"] for r in rows} == set(originals)
    assert len({r["photo_id"] for r in rows}) == mapping["unique_photo_ids"] == 110
    assert len({r["observation_id"] for r in rows}) == mapping["unique_observation_ids"] == 110
    assert mapping["duplicate_photo_id_groups"] == mapping["duplicate_observation_id_groups"] == {}
    for row in rows:
        assert originals[row["path"]] == str(row["source_row_index_zero_based"]) + ".jpg"
        assert f'/photos/{row["photo_id"]}/' in row["source_image_url"]
        assert row["observation_url"].endswith("/" + row["observation_id"])
        assert not ({"latitude", "longitude", "eventDate", "observer", "colour"} & set(row))
    for key, frame in mapping["fcp_comparison_frames"].items():
        assert frame["rows"] == 50000
        assert frame["immutable_input_sha256"] == source_audit.FRAME_HASHES[frame["git_ref"]]
        for id_type in ("photo_id", "observation_id"):
            assert sum(r["overlap"][key][id_type] for r in rows) == frame[id_type + "_matches"] == 0


def test_historical_photo_rights_not_replaced_by_provider_blanket_notice():
    mapping = read(MAPPING_PATH)
    counts = dict(Counter(r["original_photo_license_from_2025_archive"] for r in mapping["image_source_mapping"]))
    assert counts == mapping["source_license_counts"]
    assert counts == {
        "http://creativecommons.org/licenses/by-nc/4.0/": 97,
        "http://creativecommons.org/licenses/by-nc-sa/4.0/": 1,
        "http://creativecommons.org/licenses/by/4.0/": 9,
        "http://creativecommons.org/publicdomain/zero/1.0/": 3,
    }


@pytest.fixture
def artificial_archive(tmp_path, monkeypatch):
    """Opaque strings, NOT actual JPEGs: successful intake cannot decode pixels."""
    def make(change=None, extra=None, omit=None):
        payloads = {"README.dataset.txt": b"untrusted fixture", "README.roboflow.txt": b"untrusted fixture"}
        for split, count in (("train", 77), ("valid", 22), ("test", 11)):
            images, annotations = [], []
            for index in range(count):
                filename = f"{index}.jpg"
                images.append({"id": index, "file_name": filename, "width": 10, "height": 10,
                               "license": 0, "date_captured": "export", "extra": {"name": filename}})
                annotations.append({"id": index, "image_id": index, "category_id": 1,
                                    "segmentation": [[0, 0, 2, 0, 0, 2]], "bbox": [0, 0, 2, 2]})
                payloads[f"{split}/{filename}"] = f"not JPEG {split} {index}".encode()
            coco = {"images": images, "annotations": annotations, "categories": [{"id": 1, "name": "flowers"}],
                    "info": {"date_created": "export"}, "licenses": []}
            if change and split == "train":
                change(coco)
            payloads[f"{split}/_annotations.coco.json"] = json.dumps(coco).encode()
        if extra:
            payloads.update(extra)
        if omit:
            payloads.pop(omit)
        path = tmp_path / "fixture.zip"
        with zipfile.ZipFile(path, "w") as archive:
            for name, data in payloads.items():
                entry = zipfile.ZipInfo(name)
                # Windows ZipInfo otherwise normalizes the unsafe backslash
                # before the reader-under-test can inspect it.
                entry.filename = name
                archive.writestr(entry, data)
        monkeypatch.setattr(archive_audit, "EXPECTED_SHA", source_audit.sha(path))
        return path
    return make


def test_intake_can_only_inspect_opaque_bytes_without_pixel_decode(artificial_archive):
    result = archive_audit.inspect_archive(artificial_archive())
    assert result["total_images"] == result["total_annotations"] == 110
    assert result["read_depth"]["pixel_decode"] is False


def test_changed_archive_is_rejected_before_parsing(artificial_archive):
    path = artificial_archive()
    with path.open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="pinned"):
        archive_audit.inspect_archive(path)


@pytest.mark.parametrize("name", ["../escape.txt", "C:/escape.txt", "train\\escape.txt", "/escape.txt"])
def test_unsafe_member_paths_rejected(artificial_archive, name):
    with pytest.raises(ValueError, match="Unsafe archive path"):
        archive_audit.inspect_archive(artificial_archive(extra={name: b"fixture"}))


def test_case_collision_rejected(artificial_archive):
    with pytest.raises(ValueError, match="case-colliding"):
        archive_audit.inspect_archive(artificial_archive(extra={"TRAIN/0.JPG": b"fixture"}))


def test_missing_coco_image_pair_rejected(artificial_archive):
    with pytest.raises(ValueError, match="Missing"):
        archive_audit.inspect_archive(artificial_archive(omit="train/0.jpg"))


def test_unreferenced_image_member_rejected(artificial_archive):
    with pytest.raises(ValueError, match="census differ"):
        archive_audit.inspect_archive(artificial_archive(extra={"orphan.jpg": b"fixture"}))


@pytest.mark.parametrize("change,message", [
    (lambda c: c["images"][0].update(id=1), "Duplicate COCO"),
    (lambda c: c["annotations"][0].update(image_id=999), "Orphan"),
    (lambda c: c["images"][0].update(width=0), "dimensions"),
    (lambda c: c["images"].pop(), "Split count"),
])
def test_coco_integrity_errors_fail_closed(artificial_archive, change, message):
    with pytest.raises(ValueError, match=message):
        archive_audit.inspect_archive(artificial_archive(change=change))


@pytest.mark.parametrize("segmentation,check", [
    ([[0, 0, 1]], "invalid_polygon_arrays"),
    ([[0, 0, 20, 0, 0, 2]], "polygons_outside_declared_dimensions"),
    ([], "empty_polygon_segmentations"),
    ({"counts": "fixture"}, "non_polygon_segmentations_not_assessed"),
])
def test_unsupported_or_invalid_geometry_is_recorded_not_admitted(artificial_archive, segmentation, check):
    result = archive_audit.inspect_archive(artificial_archive(change=lambda c: c["annotations"][0].update(segmentation=segmentation)))
    assert result["checks"][check] == 1
    assert "not_admitted_for_model_evaluation" in result["decision"]


@pytest.mark.parametrize("value", ["123.5", "NaN", "1e4", "-1", ""])
def test_identifier_rejects_lossy_or_missing_ids(value):
    with pytest.raises(ValueError, match="Non-integral"):
        source_audit.identifier(value)


def test_identifier_retains_large_integer_precision():
    assert source_audit.identifier("9007199254740993.0") == "9007199254740993"


def test_frame_requires_pinned_hash_and_complete_denominator(monkeypatch):
    raw = b"photo_id,observation_id\n1,2\n"
    monkeypatch.setattr(source_audit.subprocess, "check_output", lambda args: raw)
    monkeypatch.setattr(source_audit, "FRAME_HASHES", {"fixture": "wrong"})
    with pytest.raises(ValueError, match="hash mismatch"):
        source_audit.frame_metadata("fixture")
    monkeypatch.setattr(source_audit, "FRAME_HASHES", {"fixture": hashlib.sha256(raw).hexdigest()})
    with pytest.raises(ValueError, match="complete 50,000"):
        source_audit.frame_metadata("fixture")


def test_mapping_joins_observations_via_gbif_not_photo_reference(tmp_path, monkeypatch):
    # Source references deliberately point to PHOTOS, not observations.
    author = tmp_path / "multimedia.tsv"
    occurrence = tmp_path / "occurrence.tsv"
    author.write_text("gbifID\tidentifier\treferences\tlicense\n" + "".join(
        f"{i}\thttps://example.test/photos/{100000+i}/original.jpg\thttps://inaturalist.org/photos/{100000+i}\tfixture\n"
        for i in range(41069)), encoding="utf-8")
    occurrence.write_text("gbifID\toccurrenceID\n" + "".join(
        f"{i}\thttps://www.inaturalist.org/observations/{300000+i}\n" for i in range(41069)), encoding="utf-8")
    monkeypatch.setattr(source_audit, "AUTHOR_SHA", source_audit.sha(author))
    monkeypatch.setattr(source_audit, "OCCURRENCE_SHA", source_audit.sha(occurrence))
    monkeypatch.setattr(source_audit, "frame_metadata", lambda ref: {
        "sha256": "fixture", "git_ref": ref, "rows": 50000, "photo_ids": set(), "observation_ids": set()})
    result = source_audit.audit(INTAKE_PATH, author, occurrence)
    assert result["reference_images_mapped"] == 110
    for row in result["image_source_mapping"]:
        index = row["source_row_index_zero_based"]
        assert row["photo_id"] == str(100000 + index)
        assert row["observation_id"] == str(300000 + index)
    broken = tmp_path / "incomplete.json"
    broken.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="complete intake receipt"):
        source_audit.audit(broken, author, occurrence)
