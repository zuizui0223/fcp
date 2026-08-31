from pathlib import Path

from PIL import Image

from scripts.data.build_inaturalist_development_review_packet import (
    PROHIBITED_REVIEW_FIELDS,
    REVIEW_FIELDS,
    assert_no_human_outcomes_before_rebuild,
    histogram_quantile,
    reusable_photos,
    reviewer_rows,
    technical_profile,
)


def test_reusable_photos_keeps_every_eligible_photo_in_stable_order():
    observation = {
        "photos": [
            {"id": 3, "license_code": "CC-BY", "url": "x/square.jpg"},
            {"id": 1, "license_code": "cc0", "url": "x/square.jpg"},
            {"id": 2, "license_code": "cc-by-nc", "url": "x/square.jpg"},
        ]
    }
    assert [photo["id"] for photo in reusable_photos(observation)] == [1, 3]


def test_reviewer_sheets_are_deterministic_independent_and_provenance_blind():
    encounters = [
        {
            "canonical_name": "Species alpha",
            "encounter_blind_id": f"FCP-{index:04d}",
            "image_files": [f"images/FCP-{index:04d}/P01.jpg"],
        }
        for index in range(30)
    ]
    a1 = reviewer_rows(encounters, "reviewer_A", "seed")
    a2 = reviewer_rows(reversed(encounters), "reviewer_A", "seed")
    b = reviewer_rows(encounters, "reviewer_B", "seed")
    assert a1 == a2
    assert [row["encounter_blind_id"] for row in a1] != [
        row["encounter_blind_id"] for row in b
    ]
    assert all(not PROHIBITED_REVIEW_FIELDS.intersection(row) for row in a1 + b)


def test_histogram_quantile_and_technical_profile(tmp_path: Path):
    histogram = [0] * 256
    histogram[10] = 1
    histogram[20] = 2
    histogram[30] = 1
    assert histogram_quantile(histogram, 0) == 10
    assert histogram_quantile(histogram, 0.5) == 20
    assert histogram_quantile(histogram, 1) == 30

    path = tmp_path / "image.png"
    Image.new("RGB", (20, 10), (10, 20, 30)).save(path)
    profile = technical_profile(path)
    assert profile["image_width"] == 20
    assert profile["image_height"] == 10
    assert profile["aspect_ratio"] == 2
    assert profile["luminance_p01"] == profile["luminance_p99"]
    assert 0 <= profile["fraction_luminance_le_5"] <= 1
    assert 0 <= profile["fraction_luminance_ge_250"] <= 1


def test_rebuild_refuses_to_overwrite_human_outcomes(tmp_path: Path):
    path = tmp_path / "reviewer_A_annotation_sheet.csv"
    row = {field: "" for field in REVIEW_FIELDS}
    row.update(
        {
            "review_order": "1",
            "canonical_name": "Species alpha",
            "encounter_blind_id": "FCP-0001",
            "image_files": "images/FCP-0001/P01.jpg",
            "n_photos": "1",
            "reviewer_id": "reviewer-1",
        }
    )
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REVIEW_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    try:
        assert_no_human_outcomes_before_rebuild(tmp_path)
    except RuntimeError as error:
        assert "refusing to overwrite" in str(error)
    else:
        raise AssertionError("rebuild accepted completed reviewer outcomes")
