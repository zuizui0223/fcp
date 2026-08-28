import pandas as pd
import pytest

from fcp_pipeline.photo_calibration import (
    build_calibration_sheet,
    calibration_summary,
    validate_calibration_sheet,
)


def make_split() -> pd.DataFrame:
    rows = []
    for s in range(6):
        species = f"Species {s+1}"
        for i in range(200):
            rows.append(
                {
                    "species": species,
                    "photo_id": f"p{s+1:02d}_{i+1:03d}",
                    "photo_url": f"https://example.org/{s+1}/{i+1}.jpg",
                    "latitude": 10 + s,
                    "longitude": 20 + i / 100,
                    "observer": f"u{i%4}",
                    "observed_on": "2026-05-01",
                    "split": "calibration" if i < 80 else "evaluation",
                }
            )
    return pd.DataFrame(rows)


def test_build_sheet_contains_only_480_calibration_rows():
    sheet = build_calibration_sheet(make_split())
    assert len(sheet) == 480
    assert sheet.groupby("species").size().eq(80).all()
    assert sheet["photo_id"].is_unique


def test_default_sheet_hides_geography_observer_and_date():
    sheet = build_calibration_sheet(make_split())
    for forbidden in ("latitude", "longitude", "observer", "observed_on"):
        assert forbidden not in sheet.columns
    assert "photo_url" in sheet.columns


def test_explicit_display_columns_reject_location_leakage():
    with pytest.raises(ValueError, match="may not expose"):
        build_calibration_sheet(
            make_split(),
            keep_columns=["species", "photo_id", "photo_url", "latitude"],
        )


def test_incomplete_sheet_is_valid_draft_but_not_complete():
    sheet = build_calibration_sheet(make_split())
    validate_calibration_sheet(sheet, require_complete=False)
    with pytest.raises(ValueError, match="flower_visibility"):
        validate_calibration_sheet(sheet, require_complete=True)


def test_completed_logical_state_is_accepted():
    sheet = build_calibration_sheet(make_split())
    sheet["flower_visibility"] = "evaluable"
    sheet["segmentation_status"] = "ok"
    sheet["colour_assignment"] = "resolved"
    sheet["colour_state"] = "state_A"
    validate_calibration_sheet(sheet, require_complete=True)


def test_not_evaluable_requires_explicit_failure_and_no_colour():
    sheet = build_calibration_sheet(make_split())
    sheet["flower_visibility"] = "not_evaluable"
    sheet["visibility_failure_code"] = "blur"
    sheet["segmentation_status"] = "not_applicable"
    sheet["colour_assignment"] = "not_applicable"
    validate_calibration_sheet(sheet, require_complete=True)


def test_summary_reports_zero_measurement_at_generation():
    sheet = build_calibration_sheet(make_split())
    summary = calibration_summary(sheet)
    assert len(summary) == 6
    assert summary["n_total"].eq(80).all()
    assert summary["n_visibility_completed"].eq(0).all()
    assert summary["n_colour_resolved"].eq(0).all()
