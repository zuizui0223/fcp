from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.measurement_validity_v2 import (
    EV_LEVELS,
    apply_exposure_ev,
    exposure_metrics,
    frozen_prompt_jitter_set,
    mask_boundary_fraction,
    mask_iou,
    neutralize_background,
    seal_technical_table,
    validate_response_blind_columns,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "docs" / "supporting" / "fcp_v2_measurement_validity_contract_v1.json"


def test_contract_is_frozen_and_measurement_first() -> None:
    x = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert x["status"] == "FROZEN_BEFORE_FRESH_COHORT_PIXEL_OPENING"
    assert x["target_total_species"] == 400
    assert x["target_total_rows"] == 40000
    assert x["technical_workers_receive_species"] is False
    assert x["technical_workers_receive_coordinates"] is False
    assert x["biological_outcomes_open_before_technical_seal"] is False
    assert tuple(x["exposure"]["ev_levels"]) == EV_LEVELS
    assert x["background"]["neutral_rgb_srgb"] == [128, 128, 128]
    assert x["roi"]["random_jitter"] is False
    assert x["ecological_predictors_allowed"] is False


def test_exposure_transform_is_monotone_and_zero_ev_is_identity() -> None:
    rgb = np.array(
        [
            [[16, 32, 64], [120, 150, 180]],
            [[220, 230, 240], [250, 250, 250]],
        ],
        dtype=np.uint8,
    )
    zero = apply_exposure_ev(rgb, 0.0)
    assert np.array_equal(zero, rgb)

    darker = apply_exposure_ev(rgb, -1.0)
    lighter = apply_exposure_ev(rgb, +1.0)
    assert np.all(darker <= rgb)
    assert np.all(lighter >= rgb)
    assert lighter.max() == 255


def test_exposure_metrics_detect_clipping_and_use_linear_luminance() -> None:
    pixels = np.array(
        [
            [255, 10, 10],
            [250, 20, 20],
            [249, 249, 249],
            [0, 0, 0],
        ],
        dtype=np.uint8,
    )
    x = exposure_metrics(pixels)
    assert x["clip_fraction"] == pytest.approx(0.25)
    assert x["near_clip_fraction"] == pytest.approx(0.50)
    assert x["clip_fraction_r"] == pytest.approx(0.25)
    assert x["black_fraction"] == pytest.approx(0.25)
    assert 0 <= x["luminance_q01"] <= x["luminance_q99"] <= 1
    assert x["dynamic_range_q99_minus_q01"] >= 0


def test_background_neutralization_preserves_flower_pixels_exactly() -> None:
    rgb = np.arange(4 * 5 * 3, dtype=np.uint8).reshape(4, 5, 3)
    mask = np.zeros((4, 5), dtype=bool)
    mask[1:3, 2:4] = True
    out = neutralize_background(rgb, mask)
    assert np.array_equal(out[mask], rgb[mask])
    assert np.all(out[~mask] == np.array([128, 128, 128], dtype=np.uint8))


def test_roi_helpers_are_deterministic() -> None:
    a = np.zeros((5, 5), dtype=bool)
    b = np.zeros((5, 5), dtype=bool)
    a[1:4, 1:4] = True
    b[1:4, 2:5] = True
    assert mask_iou(a, a) == 1.0
    assert 0.0 < mask_iou(a, b) < 1.0
    assert 0.0 < mask_boundary_fraction(a) <= 1.0

    jitters = frozen_prompt_jitter_set((10, 20, 30, 60))
    assert len(jitters) == 7
    assert jitters[0] == (10.0, 20.0, 30.0, 60.0)
    assert jitters == frozen_prompt_jitter_set((10, 20, 30, 60))


def test_technical_firewall_rejects_biological_columns() -> None:
    validate_response_blind_columns(["measurement_id", "source_sha256", "near_clip_fraction"])
    for column in ("species", "morph", "latitude", "q_white_projection", "observer_id"):
        with pytest.raises(ValueError):
            validate_response_blind_columns(["measurement_id", column])


def test_technical_seal_is_row_order_invariant_and_outcome_blind() -> None:
    frame = pd.DataFrame(
        {
            "measurement_id": ["b", "a"],
            "source_sha256": ["2" * 64, "1" * 64],
            "near_clip_fraction": [0.2, 0.1],
        }
    )
    a = seal_technical_table(frame)
    b = seal_technical_table(frame.iloc[::-1].reset_index(drop=True))
    assert a.table_sha256 == b.table_sha256
    assert a.rows == 2
    assert a.unique_measurement_ids == 2
    assert a.biological_outcomes_opened is False


def test_technical_seal_rejects_duplicate_ids() -> None:
    frame = pd.DataFrame(
        {
            "measurement_id": ["x", "x"],
            "source_sha256": ["1" * 64, "2" * 64],
            "near_clip_fraction": [0.1, 0.2],
        }
    )
    with pytest.raises(ValueError):
        seal_technical_table(frame)
