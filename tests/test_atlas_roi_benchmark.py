from __future__ import annotations

import numpy as np

from fcp_pipeline.roi_benchmark import score_flower_weights, summarize_roi_benchmark


def test_roi_metric_excludes_unlabelled_pixels_and_rewards_flower_weight() -> None:
    rgb = np.array(
        [
            [[255, 0, 0], [255, 0, 0]],
            [[0, 255, 0], [0, 0, 255]],
        ],
        dtype=np.uint8,
    )
    trimap = np.array([[1, 1], [3, 0]], dtype=np.uint8)
    perfect = score_flower_weights(
        rgb,
        np.array([[1.0, 1.0], [0.0, 1.0]]),
        trimap,
    )
    assert perfect["soft_precision"] == 1.0
    assert perfect["soft_recall"] == 1.0
    assert perfect["soft_iou"] == 1.0
    assert perfect["colour_delta_e"] < 1e-6

    contaminated = score_flower_weights(
        rgb,
        np.array([[1.0, 1.0], [1.0, 0.0]]),
        trimap,
    )
    assert contaminated["soft_precision"] < perfect["soft_precision"]
    assert contaminated["soft_iou"] < perfect["soft_iou"]
    assert contaminated["colour_delta_e"] > perfect["colour_delta_e"]


def test_roi_metric_pools_official_variable_background_palette_indices() -> None:
    rgb = np.full((2, 3, 3), 127, dtype=np.uint8)
    trimap = np.array([[1, 0, 2], [3, 4, 9]], dtype=np.uint8)
    metrics = score_flower_weights(
        rgb,
        np.array([[1.0, 1.0, 0.5], [0.5, 0.5, 1.0]]),
        trimap,
        background_labels=(2, 3, 4),
    )
    assert metrics["foreground_pixels"] == 1.0
    assert metrics["background_pixels"] == 3.0
    # Black index 0 and an unknown index are excluded, not counted as background.
    assert metrics["scored_weight"] == 2.5


def test_roi_summary_applies_every_predeclared_gate() -> None:
    rows = [
        {
            "estimator_admitted": True,
            "soft_precision": 0.9,
            "soft_recall": 0.8,
            "soft_iou": 0.75,
            "colour_delta_e": 2.0,
        }
        for _ in range(10)
    ]
    gates = {
        "minimum_scored_images": 10,
        "minimum_admitted_fraction": 0.8,
        "minimum_median_soft_precision": 0.7,
        "minimum_median_soft_recall": 0.35,
        "minimum_median_soft_iou": 0.3,
        "maximum_median_colour_delta_e": 10.0,
        "maximum_p90_colour_delta_e": 20.0,
    }
    passed = summarize_roi_benchmark(rows, gates)
    assert passed["status"] == "pass_independent_roi_benchmark"

    rows[0]["estimator_admitted"] = False
    rows[1]["estimator_admitted"] = False
    rows[2]["estimator_admitted"] = False
    stopped = summarize_roi_benchmark(rows, gates)
    assert stopped["status"] == "stop_roi_benchmark_failed"
    assert stopped["checks"]["minimum_admitted_fraction"] is False
