from __future__ import annotations

import numpy as np

from fcp_pipeline.shared_transition_surface import EqualAreaGrid
from scripts.data import run_jbi_atlas_environmental_inference_v5 as env_v5


def _flower_row(mid: str, species: str, lat: float, lon: float, z: float) -> dict[str, object]:
    return {
        "measurement_id": mid,
        "species_blind_id": species,
        "cohort_id": "C01",
        "latitude": lat,
        "longitude": lon,
        "observed_month": 6,
        "local_solar_quarter": 2,
        "flower_L_standardized": z,
        "flower_a_standardized": 0.5 * z,
        "flower_b_standardized": -0.25 * z,
    }


def test_surface_uses_already_standardized_flower_values(monkeypatch) -> None:
    # The environmental surface builder must not call the standardizer used by
    # the background diagnostic: the flower field was frozen upstream already.
    def forbidden(*args, **kwargs):
        raise AssertionError("flower field was restandardized downstream")

    monkeypatch.setattr(env_v5, "robust_standardize_lab", forbidden)
    rows = [
        _flower_row(f"m{i}", "sp1", 35.0 + 0.01 * i, 135.0 + 0.01 * i, float(i - 4))
        for i in range(9)
    ]
    surface = env_v5.build_surface_from_standardized_lab(
        rows,
        env_v5.FLOWER_FIELDS,
        grid=EqualAreaGrid(32, 16),
        scale_km=500,
        configuration="all_dates",
        knn_k=3,
        minimum_edges_per_cell=1,
        minimum_retained_edges=5,
        minimum_detectable_cells=1,
    )
    assert surface.status == "evaluable"
    assert surface.retained_edges >= 5
    assert surface.detectable_cells >= 1


def test_background_is_frozen_once_from_eligible_admitted_rows() -> None:
    flower_rows = [
        _flower_row("m1", "sp1", 35.0, 135.0, -1.0),
        _flower_row("m2", "sp1", 35.1, 135.1, 0.0),
        _flower_row("m3", "sp1", 35.2, 135.2, 1.0),
        _flower_row("m4", "sp2", 36.0, 136.0, -1.0),
        _flower_row("m5", "sp2", 36.1, 136.1, 1.0),
    ]
    measurements = {
        "m1": {
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": True,
            "background_L_mean": 10.0,
            "background_a_mean": 0.0,
            "background_b_mean": 0.0,
        },
        "m2": {
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": True,
            "background_L_mean": 20.0,
            "background_a_mean": 1.0,
            "background_b_mean": 1.0,
        },
        "m3": {
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": True,
            "background_L_mean": 30.0,
            "background_a_mean": 2.0,
            "background_b_mean": 2.0,
        },
        "m4": {
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": True,
            "background_L_mean": 5.0,
            "background_a_mean": 5.0,
            "background_b_mean": 5.0,
        },
        "m5": {
            "automated_colour_state_status": "automated_colour_state_admitted",
            "background_features_available": False,
            "background_L_mean": 8.0,
            "background_a_mean": 8.0,
            "background_b_mean": 8.0,
        },
    }
    rows, status = env_v5.freeze_background_standardized_field(flower_rows, measurements)
    assert [row["measurement_id"] for row in rows] == ["m1", "m2", "m3"]
    values = np.asarray(
        [
            [row["background_L_standardized"], row["background_a_standardized"], row["background_b_standardized"]]
            for row in rows
        ],
        dtype=float,
    )
    assert np.allclose(np.median(values, axis=0), 0.0)
    status_by_species = {row["species_blind_id"]: row for row in status}
    assert status_by_species["sp1"]["status"] == "background_field_evaluable"
    assert status_by_species["sp2"]["status"] == "not_evaluable"
    assert status_by_species["sp2"]["background_rows"] == 1


def test_background_rejects_nonadmitted_flower_support() -> None:
    flower_rows = [_flower_row("m1", "sp1", 35.0, 135.0, 0.0)]
    measurements = {
        "m1": {
            "automated_colour_state_status": "automated_colour_state_not_evaluable",
            "background_features_available": True,
            "background_L_mean": 10.0,
            "background_a_mean": 0.0,
            "background_b_mean": 0.0,
        }
    }
    try:
        env_v5.freeze_background_standardized_field(flower_rows, measurements)
    except RuntimeError as exc:
        assert "not admitted" in str(exc)
    else:
        raise AssertionError("nonadmitted measurement was accepted into frozen background support")
