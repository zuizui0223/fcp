import json
import sys

import numpy as np
import pytest

from scripts.data.run_inaturalist_automated_colour_spatial_marks import (
    benjamini_hochberg,
    flower_background_contrast_test,
    main,
    random_mark_test,
    robust_standardize,
)


def test_robust_standardize_uses_componentwise_median_and_iqr():
    values = np.array([[0.0, 10.0, 20.0], [1.0, 12.0, 24.0], [2.0, 14.0, 28.0]])
    standardized = robust_standardize(values)
    assert np.allclose(np.median(standardized, axis=0), 0.0)
    assert np.allclose(
        np.quantile(standardized, 0.75, axis=0)
        - np.quantile(standardized, 0.25, axis=0),
        1.0,
    )


def test_random_mark_test_detects_monotone_spatial_turnover():
    latitude = np.zeros(30)
    longitude = np.linspace(0.0, 10.0, 30)
    marks = np.column_stack((longitude, longitude**2, np.sin(longitude / 3.0)))
    result = random_mark_test(
        latitude, longitude, marks, species="Synthetic species", permutations=999
    )
    assert result["rho"] > 0.8
    assert result["p_greater"] < 0.05


def test_flower_background_contrast_requires_flower_specific_structure():
    rng = np.random.default_rng(7)
    latitude = np.zeros(30)
    longitude = np.linspace(0.0, 10.0, 30)
    flower = np.column_stack((longitude, longitude**2, np.sin(longitude / 3.0)))
    background = rng.normal(size=(30, 3))
    result = flower_background_contrast_test(
        latitude,
        longitude,
        flower,
        background,
        species="Synthetic species",
        permutations=999,
    )
    assert result["flower_minus_background_rho"] > 0.5
    assert result["p_greater"] < 0.05


def test_benjamini_hochberg_is_monotone_in_rank():
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.20])
    assert np.allclose(adjusted, [0.04, 0.05333333333333334, 0.05333333333333334, 0.20])


def test_main_does_not_open_locked_table_when_no_species_passes(tmp_path, monkeypatch):
    development = tmp_path / "development.json"
    development.write_text(
        json.dumps(
            {
                "protocol": "fcp-inaturalist-automated-colour-state-v2",
                "status": "complete_location_free_automated_colour_development_gate",
                "species_results": [
                    {"canonical_name": "Species A", "development_gate_status": "not_evaluable"}
                ],
            }
        ),
        encoding="utf-8",
    )
    locked = tmp_path / "must-not-be-opened.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "runner",
            "--development-gate",
            str(development),
            "--locked-table",
            str(locked),
            "--output-dir",
            str(tmp_path / "out"),
            "--public-manifest",
            str(tmp_path / "public.json"),
        ],
    )
    with pytest.raises(RuntimeError, match="coordinates must remain unopened"):
        main()
    assert not locked.exists()
