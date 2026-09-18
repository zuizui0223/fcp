import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from disttrait import (
    matched_difference_spatial_permutation_null,
    observer_disjoint_reliability,
    one_vs_rest_contrast,
    spatial_permutation_null,
    structured_alignment_null,
    two_mode_axis,
)


ROOT = Path(__file__).resolve().parents[3]


def _h1_fixture() -> pd.DataFrame:
    rows = []
    patterns = {
        "spA": {
            "o1": ["white"] * 8 + ["yellow_orange"] * 2,
            "o2": ["white"] * 6 + ["yellow_orange"] * 4,
            "o3": ["white"] * 4 + ["yellow_orange"] * 6,
            "o4": ["white"] * 2 + ["yellow_orange"] * 8,
        },
        "spB": {
            "o1": ["white"] * 5 + ["red_pink"] * 5,
            "o2": ["white"] * 8 + ["red_pink"] * 2,
            "o3": ["white"] * 2 + ["red_pink"] * 8,
            "o4": ["white"] * 5 + ["red_pink"] * 5,
        },
        "spC": {
            "o1": ["blue_purple"] * 10,
            "o2": ["blue_purple"] * 8 + ["white"] * 2,
            "o3": ["blue_purple"] * 6 + ["white"] * 4,
            "o4": ["blue_purple"] * 4 + ["white"] * 6,
        },
        "spD": {
            "o1": ["white"] * 3 + ["yellow_orange"] * 3 + ["red_pink"] * 4,
            "o2": ["white"] * 4 + ["yellow_orange"] * 3 + ["red_pink"] * 3,
            "o3": ["white"] * 3 + ["yellow_orange"] * 4 + ["red_pink"] * 3,
            "o4": ["white"] * 4 + ["yellow_orange"] * 4 + ["red_pink"] * 2,
        },
        "spE": {
            "o1": ["white"] * 7 + ["yellow_orange", "red_pink", "blue_purple"],
            "o2": ["yellow_orange"] * 7 + ["white", "red_pink", "blue_purple"],
            "o3": ["red_pink"] * 7 + ["white", "yellow_orange", "blue_purple"],
            "o4": ["blue_purple"] * 7 + ["white", "yellow_orange", "red_pink"],
        },
    }
    for species, observers in patterns.items():
        for observer, values in observers.items():
            for state in values:
                rows.append(
                    {
                        "species": species,
                        "observer_id": observer,
                        "morph": state,
                        "global_classifiable": True,
                    }
                )
            # An allowed-looking state that must be excluded by the explicit
            # classifiability flag.
            rows.append(
                {
                    "species": species,
                    "observer_id": observer,
                    "morph": "white",
                    "global_classifiable": False,
                }
            )
        # Blank observers are excluded by the frozen FCP H1 rule.
        rows.append(
            {
                "species": species,
                "observer_id": "  ",
                "morph": "white",
                "global_classifiable": True,
            }
        )
    return pd.DataFrame(rows)


def test_h1_observer_split_matches_frozen_fcp_algorithm_fixture() -> None:
    result = observer_disjoint_reliability(
        _h1_fixture(),
        species_col="species",
        observer_col="observer_id",
        state_col="morph",
        classifiable_col="global_classifiable",
        states=["white", "yellow_orange", "red_pink", "blue_purple"],
        n_partitions=7,
        base_seed=500,
        min_full_classifiable=8,
        min_half_classifiable=4,
    )

    expected_rho = np.array(
        [
            0.8207826816681233,
            0.8999999999999998,
            0.9999999999999999,
            0.9746794344808964,
            0.7999999999999999,
            0.8207826816681233,
            0.8947368421052632,
        ]
    )
    expected_ccc = np.array(
        [
            0.5241451419277388,
            0.8125057218712805,
            1.0,
            0.8017680841453182,
            0.8248092681548456,
            0.5613465176882159,
            0.9992429977289935,
        ]
    )
    np.testing.assert_allclose(result.partitions["rho"], expected_rho, rtol=0, atol=1e-14)
    np.testing.assert_allclose(result.partitions["ccc"], expected_ccc, rtol=0, atol=1e-14)
    assert result.summary["eligible_species"] == 5.0
    assert result.summary["paired_n_median"] == 5.0
    assert result.summary["rho_median"] == pytest.approx(0.8947368421052632, abs=1e-14)
    assert result.summary["rho_q05"] == pytest.approx(0.8062348045004369, abs=1e-14)
    assert result.summary["rho_q95"] == pytest.approx(0.9924038303442688, abs=1e-14)
    assert result.summary["ccc_median"] == pytest.approx(0.8125057218712805, abs=1e-14)
    assert result.summary["spearman_brown_median"] == pytest.approx(0.9444444444444444, abs=1e-14)
    assert result.summary["mae_median"] == pytest.approx(0.03200000000000003, abs=1e-14)


def test_h2_two_means_matches_frozen_fcp_algorithm_fixture() -> None:
    compositions = np.array(
        [
            [0.95, 0.03, 0.02],
            [0.90, 0.05, 0.05],
            [0.92, 0.04, 0.04],
            [0.05, 0.50, 0.45],
            [0.02, 0.48, 0.50],
            [0.04, 0.46, 0.50],
        ],
        dtype=float,
    )
    result = two_mode_axis(compositions)
    np.testing.assert_array_equal(result.labels, [0, 0, 0, 1, 1, 1])
    assert result.cluster_sizes == (3, 3)
    assert result.minor_fraction == 0.5
    assert result.separation_ratio == pytest.approx(25.56850550988829, abs=1e-12)
    np.testing.assert_allclose(
        result.unit_axis,
        [0.81648889, -0.40517494, -0.41131395],
        rtol=0,
        atol=5e-9,
    )


def test_structured_alignment_null_matches_frozen_fcp_construction_fixture() -> None:
    rows = []
    species = []
    strata = []
    fixture = {
        "sp1": [
            ("A", [0.90, 0.07, 0.03]),
            ("A", [0.85, 0.10, 0.05]),
            ("A", [0.88, 0.08, 0.04]),
            ("B", [0.10, 0.45, 0.45]),
            ("B", [0.08, 0.50, 0.42]),
            ("B", [0.12, 0.40, 0.48]),
        ],
        "sp2": [
            ("A", [0.75, 0.15, 0.10]),
            ("A", [0.78, 0.12, 0.10]),
            ("A", [0.72, 0.18, 0.10]),
            ("B", [0.20, 0.60, 0.20]),
            ("B", [0.18, 0.62, 0.20]),
            ("B", [0.22, 0.58, 0.20]),
        ],
        "sp3": [
            ("A", [0.65, 0.25, 0.10]),
            ("A", [0.68, 0.22, 0.10]),
            ("A", [0.62, 0.28, 0.10]),
            ("B", [0.30, 0.20, 0.50]),
            ("B", [0.28, 0.22, 0.50]),
            ("B", [0.32, 0.18, 0.50]),
        ],
    }
    for sp, values in fixture.items():
        for stratum, composition in values:
            species.append(sp)
            strata.append(stratum)
            rows.append(composition)

    result = structured_alignment_null(
        np.asarray(rows, dtype=float),
        species,
        strata,
        one_vs_rest_contrast(3, focal_index=0),
        n_permutations=7,
        seed=1234,
        strata_order=["A", "B"],
    )
    assert result.n_species == 3
    assert result.observed == pytest.approx(0.8415894074641564, abs=1e-14)
    np.testing.assert_allclose(
        result.null,
        [
            0.9795760222919392,
            0.9546391137953703,
            0.9761887862595559,
            0.9954539526745961,
            0.9970410115104974,
            0.9887474330894296,
            0.989360775244641,
        ],
        rtol=0,
        atol=1e-14,
    )
    assert result.p_upper == 1.0


def test_generic_one_vs_rest_contrast_matches_frozen_q_white() -> None:
    frozen = json.loads(
        (
            ROOT
            / "results"
            / "polymorphism_white_axis_targeted_test_20260912"
            / "result.json"
        ).read_text(encoding="utf-8")
    )
    order = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
    expected = np.array([frozen["fixed_axis_palette_loadings"][x] for x in order], dtype=float)
    observed = one_vs_rest_contrast(9, focal_index=0)
    np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-15)


def test_spatial_vertex_null_matches_frozen_rgfca_algorithm_fixture() -> None:
    latitude = [0.0] * 5
    longitude = [0.0, 1.0, 2.0, 3.0, 4.0]
    traits = np.array(
        [
            [1.0, 0.0],
            [0.85, 0.15],
            [0.50, 0.50],
            [0.15, 0.85],
            [0.0, 1.0],
        ]
    )
    observed, null = spatial_permutation_null(
        latitude,
        longitude,
        traits,
        n_permutations=7,
        seed=202609070901,
        key="fixture_sp",
    )
    assert observed == pytest.approx(0.9218930058986121, abs=1e-14)
    np.testing.assert_allclose(
        null,
        [
            -0.5031348269480561,
            0.4687591555416672,
            -0.10312701421916678,
            -0.18750366221666687,
            0.15937811288416684,
            0.4687591555416672,
            0.3125061036944448,
        ],
        rtol=0,
        atol=1e-14,
    )


def test_matched_background_joint_null_matches_frozen_reserve_algorithm_fixture() -> None:
    latitude = [0.0] * 5
    longitude = [0.0, 1.0, 2.0, 3.0, 4.0]
    focal = np.array(
        [[1.0, 0.0], [0.85, 0.15], [0.50, 0.50], [0.15, 0.85], [0.0, 1.0]]
    )
    background = np.array(
        [[0.6, 0.4], [0.5, 0.5], [0.4, 0.6], [0.5, 0.5], [0.6, 0.4]]
    )
    observed, null = matched_difference_spatial_permutation_null(
        latitude,
        longitude,
        focal,
        background,
        n_permutations=7,
        seed=202609071503,
        key="fixture_sp",
    )
    assert observed == pytest.approx(0.9218930058986121, abs=1e-14)
    np.testing.assert_allclose(
        null,
        [
            0.3718822633963893,
            -0.1375026856255557,
            0.31250610369444476,
            0.05000097659111116,
            -0.18750366221666687,
            0.28438055436194476,
            0.3375065919900004,
        ],
        rtol=0,
        atol=1e-14,
    )
