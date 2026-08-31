from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.atlas_environment import (
    categorical_boundary,
    composition_boundary,
    continuous_boundary,
    environmental_boundary_surfaces,
    rook_adjacency_without_repair,
    validate_environment_contract,
)


CONTRACT = Path("docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json")


def test_environment_contract_is_precolour_and_fixed() -> None:
    validate_environment_contract(json.loads(CONTRACT.read_text(encoding="utf-8")))


def test_rook_environment_boundaries_do_not_join_distant_islands() -> None:
    # Cells 0,1 are neighbours; cell 7 is isolated in a 4 x 2 grid.
    adjacency = rook_adjacency_without_repair([0, 1, 7], n_lon=4, n_sinlat=2)
    assert adjacency.tolist() == [
        [False, True, False],
        [True, False, False],
        [False, False, False],
    ]
    surface, scaling = continuous_boundary(
        np.array([[0.0], [2.0], [4.0]]), adjacency
    )
    assert surface[0] == surface[1]
    assert np.isnan(surface[2])
    assert scaling["sd"][0] > 0


def test_composition_and_categorical_boundaries_are_zero_for_matches() -> None:
    adjacency = np.array([[False, True], [True, False]])
    assert composition_boundary(np.array([[1.0, 0.0], [1.0, 0.0]]), adjacency).tolist() == [
        0.0,
        0.0,
    ]
    assert categorical_boundary(["A", "A"], adjacency).tolist() == [0.0, 0.0]
    assert categorical_boundary(["A", "B"], adjacency).tolist() == [1.0, 1.0]


def test_full_environment_builder_rejects_colour_and_builds_four_families() -> None:
    rows = []
    for index in range(4):
        row = {
            "bio1": index,
            "bio4": index * 2,
            "bio12": index * 3,
            "bio15": index * 4,
            "elevation": index,
            "slope": index * 2,
            "terrain_ruggedness": index * 3,
            "realm": f"R{index // 2}",
            "biome": f"B{index // 2}",
            "ecoregion": f"E{index}",
        }
        for code in (10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100):
            row[f"worldcover_{code}"] = 1.0 if code == (10 if index < 2 else 20) else 0.0
        rows.append(row)
    adjacency = np.array(
        [
            [False, True, False, False],
            [True, False, True, False],
            [False, True, False, True],
            [False, False, True, False],
        ]
    )
    result = environmental_boundary_surfaces(rows, adjacency)
    assert {"macroclimate", "terrain", "land_cover", "ecoregion"}.issubset(result)
    with pytest.raises(ValueError, match="flower-colour"):
        environmental_boundary_surfaces([dict(rows[0], flower_colour=1)], np.zeros((1, 1)))
