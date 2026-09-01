from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fcp_pipeline.atlas_real_inference_v5 import (
    COHORTS,
    FrozenSpeciesColourState,
    run_shared_transition_test,
    run_spatial_organization_test,
    validate_real_inference_amendment,
)
from fcp_pipeline.shared_transition_surface import EdgeCellGeometry


REAL = Path("docs/supporting/jbi_atlas_real_colour_inference_amendment_v5.json")
INFERENCE = Path("docs/supporting/jbi_image_first_atlas_inference_contract_v5.json")


def _contracts():
    return (
        json.loads(REAL.read_text(encoding="utf-8")),
        json.loads(INFERENCE.read_text(encoding="utf-8")),
    )


def _geometry(*, n_cells: int = 20, edges_per_cell: int = 2) -> EdgeCellGeometry:
    n_edges = n_cells * edges_per_cell
    edges = np.arange(n_edges * 2, dtype=int).reshape(n_edges, 2)
    cell_id = np.repeat(np.arange(n_cells, dtype=int), edges_per_cell)
    counts = np.bincount(cell_id, minlength=n_cells)
    return EdgeCellGeometry(
        retained_edge_indices=np.arange(n_edges, dtype=int),
        retained_edges=edges,
        retained_edge_distance_km=np.full(n_edges, 20.0),
        edge_cell_id=cell_id,
        cell_edge_count=counts,
        detectable=counts >= edges_per_cell,
    )


def _state(species: str, cohort: str, *, signal_cell: int = 0, seed: int = 0) -> FrozenSpeciesColourState:
    geometry = _geometry()
    rng = np.random.default_rng(seed)
    lab = np.zeros((80, 3), dtype=float)
    # Give every cell a small, species-specific continuous edge score so the
    # top-tail is not driven by deterministic score ties. Cell 0 receives a
    # much larger discontinuity shared across species.
    for edge_index, (a, b) in enumerate(geometry.retained_edges):
        cell = int(geometry.edge_cell_id[edge_index])
        base = 0.02 * (cell + 1) + rng.uniform(-0.003, 0.003)
        if cell == signal_cell:
            base += 8.0
        direction = np.asarray([1.0, 0.37, -0.21])
        lab[a] = -0.5 * base * direction
        lab[b] = 0.5 * base * direction
    return FrozenSpeciesColourState(
        species_id=species,
        cohort_id=cohort,
        standardized_lab=lab,
        geometry=geometry,
    )


def test_real_colour_inference_amendment_is_frozen_pre_pixel() -> None:
    real, _ = _contracts()
    validate_real_inference_amendment(real)
    assert real["species_conditioned_spatial_organization"]["randomizations"] == 9999
    assert real["shared_transition"]["randomizations"] == 9999
    assert real["environmental_concordance"]["randomizations"] == 9999
    assert real["cascade"]["not_evaluable_never_advances_confirmatory"] is True


def test_spatial_test_is_reproducible_and_keeps_all_eight_cohorts() -> None:
    _, inference = _contracts()
    states = [
        _state(f"sp-{i}", cohort, seed=100 + i)
        for i, cohort in enumerate(COHORTS)
    ]
    first = run_spatial_organization_test(
        states,
        inference_v5=inference,
        randomizations=29,
        seed=1234,
        require_terminal=False,
    )
    second = run_spatial_organization_test(
        states,
        inference_v5=inference,
        randomizations=29,
        seed=1234,
        require_terminal=False,
    )
    assert first == second
    assert set(first["cohort_directions"]) == set(COHORTS)
    assert 0.0 < first["pooled_p_value"] <= 1.0
    assert first["randomizations"] == 29


def test_shared_transition_stops_without_preimage_qualification() -> None:
    _, inference = _contracts()
    states = [
        _state(f"sp-{i}", cohort, seed=200 + i)
        for i, cohort in enumerate(COHORTS)
    ]
    result = run_shared_transition_test(
        states,
        inference_v5=inference,
        qualification_passed=False,
        randomizations=9,
        require_terminal=False,
    )
    assert result["outcome"] == "not_evaluable"
    assert result["reason"] == "preimage_qualification_not_passed"


def test_shared_transition_real_null_is_lab_vector_permutation() -> None:
    _, inference = _contracts()
    states = []
    index = 0
    for cohort in COHORTS:
        for _ in range(4):
            states.append(_state(f"sp-{index}", cohort, seed=500 + index))
            index += 1
    result = run_shared_transition_test(
        states,
        inference_v5=inference,
        qualification_passed=True,
        randomizations=19,
        seed=4321,
        min_detectable_species=4,
        require_terminal=False,
    )
    assert result["qualification_passed"] is True
    assert result["real_colour_null"] == "complete_within_species_lab_vector_permutation"
    assert result["valid_pooled_cells"] == 20
    assert set(result["cohort_directions"]) == set(COHORTS)
    assert 0.0 < result["pooled_p_value"] <= 1.0
    assert result["randomizations"] == 19
