from __future__ import annotations

import numpy as np

from fcp_pipeline.global_community_colour_assembly import matched_sympatry_colour_assembly_test


def _fixture(direction: str, n_focal: int = 40, sets_per_focal: int = 3, controls: int = 4):
    rng = np.random.default_rng(20260904)
    focal = np.repeat(np.arange(n_focal), sets_per_focal)
    n_sets = len(focal)
    if direction == "convergence":
        sympatric = rng.normal(0.15, 0.02, n_sets)
        allopatric = rng.normal(0.75, 0.04, (n_sets, controls))
    elif direction == "divergence":
        sympatric = rng.normal(0.80, 0.03, n_sets)
        allopatric = rng.normal(0.20, 0.03, (n_sets, controls))
    else:
        sympatric = rng.normal(0.50, 0.08, n_sets)
        allopatric = rng.normal(0.50, 0.08, (n_sets, controls))
    return focal, sympatric, allopatric


def test_strong_sympatric_convergence_is_detected():
    focal, sympatric, controls = _fixture("convergence")
    payload = matched_sympatry_colour_assembly_test(
        focal_species=focal,
        sympatric_colour_distance=sympatric,
        allopatric_control_colour_distance=controls,
        minimum_controls_per_set=2,
        minimum_sets_per_focal=2,
        minimum_focal_species=30,
        permutations=199,
        seed=17,
    )
    result = payload["global"]
    assert result.status == "evaluated"
    assert result.mean_focal_delta < -0.45
    assert result.convergent_focal_fraction > 0.95
    assert result.p_lower <= 0.01
    assert payload["directional_interpretation"]["convergence"] is True
    assert payload["directional_interpretation"]["divergence"] is False


def test_strong_sympatric_divergence_is_detected():
    focal, sympatric, controls = _fixture("divergence")
    payload = matched_sympatry_colour_assembly_test(
        focal_species=focal,
        sympatric_colour_distance=sympatric,
        allopatric_control_colour_distance=controls,
        minimum_controls_per_set=2,
        minimum_sets_per_focal=2,
        minimum_focal_species=30,
        permutations=199,
        seed=19,
    )
    result = payload["global"]
    assert result.status == "evaluated"
    assert result.mean_focal_delta > 0.45
    assert result.divergent_focal_fraction > 0.95
    assert result.p_upper <= 0.01
    assert payload["directional_interpretation"]["divergence"] is True
    assert payload["directional_interpretation"]["convergence"] is False


def test_species_level_characterization_uses_fdr_labels():
    focal, sympatric, controls = _fixture("convergence")
    payload = matched_sympatry_colour_assembly_test(
        focal_species=focal,
        sympatric_colour_distance=sympatric,
        allopatric_control_colour_distance=controls,
        minimum_controls_per_set=2,
        minimum_sets_per_focal=2,
        minimum_focal_species=30,
        permutations=199,
        seed=23,
        species_fdr_alpha=0.10,
    )
    labels = [entry["label"] for entry in payload["species"].values()]
    assert labels.count("convergent_in_sympatry") >= 35
    assert all(entry["delta_sympatric_minus_allopatric"] < 0 for entry in payload["species"].values())


def test_too_few_focal_species_is_not_evaluable():
    focal, sympatric, controls = _fixture("convergence", n_focal=10)
    payload = matched_sympatry_colour_assembly_test(
        focal_species=focal,
        sympatric_colour_distance=sympatric,
        allopatric_control_colour_distance=controls,
        minimum_controls_per_set=2,
        minimum_sets_per_focal=2,
        minimum_focal_species=30,
        permutations=99,
    )
    result = payload["global"]
    assert result.status == "not_evaluable_focal_species_coverage"
    assert result.n_focal_species == 10


def test_rows_without_enough_controls_are_dropped_before_focal_gate():
    focal, sympatric, controls = _fixture("convergence")
    controls = controls.copy()
    controls[:9, 1:] = np.nan
    payload = matched_sympatry_colour_assembly_test(
        focal_species=focal,
        sympatric_colour_distance=sympatric,
        allopatric_control_colour_distance=controls,
        minimum_controls_per_set=2,
        minimum_sets_per_focal=2,
        minimum_focal_species=30,
        permutations=99,
    )
    result = payload["global"]
    assert result.status == "evaluated"
    assert result.n_matched_sets == len(focal) - 9
