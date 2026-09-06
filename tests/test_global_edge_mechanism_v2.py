from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import spearmanr

from fcp_pipeline.global_edge_mechanism_v2 import (
    edge_mechanism_permutation_test,
    evaluate_edge_mechanism_family,
    partial_spearman_distance,
)


def test_distance_confound_is_removed_by_partial_spearman():
    rng = np.random.default_rng(20260904)
    n = 500
    distance = np.linspace(1.0, 1000.0, n)
    colour = distance + rng.normal(0.0, 80.0, n)
    external = distance + rng.normal(0.0, 80.0, n)
    raw = float(spearmanr(colour, external).statistic)
    partial = partial_spearman_distance(colour, external, distance)
    assert raw > 0.85
    assert abs(partial) < 0.12


def test_beyond_distance_shared_signal_is_recovered():
    rng = np.random.default_rng(17)
    n = 500
    distance = np.linspace(1.0, 1000.0, n)
    shared = rng.normal(0.0, 180.0, n)
    colour = distance + shared + rng.normal(0.0, 25.0, n)
    external = distance + shared + rng.normal(0.0, 25.0, n)
    partial = partial_spearman_distance(colour, external, distance)
    assert partial > 0.85


def _family_fixture(n_species: int = 40, edges_per_species: int = 10):
    rng = np.random.default_rng(31)
    species = np.repeat(np.arange(n_species), edges_per_species)
    distance = np.tile(np.linspace(50.0, 950.0, edges_per_species), n_species)
    shared = rng.normal(0.0, 300.0, len(species))
    colour = distance + shared + rng.normal(0.0, 15.0, len(species))
    external = distance + shared + rng.normal(0.0, 15.0, len(species))
    null_rows = []
    for _ in range(199):
        row = colour.copy()
        for sid in range(n_species):
            idx = np.flatnonzero(species == sid)
            row[idx] = row[rng.permutation(idx)]
        null_rows.append(row)
    return colour, np.vstack(null_rows), external, distance, species


def test_permutation_test_detects_repeated_beyond_distance_alignment():
    colour, null, external, distance, species = _family_fixture()
    result = edge_mechanism_permutation_test(
        predictor_name="climate_turnover",
        colour_scores=colour,
        null_colour_scores=null,
        external_edge_scores=external,
        geographic_edge_distance_km=distance,
        species_index=species,
        minimum_edges_per_species=5,
        minimum_species=30,
    )
    assert result.status == "evaluated"
    assert result.mean_species_partial_rho > 0.75
    assert result.positive_species_fraction > 0.95
    assert result.p_upper <= 0.01


def test_too_few_species_is_not_evaluable():
    colour, null, external, distance, species = _family_fixture(n_species=10)
    result = edge_mechanism_permutation_test(
        predictor_name="pollinator_turnover",
        colour_scores=colour,
        null_colour_scores=null,
        external_edge_scores=external,
        geographic_edge_distance_km=distance,
        species_index=species,
        minimum_edges_per_species=5,
        minimum_species=30,
    )
    assert result.status == "not_evaluable_species_coverage"
    assert result.n_species == 10


def test_negative_geographic_distance_is_rejected():
    colour, null, external, distance, species = _family_fixture()
    distance = distance.copy()
    distance[0] = -1.0
    with pytest.raises(ValueError, match="distance"):
        edge_mechanism_permutation_test(
            predictor_name="terrain",
            colour_scores=colour,
            null_colour_scores=null,
            external_edge_scores=external,
            geographic_edge_distance_km=distance,
            species_index=species,
        )


def test_family_evaluation_uses_holm_and_is_independent_of_g1():
    colour, null, external, distance, species = _family_fixture()
    rng = np.random.default_rng(99)
    noise = rng.normal(size=len(colour))
    payload = evaluate_edge_mechanism_family(
        colour_scores=colour,
        null_colour_scores=null,
        external_edge_scores={
            "climate_turnover": external,
            "noise": noise,
        },
        geographic_edge_distance_km=distance,
        species_index=species,
        minimum_edges_per_species=5,
        minimum_species=30,
        alpha=0.05,
    )
    assert payload["status"] == "evaluated"
    assert payload["does_not_require_G1"] is True
    assert payload["cannot_rescue_null_G1"] is True
    assert "partial Spearman" in payload["statistic"]
    assert payload["results"]["climate_turnover"]["supported"] is True
    assert payload["results"]["climate_turnover"]["p_holm"] <= 0.05
