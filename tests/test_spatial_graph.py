import numpy as np
import pytest

from fcp_pipeline.spatial_graph import (
    SpeciesGraph,
    equal_species_graph_discontinuity,
    species_conditioned_graph_permutation_null,
    species_conditioned_knn_graph,
    spherical_knn_edges,
)


def test_spherical_knn_graph_is_undirected_and_respects_label_blind_truncation():
    lat = np.array([0.0, 0.0, 0.0, 10.0])
    lon = np.array([0.0, 0.1, 0.2, 10.0])

    edges, distance = spherical_knn_edges(lat, lon, k=1, max_edge_km=50.0)

    assert edges.ndim == 2 and edges.shape[1] == 2
    assert np.all(edges[:, 0] < edges[:, 1])
    assert len({tuple(edge) for edge in edges.tolist()}) == len(edges)
    assert np.all(distance <= 50.0)
    # The far point has no retained edge after the geometry-only distance filter.
    assert not any(3 in edge for edge in edges.tolist())


def test_species_conditioned_graph_never_connects_different_species():
    lat = np.array([0.0, 0.0, 0.1, 20.0, 20.0, 20.1])
    lon = np.array([0.0, 0.1, 0.0, 20.0, 20.1, 20.0])
    species = np.array(["a", "a", "a", "b", "b", "b"], dtype=object)

    graph = species_conditioned_knn_graph(lat, lon, species, k=1)

    assert len(graph.edges) == len(graph.edge_species) == len(graph.edge_distance_km)
    for edge, edge_sp in zip(graph.edges, graph.edge_species, strict=True):
        assert species[edge[0]] == edge_sp
        assert species[edge[1]] == edge_sp


def test_equal_species_global_statistic_gives_one_vote_per_species_not_per_edge():
    # Species a has one edge with Q=0; species b has three edges with Q=1.
    # A pooled edge mean would be .75, but the equal-species mean must be .5.
    edge_scores = np.array([0.0, 1.0, 1.0, 1.0])
    edge_species = np.array(["a", "b", "b", "b"], dtype=object)

    species_order, q, global_q = equal_species_graph_discontinuity(
        edge_scores,
        edge_species,
        expected_species=np.array(["a", "b"], dtype=object),
    )

    np.testing.assert_array_equal(species_order, np.array(["a", "b"], dtype=object))
    np.testing.assert_allclose(q, np.array([0.0, 1.0]))
    assert global_q == pytest.approx(0.5)
    assert global_q != pytest.approx(edge_scores.mean())


def test_equal_species_global_statistic_refuses_silent_species_dropout():
    with pytest.raises(ValueError, match="no species may disappear"):
        equal_species_graph_discontinuity(
            np.array([0.2, 0.4]),
            np.array(["a", "a"], dtype=object),
            expected_species=np.array(["a", "b"], dtype=object),
        )


def test_graph_permutation_null_keeps_species_and_graph_contract_fixed():
    species = np.array(["a", "a", "a", "b", "b", "b"], dtype=object)
    values = np.array([[0.0], [0.0], [1.0], [0.0], [1.0], [1.0]])
    graph = SpeciesGraph(
        edges=np.array([[0, 1], [1, 2], [3, 4], [4, 5]], dtype=int),
        edge_species=np.array(["a", "a", "b", "b"], dtype=object),
        edge_distance_km=np.array([1.0, 1.0, 1.0, 1.0]),
    )

    result = species_conditioned_graph_permutation_null(
        values,
        species,
        graph,
        n_permutations=50,
        rng=np.random.default_rng(20260830),
    )

    np.testing.assert_array_equal(result["species"], np.array(["a", "b"], dtype=object))
    assert result["observed_species_q"].shape == (2,)
    assert result["null_species_q"].shape == (50, 2)
    assert result["null_global_equal_species_mean"].shape == (50,)
    assert result["observed_global_equal_species_mean"] == pytest.approx(
        result["observed_species_q"].mean()
    )


def test_graph_null_rejects_cross_species_edge_even_if_edge_species_tag_is_forged():
    species = np.array(["a", "a", "b", "b"], dtype=object)
    values = np.array([[0.0], [1.0], [0.0], [1.0]])
    graph = SpeciesGraph(
        edges=np.array([[0, 2], [1, 0], [2, 3]], dtype=int),
        edge_species=np.array(["a", "a", "b"], dtype=object),
        edge_distance_km=np.ones(3),
    )

    with pytest.raises(ValueError, match="crosses species"):
        species_conditioned_graph_permutation_null(
            values,
            species,
            graph,
            n_permutations=2,
            rng=np.random.default_rng(1),
        )
