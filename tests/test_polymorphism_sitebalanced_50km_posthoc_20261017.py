"""Synthetic null/sampling-density controls; these rows are not flower data."""
import numpy as np
import pandas as pd

from scripts.analysis.run_polymorphism_sitebalanced_50km_posthoc_20261017 import (
    PERMUTATIONS, photo_cells, analyse_sitebalanced_species, summarise,
)


def sample(na=12, nb=12, nc=36, ca="white", cb="blue_purple", cc="white"):
    return pd.DataFrame({
        "inat_taxon_id": [123] * (na + nb + nc),
        "species": ["Synthetic testus"] * (na + nb + nc),
        "photo_id": [str(k).zfill(4) for k in range(na + nb + nc)],
        "morph": [ca] * na + [cb] * nb + [cc] * nc,
        "latitude": [0.0] * (na + nb + nc),
        "longitude": [0.0] * na + [0.30] * nb + [10.0] * nc,
    })


def test_balancing_blocks_not_photo_pairs():
    d = sample()
    row, nulls = analyse_sitebalanced_species(d, "discovery")
    assert row is not None
    assert row["n_all_50km_photo_pairs"] == 66 + 144 + 66 + 630
    # 0.25-degree blocks A-A, B-B, A-B, C-C have equal weight.
    assert row["n_cell_pair_blocks"] == 4
    assert np.isclose(row["D_local_sitebalanced"], 0.25)
    assert np.isclose(row["D_local_unweighted"], 144 / (66 + 144 + 66 + 630))
    assert nulls.shape == (PERMUTATIONS,)
    assert np.isfinite(nulls).all()


def test_invariance_to_photo_row_order():
    d = sample()
    a, n1 = analyse_sitebalanced_species(d, "validation")
    b, n2 = analyse_sitebalanced_species(d.sample(frac=1, random_state=19), "validation")
    assert a == b
    np.testing.assert_array_equal(n1, n2)


def test_one_geographic_cell_is_underidentified_not_monomorphic():
    d = sample(na=50, nb=0, nc=0)
    a, nulls = analyse_sitebalanced_species(d, "third")
    assert a is nulls is None


def test_require_supported_cell_pair_blocks():
    d = sample(na=39, nb=1, nc=0)
    a, nulls = analyse_sitebalanced_species(d, "third")
    assert a is nulls is None


def test_invalid_coordinates_are_rejected():
    with np.testing.assert_raises_regex(ValueError, "outside"):
        photo_cells(np.array([92.]), np.array([5.]))
    assert photo_cells([0.0], [180.0]).item() == photo_cells([0.0], [-180.0]).item()


def test_null_result_is_not_called_significant_without_floral_structure():
    d = sample(ca="white", cb="white", cc="white")
    a, n = analyse_sitebalanced_species(d, "discovery")
    assert a["depletion_sitebalanced"] == 0
    s = summarise([a], [n])
    assert s["p_upper"] == 1
    assert s["mean_sitebalanced_depletion"] == 0


def test_stability_of_species_denominator_for_paired_comparator():
    d = sample()
    a, n = analyse_sitebalanced_species(d, "discovery")
    s = summarise([a], [n])
    assert s["n_species"] == 1
    assert np.isclose(s["mean_unweighted_depletion_same_species"], a["depletion_unweighted"])
    assert s["species_bootstrap_mean_95"] == [a["depletion_sitebalanced"]] * 2
