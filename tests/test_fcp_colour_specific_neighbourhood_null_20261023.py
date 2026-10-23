"""Synthetic guards; do not mistake RGB cluster labels for inherited morphs."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
from fcp_colour_specific_neighbourhood_null_20261023 import (
    COLORS, CASES, Accumulator, prepare, permute_labels, self_test
)


def fake_two_colour_neighbourhoods():
    rows = []
    for i in range(40):
        is_white = i < 20
        rows.append({
            "morph": "white" if is_white else "blue_purple",
            "latitude": 0. if is_white else 20.,
            "longitude": .001 * (i % 20),
            "positional_accuracy_m": 100.,
            "observer_id": i,
            "photo_id": i + 1,
            "month_of_year": 3 if is_white else 5,
            "primary_colour_pixel_fraction": .85,
            "other_max_colour_pixel_fraction": .10,
        })
    return pd.DataFrame(rows)


def test_hue_components_sum_to_original_global_minus_local_mismatch():
    self_test()
    g = fake_two_colour_neighbourhoods()
    p = prepare(g, min_each=5)
    assert p is not None and p["n_photos"] == 40
    assert p["n_pairs"] == 380 and p["richness"] == 2
    expected = 2 * 20 * 20 / (40 * 39)
    assert abs(p["observed"].sum() - expected) < 1e-12
    assert np.isclose(p["observed"][0], expected / 2)
    assert np.isclose(p["observed"][3], expected / 2)
    assert np.isclose(p["observed"][1], 0)
    assert np.isclose(p["observed"][2], 0)


def test_month_null_preserves_season_confounded_geography():
    g = fake_two_colour_neighbourhoods()
    p = prepare(g, min_each=5)
    a = Accumulator()
    a.add(p, seed=20261023, min_each=5, month_preserving=True)
    r = a.finish("month", month_preserving=True)
    assert r["n_species"] == 1
    # All the white photographs are in March and blue in May: a month
    # label-null has no exchangeability and must not call this selection.
    assert r["overall_upper_p_199"] == 1
    assert r["n_month_permutable_species"] == 0
    assert abs(r["overall_observed_mean_global_minus_local_mismatch"] -
               r["overall_matched_null_mean"]) < 1e-12
    # Single synthetic taxon cannot establish aggregate hue-specific p.
    assert all(x["status"] == "not_estimable_less_than_20_species"
               for x in r["colour_specific_supported_purity_not_genotype"].values())


def test_unconditional_colour_null_does_not_preserve_seasonal_clusters():
    p = prepare(fake_two_colour_neighbourhoods(), min_each=5)
    a = Accumulator()
    a.add(p, seed=20261023, min_each=5, month_preserving=False)
    r = a.finish("unrestricted", month_preserving=False)
    assert r["overall_observed_mean_global_minus_local_mismatch"] > 0.5
    assert r["overall_upper_p_199"] == 0.005


def test_observer_and_dominance_qc_are_explicit_not_genetic():
    g = fake_two_colour_neighbourhoods()
    assert prepare(g, min_each=5, one_per_observer=True)["richness"] == 2
    bad = g.copy()
    bad["primary_colour_pixel_fraction"] = .40
    assert prepare(bad, min_each=5, dominant=True) is None
    assert len(COLORS) == 4
    assert len(CASES) >= 6
    assert all(x[0] != "genetically_identified_morph" for x in CASES)
