"""Numerical and biological null-separation regression checks; no source data invented."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analysis" / "fcp_all4475_continuous_9palette_geography_20261023.py"
spec = importlib.util.spec_from_file_location("fcp_cont9", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def toy(n=24, mixed=False):
    rows = []
    for patch in (0, 1):
        for j in range(n):
            x = {"inat_taxon_id": 111, "species": "Example alpha",
                 "photo_id": patch * 100 + j, "observer_id": patch * 100 + j,
                 "calendar_month": 6,
                 "morph": mod.COARSE[0] if mixed else mod.COARSE[patch],
                 "latitude": 35. + 4 * patch, "longitude": 135. + .001 * j,
                 "positional_accuracy_m": 25.}
            x.update({col: float(i == (0 if patch == 0 else 6))
                      for i, col in enumerate(mod.FRACTIONS)})
            rows.append(x)
    return pd.DataFrame(rows)


def test_jsd_unit_and_no_missing_zero_mass():
    z = np.zeros((2, 9)); z[0, 0] = 1; z[1, 6] = 1
    m = mod.jsd_base2_matrix(z)
    assert np.allclose(m, [[0, 1], [1, 0]])
    with pytest.raises(ValueError, match="zero-mass"):
        mod.jsd_base2_matrix(np.zeros((3, 9)))


def test_correct_two_stage_nulls():
    d = toy()
    r, n = mod.taxon_result(d, 111, permutations=49)
    assert r["n_pairs_50km"] == 2 * 24 * 23 // 2
    assert r["excess_50km_within_month"] > 0.2
    assert abs(r["excess_50km_within_month_and_coarse_label"]) < 1e-12
    assert n.shape == (3, 49, 5)
    d.loc[d.morph == mod.COARSE[0], "calendar_month"] = 5
    d.loc[d.morph == mod.COARSE[1], "calendar_month"] = 6
    r, n = mod.taxon_result(d, 111, permutations=49)
    assert abs(r["excess_50km_within_month"]) < 1e-12
    assert np.ptp(n[1, :, 2]) < 1e-12


def test_within_one_coarse_label_gradient_is_detected():
    d = toy(mixed=True)
    r, _ = mod.taxon_result(d, 111, permutations=49)
    assert r["coarse_label_richness_observed"] == 1
    assert r["coarse_delta_50km"] == 0
    assert r["excess_50km_within_month_and_coarse_label"] > 0.2


def test_photo_order_only_affects_monte_carlo_world_not_statistic():
    d = toy()
    original, _ = mod.taxon_result(d, 111, permutations=19)
    shuffled, _ = mod.taxon_result(d.sample(frac=1, random_state=121), 111, permutations=19)
    assert original == shuffled


def test_eligibility_zero_pairs_and_observer_thinning():
    d = toy()
    x, _ = mod.taxon_result(d, 111, one_per_observer=True, permutations=19)
    assert x["n_after_thinning"] == 48
    d.loc[d.photo_id >= 100, "latitude"] = -12
    y, _ = mod.taxon_result(d, 111, permutations=19)
    assert y["n_pairs_50km"] > 0
    d["latitude"] = np.linspace(-85, 85, len(d))
    assert mod.taxon_result(d, 111, permutations=19) is None


def test_holm_and_nonestimable_radius():
    assert np.allclose(mod.holm_adjust([0.01, 0.02, 0.03]), [0.03, 0.04, 0.04])
    r, arr = mod.taxon_result(toy(), 111, permutations=19)
    out = mod.make_summary([(r, arr)], 19, strict=False, min_pairs=1, do_bootstrap=False)
    assert out["statistics"]["full_species_shuffle"][2]["n_species"] == 1
    assert out["statistics"]["within_month"][0]["p_upper_199"] is not None
    assert out["statistics"]["within_month"][2]["genus_bootstrap_excess_95ci"] is None
