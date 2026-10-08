"""Synthetic regression tests for post hoc white/chromatic gradients."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE = Path(__file__).resolve().parents[1] / "scripts/analysis/run_polymorphism_white_chromatic_clines_20261008.py"
spec = importlib.util.spec_from_file_location("fcp_chromatic_clines", SOURCE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_standardized_direction_and_label_permutation():
    rng = np.random.default_rng(10)
    x = np.r_[np.linspace(10, 12, 20), np.linspace(30, 32, 20)]
    morph = np.r_[np.zeros(20, dtype=bool), np.ones(20, dtype=bool)]
    obs, null = mod.species_delta(x, morph, 1.0, rng=rng, n_perm=199)
    assert obs > 1.0
    assert len(null) == 199
    assert abs(null.mean()) < 0.3
    assert float((null >= obs).mean()) < 0.02


def test_insufficient_morphs_and_axis_span_fail_closed():
    rng = np.random.default_rng(1)
    assert mod.species_delta(np.r_[np.zeros(35), np.ones(4)],
                             np.r_[np.zeros(35, bool), np.ones(4, bool)], 0,
                             rng=rng) is None
    assert mod.species_delta(np.zeros(40), np.r_[np.zeros(20, bool), np.ones(20, bool)],
                             1, rng=rng) is None


def test_local_geographic_segregation_and_cooccurrence():
    # One pure white site and one pure coloured site, >50 km apart.
    d = pd.DataFrame({
        "latitude": np.r_[np.full(20, 10.0), np.full(20, 12.0)],
        "longitude": np.zeros(40),
        "morph": ["white"]*20 + ["red_pink"]*20,
    })
    result = mod.local_white_chromatic(d)
    assert result["eligible"]
    assert result["white_chromatic_ratio"] == 0
    assert not result["has_local_white_chromatic_pair"]
    # Each site now has ten of each, so white/colour pairs are locally possible.
    d["morph"] = (["white"]*10+["red_pink"]*10)*2
    mixed = mod.local_white_chromatic(d)
    assert mixed["white_chromatic_ratio"] > 0.9
    assert mixed["has_local_white_chromatic_pair"]


def test_holm_two_axis_and_residualization():
    p_lat, p_alt = mod.holm_two(0.01, 0.03)
    assert np.isclose(p_lat, 0.02)
    assert np.isclose(p_alt, 0.03)
    assert mod.holm_two(0.03, 0.01) == (0.03, 0.02)
    x = np.arange(40, dtype=float)
    y = 20 + 3*x + np.sin(x)
    resid = mod.residualize(y, x)
    assert np.isclose(np.corrcoef(resid, x)[0, 1], 0.0, atol=1e-10)
