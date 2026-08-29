import importlib.util
from pathlib import Path

import numpy as np

PATH = Path(__file__).parents[1] / "scripts" / "data" / "analyze_jbi_ch1_calibration_feature_geometry.py"
SPEC = importlib.util.spec_from_file_location("calibration_feature_geometry", PATH)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_raphanus_uses_two_visual_pigment_proxies():
    row = {
        "visual_colour_axes": {
            "white_signal": 0.4,
            "anthocyanin_like_signal": 0.3,
            "carotenoid_like_signal": 0.2,
            "blue_purple_signal": 0.1,
        }
    }
    assert mod.feature_vector(row, "Raphanus sativus") == [0.3, 0.2]


def test_binary_species_uses_predeclared_candidate_scores():
    row = {"candidate_scores": {"blue": 0.8, "red": 0.2}}
    assert mod.feature_vector(row, "Lysimachia arvensis") == [0.8, 0.2]


def test_clear_two_cluster_synthetic_data_prefers_two_components():
    rng = np.random.default_rng(12)
    a = rng.normal(loc=(-4, -4), scale=0.25, size=(40, 2))
    b = rng.normal(loc=(4, 4), scale=0.25, size=(40, 2))
    x = np.vstack([a, b])
    grid = mod.fit_bic_grid(x, max_components=2, seed=123)
    assert mod.best_bic_components(grid) == 2


def test_bootstrap_support_is_a_probability_distribution():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(30, 2))
    support = mod.bootstrap_component_support(x, max_components=2, n_bootstraps=10, seed=9)
    assert set(support) == {"1", "2"}
    assert abs(sum(support.values()) - 1.0) < 1e-12
