import importlib.util
import json
from pathlib import Path

import pytest


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "nonflower_binary_trait.py"
ROOT = Path(__file__).resolve().parents[3]
RECEIPT = ROOT / "results" / "disttrait_nonflower_synthetic_v0_2_20260918" / "result.json"


def _load_example():
    spec = importlib.util.spec_from_file_location("disttrait_nonflower_example", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_nonflower_signal_demo_matches_frozen_receipt() -> None:
    module = _load_example()
    out = module.run_signal_demo()
    frozen = json.loads(RECEIPT.read_text(encoding="utf-8"))["signal_demo"]
    assert out["n_species"] == frozen["species"] == 80
    assert out["rho_diversity_spatial"] == pytest.approx(
        frozen["rho_diversity_spatial"], abs=1e-12, rel=0
    )
    assert out["rho_latent_strength_spatial"] == pytest.approx(
        frozen["rho_latent_strength_spatial"], abs=1e-12, rel=0
    )
    assert out["rho_diversity_spatial"] > 0.70
    assert out["rho_latent_strength_spatial"] > 0.70


def test_species_conditioning_matches_frozen_confounding_receipt() -> None:
    module = _load_example()
    out = module.run_pooled_confounding_demo()
    frozen = json.loads(RECEIPT.read_text(encoding="utf-8"))["pooled_confounding_demo"]
    assert out["n_species"] == frozen["species"] == 40
    assert out["naive_pooled_rho"] == pytest.approx(
        frozen["naive_pooled_rho"], abs=1e-12, rel=0
    )
    assert out["species_conditioned_mean_rho"] == pytest.approx(
        frozen["species_conditioned_mean_rho"], abs=1e-12, rel=0
    )
    assert out["species_conditioned_median_rho"] == pytest.approx(
        frozen["species_conditioned_median_rho"], abs=1e-12, rel=0
    )
    assert out["naive_pooled_rho"] > 0.40
    assert abs(out["species_conditioned_mean_rho"]) < 0.08
    assert abs(out["species_conditioned_median_rho"]) < 0.05
