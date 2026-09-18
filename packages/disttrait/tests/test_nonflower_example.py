import importlib.util
from pathlib import Path


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "nonflower_binary_trait.py"


def _load_example():
    spec = importlib.util.spec_from_file_location("disttrait_nonflower_example", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_nonflower_signal_demo_recovers_induced_species_level_structure() -> None:
    module = _load_example()
    out = module.run_signal_demo()
    assert out["n_species"] == 80
    assert out["rho_diversity_spatial"] > 0.70
    assert out["rho_latent_strength_spatial"] > 0.70


def test_species_conditioning_rejects_naive_between_species_geographic_confounding() -> None:
    module = _load_example()
    out = module.run_pooled_confounding_demo()
    assert out["n_species"] == 40
    assert out["naive_pooled_rho"] > 0.40
    assert abs(out["species_conditioned_mean_rho"]) < 0.08
    assert abs(out["species_conditioned_median_rho"]) < 0.05
