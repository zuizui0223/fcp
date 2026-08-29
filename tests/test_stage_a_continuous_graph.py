import importlib.util
from pathlib import Path

import numpy as np


def load_stage_a_module():
    path = Path(__file__).parents[1] / "scripts" / "analysis" / "run_jbi_ch1_stage_a_continuous_graph.py"
    spec = importlib.util.spec_from_file_location("stage_a", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_run_k_supports_species_specific_vector_dimensions_and_equal_species_global_mean():
    module = load_stage_a_module()
    rows = []
    coords = {}
    representation = {"per_species": {}}

    for j, species in enumerate(module.EXPECTED_SPECIES):
        dim = 3 if species == "Ipomoea purpurea" else 2
        representation["per_species"][species] = {"feature_names": [f"f{i}" for i in range(dim)]}
        for i in range(6):
            photo_id = f"{j}-{i}"
            vector = [float((i + axis) % 3) for axis in range(dim)]
            rows.append(
                {
                    "species": species,
                    "photo_id": photo_id,
                    "feature_status": "ok",
                    "continuous_colour_vector_z": vector,
                }
            )
            coords[photo_id] = (species, float(j * 10) + i * 0.01, float(j * 10) + i * 0.01)

    result, null_matrix, species_order = module.run_k(
        rows,
        coords,
        representation,
        k=1,
        n_permutations=20,
        seed=123,
        min_rows=6,
    )

    assert species_order == module.EXPECTED_SPECIES
    assert null_matrix.shape == (20, 7)  # global + six species
    assert set(result["species"]) == set(module.EXPECTED_SPECIES)
    observed_species = [result["species"][sp]["q"]["observed"] for sp in module.EXPECTED_SPECIES]
    assert result["global_equal_species_mean_q"]["observed"] == np.mean(observed_species)
    assert all(result["species"][sp]["n_measurement_evaluable"] == 6 for sp in module.EXPECTED_SPECIES)


def test_mc_summary_lower_tail_is_monte_carlo_corrected():
    module = load_stage_a_module()
    summary = module.mc_summary(0.0, np.array([1.0, 1.5, 2.0, 2.5]))
    assert summary["p_lower_tail"] == 1 / 5
    assert summary["clustering_deficit"] > 0
