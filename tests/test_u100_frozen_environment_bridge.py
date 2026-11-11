"""Synthetic checks of species-matched 18-test environmental bridge."""
import numpy as np
import pandas as pd
import pytest
from fcp_pipeline.u100_frozen_environment_bridge import holm_adjust, compare_species, verified_run


def test_holm_handles_all_eighteen_including_missing():
    p = np.array([.001, .03, .06, np.nan] + [1.] * 14)
    adj = holm_adjust(p)
    assert len(adj) == 18
    assert adj[0] == pytest.approx(.018)
    assert adj[1] >= .05
    assert adj[-1] == 1


def test_requires_unique_species_and_correct_cohort():
    registry = pd.DataFrame([{"inat_taxon_id": i + 1, "species": f"S{i}", "cohort": "prospective",
                              "D": i / 5000, "w_alignment_sq": i / 5000,
                              "depletion_50km": i / 5000} for i in range(4730)])
    env = pd.DataFrame([{"species": f"S{i}", "variable": k,
                         "delta_white_minus_nonwhite_SD": (i - 30) / 100}
                        for k in ("bio5", "bio14", "srad_mean") for i in range(50)])
    result, joined = compare_species(registry, env)
    assert len(result) == 18
    assert result["n_species"].eq(50).all()
    assert len(joined) == 4730
    assert result["holm_18_p"].between(0, 1).all()
    registry.loc[0, "species"] = "S1"
    with pytest.raises(ValueError, match="duplicated species"):
        compare_species(registry, env)


def test_unpinned_environment_source_is_rejected(tmp_path):
    source = tmp_path / "fake.csv"
    source.write_text("species,variable,delta_white_minus_nonwhite_SD\nX,bio5,0.1\n")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        verified_run(source, tmp_path / "ref.json", tmp_path / "reg.csv", tmp_path / "out")
