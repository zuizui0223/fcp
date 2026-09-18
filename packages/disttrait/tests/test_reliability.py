import pandas as pd

from disttrait import observer_disjoint_reliability


def _synthetic_frame() -> pd.DataFrame:
    rows = []
    patterns = {
        "sp0": ["a", "a", "a", "a"],
        "sp1": ["a", "a", "a", "b"],
        "sp2": ["a", "a", "b", "b"],
        "sp3": ["a", "b", "b", "b"],
        "sp4": ["b", "b", "b", "b"],
    }
    for species, pattern in patterns.items():
        for observer in range(4):
            for state in pattern:
                rows.append(
                    {
                        "species": species,
                        "observer": f"o{observer}",
                        "state": state,
                    }
                )
    return pd.DataFrame(rows)


def test_observer_disjoint_reliability_recovers_identical_species_ordering():
    result = observer_disjoint_reliability(
        _synthetic_frame(),
        species_col="species",
        observer_col="observer",
        state_col="state",
        states=["a", "b"],
        n_partitions=25,
        base_seed=100,
        min_full_classifiable=8,
        min_half_classifiable=4,
    )
    assert result.summary["eligible_species"] == 5
    assert result.summary["paired_n_median"] == 5
    assert result.summary["rho_median"] == 1.0
    assert result.summary["ccc_median"] == 1.0
