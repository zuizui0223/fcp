import json
import os
from pathlib import Path

import pandas as pd
import pytest

from disttrait import observer_disjoint_reliability


ROOT = Path(__file__).resolve().parents[3]
STATES = ["white", "yellow_orange", "red_pink", "blue_purple"]


def _path_from_env(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        pytest.skip(f"{name} is not set; full raw-artifact replay is optional")
    path = Path(value)
    if not path.exists():
        pytest.skip(f"{name} does not exist: {path}")
    return path


def _replay(path: Path):
    frame = pd.read_csv(path, low_memory=False)
    return observer_disjoint_reliability(
        frame,
        species_col="species",
        observer_col="observer_id",
        state_col="morph",
        classifiable_col="global_classifiable",
        states=STATES,
        n_partitions=200,
        base_seed=20260913,
        min_full_classifiable=40,
        min_half_classifiable=20,
        minimum_distinct_observers=2,
    )


def test_optional_full_fcp_h1_raw_artifact_replay() -> None:
    discovery_path = _path_from_env("DISTTRAIT_FCP_DISCOVERY")
    reserve_path = _path_from_env("DISTTRAIT_FCP_RESERVE")
    expected_path = ROOT / "results" / "polymorphism_h1_observer_disjoint_reliability_20260913" / "result.json"
    expected = json.loads(expected_path.read_text(encoding="utf-8"))

    observed = {
        "discovery": _replay(discovery_path).summary,
        "reserve": _replay(reserve_path).summary,
    }
    fields = [
        "paired_n_median",
        "rho_median",
        "rho_q05",
        "rho_q95",
        "ccc_median",
        "spearman_brown_median",
        "mae_median",
    ]
    for cohort in ("discovery", "reserve"):
        frozen = expected[cohort]["primary20"]
        for field in fields:
            assert observed[cohort][field] == pytest.approx(
                frozen[field], abs=1e-12, rel=0
            )
