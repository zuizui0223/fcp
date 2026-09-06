from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/acquisition/freeze_global_monte_carlo_species_discovery.py"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_species_discovery_contract_v1.json"


def load_script_module():
    spec = importlib.util.spec_from_file_location("global_species_discovery", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_discovery_contract_has_fixed_nonadaptive_request_count():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    discovery = contract["fresh_discovery"]
    assert discovery["rounds"] == 20
    assert discovery["grid"]["cells"] == 162
    assert discovery["total_fresh_request_attempts"] == 20 * 162
    assert discovery["early_stopping"] is False
    assert discovery["replacement_after_empty_or_failed_request"] is False
    assert discovery["request_retries"] == 0
    assert contract["outcome_firewall"]["flower_colour_used"] is False


def test_jaccard_and_compact_frame_are_outcome_free():
    module = load_script_module()
    assert module.jaccard({"a", "b"}, {"b", "c"}) == 1 / 3
    frame = pd.DataFrame(
        {
            "cell_id": [1, 2],
            "observation_id": [10, 11],
            "photo_id": [20, 21],
            "species": ["Alpha beta", "Gamma delta"],
            "inat_taxon_id": [30, 31],
            "photo_url_large": ["unused-a", "unused-b"],
            "latitude": [1.0, 2.0],
            "longitude": [3.0, 4.0],
        }
    )
    compact = module.compact_frame(frame, round_id=7)
    assert compact.columns.tolist() == [
        "round_id",
        "cell_id",
        "observation_id",
        "photo_id",
        "species",
        "inat_taxon_id",
    ]
    assert compact["round_id"].tolist() == [7, 7]
    assert not any("colour" in column.casefold() for column in compact.columns)


def test_failure_policy_is_fail_closed_not_rerun_for_favourable_species_set():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    failure = contract["failure_policy"]
    assert failure["maximum_error_fraction_for_coverage_claim"] == 0.05
    assert failure["rerun_to_get_a_more_favourable_random_species_set"] is False
