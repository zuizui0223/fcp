from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/acquisition/run_global_monte_carlo_capacity_scan.py"
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json"


def load_script_module():
    spec = importlib.util.spec_from_file_location("global_capacity_scan", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_capacity_target_rule_is_fixed_largest_first_and_colour_blind():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert contract["target_rule"]["candidate_raw_photos_per_species"] == [100, 80, 60]
    assert contract["target_rule"]["minimum_metadata_eligible_species"] == 300
    assert contract["target_rule"]["post_outcome_target_relaxation"] is False
    assert contract["query"]["observer_cap"] == 2
    assert contract["query"]["request_retries"] == 0
    assert contract["outcome_firewall"]["flower_colour_used"] is False
    assert contract["actual_image_acquisition_boundary"]["capacity_scan_rows_are_not_image_measurement_candidates"] is True


def test_maximum_span_is_zero_for_singleton_and_large_for_antipodal_points():
    module = load_script_module()
    singleton = pd.DataFrame({"latitude": [0.0], "longitude": [0.0]})
    assert module.maximum_span_km(singleton) == 0.0
    two = pd.DataFrame({"latitude": [0.0, 0.0], "longitude": [0.0, 180.0]})
    assert module.maximum_span_km(two) > 20000


def test_row_id_digest_is_order_invariant():
    module = load_script_module()
    a = pd.DataFrame({"observation_id": [2, 1], "photo_id": [20, 10]})
    b = pd.DataFrame({"observation_id": [1, 2], "photo_id": [10, 20]})
    assert module.row_id_digest(a) == module.row_id_digest(b)


def test_capacity_scan_excludes_prior_experiment_ids_before_candidate_acquisition():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    sources = contract["prior_experiment_exclusion_sources"]
    assert "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv" in sources
    assert "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv" in sources
    assert contract["actual_image_acquisition_boundary"]["actual_candidate_acquisition_is_a_separate_fresh_frozen_draw"] is True
