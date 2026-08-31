import hashlib
import json
from pathlib import Path

from scripts.data.extract_inaturalist_automated_colour_states import PROTOCOL_VERSION


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "docs/supporting/jbi_inaturalist_automated_colour_plan_manifest_v2.json"


def test_automated_colour_plan_is_fail_closed_and_location_sealed():
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    assert plan["protocol"] == PROTOCOL_VERSION
    assert plan["human_judgement_used"] is False
    assert plan["development_spatial_fields_opened"] is False
    assert plan["locked_partition_opened"] is False
    assert plan["spatial_colour_outcome_opened"] is False
    assert plan["pre_freeze_visibility"]["partial_admission_counts_used_for_retuning"] is False
    assert plan["development_gate"]["total_encounters"] == 480
    assert plan["development_gate"]["total_attached_photos"] == 886
    assert plan["locked_gate"]["encounters_per_species"] == 120
    assert plan["spatial_test"]["whole_vector_permutations_within_species"] == 9999
    assert plan["spatial_test"]["universality_claim_allowed"] is False
    assert plan["spatial_test"]["mechanism_claim_allowed"] is False


def test_automated_colour_plan_source_hashes_match_exact_files():
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    for relative_path, expected in plan["source_sha256"].items():
        observed = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert observed == expected, relative_path
