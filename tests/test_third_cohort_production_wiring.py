from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIREWALL = ROOT / "scripts" / "analysis" / "build_polymorphism_h2_third_cohort_measurement_firewall_20260917.py"
REASSEMBLER = ROOT / "scripts" / "analysis" / "reassemble_polymorphism_h2_third_cohort_measurement_20260917.py"
H2 = ROOT / "scripts" / "analysis" / "run_polymorphism_h2_third_cohort_prospective_white_axis_20260917.py"


def test_third_cohort_production_scripts_exist_and_do_not_import_p500_gate():
    for path in (FIREWALL, REASSEMBLER, H2):
        assert path.exists(), path
        text = path.read_text(encoding="utf-8")
        assert "p500_prospective_execution_gate" not in text
        assert "third_cohort_prospective_execution_gate" in text
        assert "POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md" in text


def test_h2_executor_keeps_frozen_target_null_and_seeds():
    text = H2.read_text(encoding="utf-8")
    for needle in (
        "N_NULL = 999",
        '"primary_0_10": 20260915',
        '"strict_0_20": 20261015',
        '"primary_0_10": 0.10',
        '"strict_0_20": 0.20',
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
        "H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED",
        "H2_PROSPECTIVE_NOT_EVALUABLE_MEASUREMENT_SUPPORT",
        "H2_PROSPECTIVE_NOT_EVALUABLE_VECTOR_SUPPORT",
    ):
        assert needle in text


def test_reassembler_writes_durable_support_stage_receipt():
    text = REASSEMBLER.read_text(encoding="utf-8")
    for needle in (
        '"schema": "third_cohort_prospective_measurement_result_v1"',
        '"stage": "SUPPORT_GATE_COMPLETE"',
        '"persisted_image_pixels": False',
        '"replacement_rows": 0',
        '"replacement_species": 0',
        '"H2_opened": False',
    ):
        assert needle in text
