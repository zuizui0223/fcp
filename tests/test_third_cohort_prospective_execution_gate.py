import pandas as pd
import pytest

from fcp_pipeline.third_cohort_prospective_execution_gate import (
    AUTHORIZED_METADATA_SHA256,
    AUTHORIZED_SPECIES_SHA256,
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    MIN_H2_VECTOR_SPECIES,
    MIN_MEASUREMENT_EVALUABLE_SPECIES,
    build_start_record,
    measurement_batch,
    measurement_id,
    validate_authorized_denominator,
    validate_candidate_metadata,
    validate_stage_transition,
)


def candidate_receipt():
    return {
        "analysis": "polymorphism_h2_third_cohort_candidate_metadata_acquisition",
        "status": "complete_fresh_metadata_draw_before_pixels",
        "execution_head": "f8ee9f65b88e5608a07cd9008f34e48a09abff88",
        "decision": {
            "verdict": "THIRD_COHORT_METADATA_GATE_PASS",
            "pixels_may_be_authorized_in_separate_next_step": True,
            "biological_pixel_opening_authorized_by_this_gate": False,
            "no_species_replacement": True,
            "no_target_relaxation": True,
            "one_bounded_metadata_draw": True,
        },
        "capacity": {
            "authorized_species": 499,
            "authorized_metadata_rows": 49900,
        },
        "transport": {"request_errors": 0},
        "files": {
            "authorized_metadata_sha256": AUTHORIZED_METADATA_SHA256,
            "authorized_species_sha256": AUTHORIZED_SPECIES_SHA256,
        },
        "outcome_firewall": {
            "image_pixels_opened": False,
            "flower_colour_opened": False,
            "morph_opened": False,
            "palette_opened": False,
            "D_opened": False,
            "H2_vectors_opened": False,
            "H2_W_opened": False,
            "structured_null_opened": False,
            "P500_recovered_H2_read": False,
        },
    }


def test_candidate_metadata_contract_is_exact_and_preoutcome():
    assert EXPECTED_SPECIES == 499
    assert EXPECTED_ROWS == 49900
    validate_candidate_metadata(candidate_receipt())
    bad = candidate_receipt()
    bad["files"]["authorized_metadata_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="authorized metadata hash"):
        validate_candidate_metadata(bad)


def test_start_record_opens_nothing_before_authorization():
    r = build_start_record(
        branch="analysis/h2-third-cohort-preopening-20260916",
        head_sha="1" * 40,
        known_output_paths_absent=True,
        working_tree_inputs_verified=True,
    )
    assert r["stage"] == "PREOPENING"
    assert r["species"] == 499 and r["rows"] == 49900
    assert r["image_pixels_opened_before_start"] is False
    assert r["H2_opened_before_start"] is False
    assert r["opening_authorized"] is True


def test_stage_contract_reaches_h2_complete_only_from_passed_support():
    support = {
        "stage": "SUPPORT_GATE_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "measurement_evaluable_species": MIN_MEASUREMENT_EVALUABLE_SPECIES,
        "support_decision": "PASS",
    }
    h2 = {
        "stage": "H2_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "primary_h2_vector_species": MIN_H2_VECTOR_SPECIES,
        "primary_structured_null_p": 0.049,
        "h2_decision": "CONFIRMED",
    }
    validate_stage_transition(support, h2)
    bad = dict(h2)
    bad["h2_decision"] = "NOT_CONFIRMED"
    with pytest.raises(RuntimeError, match="H2 decision"):
        validate_stage_transition(support, bad)


def test_blind_ids_are_deterministic_and_batch_only_0_or_1():
    a = measurement_id(12345)
    assert a == measurement_id(12345)
    assert a != measurement_id(12346)
    assert a.startswith("FCPH2T3-")
    assert measurement_batch(a) in {0, 1}


def test_authorized_denominator_requires_exact_species_sets_counts_and_unique_ids():
    species = pd.DataFrame({
        "species": ["A a", "B b"],
        "inat_taxon_id": [1, 2],
    })
    metadata = pd.DataFrame({
        "query_species": ["A a", "A a", "B b", "B b"],
        "inat_taxon_id": [1, 1, 2, 2],
        "observation_id": [11, 12, 21, 22],
        "photo_id": [101, 102, 201, 202],
        "photo_url_large": ["https://x/a.jpg"] * 4,
        "photo_license": ["cc-by"] * 4,
        "latitude": [1.0, 2.0, 3.0, 4.0],
        "longitude": [10.0, 20.0, 30.0, 40.0],
    })
    normalized = validate_authorized_denominator(
        metadata,
        species,
        expected_species=2,
        target_per_species=2,
    )
    assert len(normalized) == 4
    assert normalized["species"].tolist() == ["A a", "A a", "B b", "B b"]

    bad = metadata.copy()
    bad.loc[3, "photo_id"] = 201
    with pytest.raises(RuntimeError, match="duplicate photo"):
        validate_authorized_denominator(bad, species, expected_species=2, target_per_species=2)
