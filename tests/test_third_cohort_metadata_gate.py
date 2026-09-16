import pandas as pd
import pytest

from fcp_pipeline.third_cohort_metadata_gate import (
    metadata_gate_verdict,
    seal_exact_denominator,
    validate_fresh_rows,
    validate_selected_manifest,
)


def test_manifest_requires_exact_500_ranked_unique_species():
    ok = pd.DataFrame(
        {
            "prospective_rank": range(1, 501),
            "inat_taxon_id": range(1001, 1501),
            "species": [f"Species {i}" for i in range(500)],
        }
    )
    validate_selected_manifest(ok, expected_n=500)
    bad = ok.copy()
    bad.loc[499, "prospective_rank"] = 499
    with pytest.raises(RuntimeError, match="prospective_rank"):
        validate_selected_manifest(bad, expected_n=500)


def test_validate_fresh_rows_rejects_prior_or_duplicate_ids():
    rows = pd.DataFrame({"observation_id": [1, 2], "photo_id": [11, 12]})
    with pytest.raises(RuntimeError, match="prior experiment"):
        validate_fresh_rows(rows, {2}, set())
    dup = pd.DataFrame({"observation_id": [1, 1], "photo_id": [11, 12]})
    with pytest.raises(RuntimeError, match="duplicate observation"):
        validate_fresh_rows(dup, set(), set())


def test_seal_exact_denominator_keeps_only_full100_species_and_exact_rows():
    audit = pd.DataFrame(
        {
            "prospective_rank": [1, 2, 3],
            "inat_taxon_id": [101, 102, 103],
            "species": ["A a", "B b", "C c"],
            "full_fixed_n": [True, False, True],
            "retained": [100, 98, 100],
            "request_error": ["", "", ""],
        }
    )
    rows = pd.DataFrame(
        {
            "query_species": ["A a"] * 100 + ["B b"] * 98 + ["C c"] * 100,
            "inat_taxon_id": [101] * 100 + [102] * 98 + [103] * 100,
            "observation_id": range(1, 299),
            "photo_id": range(1001, 1299),
        }
    )
    species, authorized = seal_exact_denominator(audit, rows, target_n=100)
    assert species["inat_taxon_id"].tolist() == [101, 103]
    assert len(authorized) == 200
    assert authorized.groupby("inat_taxon_id").size().to_dict() == {101: 100, 103: 100}
    assert 102 not in set(authorized["inat_taxon_id"])


def test_metadata_gate_verdict_is_fail_closed():
    assert metadata_gate_verdict(0, 500, 420, min_full_species=300) == (
        "THIRD_COHORT_METADATA_GATE_PASS",
        True,
    )
    assert metadata_gate_verdict(26, 500, 420, min_full_species=300) == (
        "THIRD_COHORT_METADATA_TRANSPORT_NOT_EVALUABLE",
        False,
    )
    assert metadata_gate_verdict(0, 500, 299, min_full_species=300) == (
        "THIRD_COHORT_METADATA_CAPACITY_NOT_EVALUABLE",
        False,
    )
