"""Do not generate final-FCP biological prose from source-only or partial CI success."""
from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

import pytest

from scripts.analysis import build_fcp_u100_final_publication_packet_20261017 as m


def _row(n: int, effect: float, p: float = 0.005) -> dict:
    if n == 0:
        return {"evaluable": False, "n_species": 0}
    return {
        "evaluable": True,
        "radius_km": 50,
        "n_species": n,
        "mean_D_pair": 0.25,
        "mean_D_local": 0.25 - effect,
        "mean_depletion": effect,
        "mean_depletion_p_upper": p,
    }


def _fixture(tmp_path: Path):
    photos = tmp_path / "full_447500_measured.csv.gz"
    photos.write_bytes(b"synthetic fixture for hashes, NOT actual photo data\n")
    rows = {
        name: _row(n, effect, p)
        for name, (n, effect, p) in m.HISTORICAL.items()
    }
    rows["new_U100"] = _row(400, 0.0135)
    weighted = sum(r["n_species"] * r["mean_depletion"] for r in rows.values())
    weighted /= sum(r["n_species"] for r in rows.values())
    pooled = _row(sum(r["n_species"] for r in rows.values()), weighted)
    joined = {
        "schema": "fcp_all4730_original_149900_plus_all_new_pixel_measurements_v1",
        "total_original_100photo_measured_species": 4475,
        "total_original_pixel_terminal_photo_rows": 447500,
        "original_149900_photo_rows_unchanged": True,
        "old_50km_and_W_estimators_unchanged": True,
        "same_original_iNaturalist_source_not_independent_provider": True,
        "new_U100_source_reuse_excluded_taxa": 3,
        "new_U100_source_reuse_excluded_taxon_IDs": m.EXCLUDED_TAXON_IDS[:],
        "new_U100_source_reuse_excluded_original_photo_rows": 300,
        "remaining_U100_under_100_original_H9_photos": 252,
        "source_sha256": m.ORIGINAL_SOURCE_SHAS.copy(),
        "merged_photo_table_sha256": m.digest(photos),
        "cohorts": [
            {"cohort": name, "named_species": spp, "photo_rows": 100*spp,
             "four_state_classifiable_rows": spp * 60}
            for name, spp in zip(m.COHORTS, (500, 500, 499, 498, 2478))
        ],
    }
    spatial = {
        "schema": "fcp_all4730_same_original_50km_vertex_null_v1",
        "old_three_original_149900_effects_exactly_reproduced": True,
        "original_50km_radius_changed": False,
        "original_species_equal_weighting_changed": False,
        "original_199_vertex_null_changed": False,
        "still_posthoc_50km": True,
        "all_sources": {name: {"50": row} for name, row in rows.items()},
        "pooled_original_50km": pooled,
        "primary_50km_evaluable_species": pooled["n_species"],
        "primary_50km_evaluable_species_by_cohort":
            {name: row["n_species"] for name, row in rows.items()},
    }
    h1w = {
        "schema": "fcp_all4730_unchanged_historical_H1_D_H2_W_v1",
        "all_real_photo_rows": 447500,
        "existing_method_called_directly": True,
        "old_frozen_149900_inference_not_overwritten": True,
        "same_iNaturalist_source": True,
        "H1_original_observer_disjoint": {"primary": {"some_reliability_value": 0.8}},
        "H2_W_original_fixed_white_axis": {"primary_0_10": {"observed_W": 0.5}},
    }
    paths = {}
    for name, document in (("combine", joined), ("spatial", spatial), ("h1w", h1w)):
        path = tmp_path / (name + ".json")
        path.write_text(json.dumps(document), encoding="utf-8")
        paths[name] = path
    paths["photos"] = photos
    receipt = {
        "schema": "fcp_all4730_permanent_measured_original_photo_results_v1",
        "final_run_id": 39000000000,
        "full_original_eligible_species": 4730,
        "genuinely_measured_100photo_species": 4475,
        "genuinely_measured_photo_rows": 447500,
        "actual_100photo_species_analysis_omissions_total": 255,
        "old_149900_original_photos_unchanged": True,
        "original_149900_50km_three_effects_exactly_reproduced": True,
        "no_new_ecological_statistic_or_classifier": True,
        "source_archive_sha256": m.digest(photos),
        "source_join_receipt_sha256": m.digest(paths["combine"]),
        "original_spatial_result_sha256": m.digest(paths["spatial"]),
        "original_H1_D_W_result_sha256": m.digest(paths["h1w"]),
    }
    return receipt, joined, spatial, h1w, paths


def _check(data):
    return m.verify_packet(*data)


def _rehash(receipt, paths):
    """For tests about statistical contracts, make a newly tampered file SHA-consistent."""
    for key, receipt_key in (
        ("spatial", "original_spatial_result_sha256"),
        ("combine", "source_join_receipt_sha256"),
        ("h1w", "original_H1_D_W_result_sha256"),
    ):
        receipt[receipt_key] = m.digest(paths[key])


def test_valid_full_five_cohort_packet_without_claiming_h1w_support(tmp_path):
    data = _fixture(tmp_path)
    packet = _check(data)
    assert packet["actual_measured_species"] == 4475
    assert packet["actual_measured_photos"] == 447500
    assert packet["primary_spatial_species"] == 1141
    assert packet["cohorts"]["wave4"]["n_species"] == 190
    assert packet["cohorts"]["new_U100"]["n_species"] == 400
    assert packet["pooled_50km_positive_under_null"] is True
    assert packet["h1w_requires_separate_reporting_of_actual_full_numeric_results"] is True
    text = m.manuscript_text(packet)
    assert "4,475 distinct species" in text
    assert "| Fifth cohort (new U100) | 2,478 | 247,800 | 400 |" in text
    assert "post hoc 50-km" in text
    assert "not natural-population morph-frequency estimates" in text
    assert "no H1/W support verdict is inferred" in text


@pytest.mark.parametrize("key", [
    "source_archive_sha256",
    "source_join_receipt_sha256",
    "original_spatial_result_sha256",
    "original_H1_D_W_result_sha256",
])
def test_fail_closed_on_any_unverified_source_byte(tmp_path, key):
    data = _fixture(tmp_path)
    data[0][key] = "0" * 64
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        _check(data)


def test_incomplete_source_only_receipt_cannot_make_publication_packet(tmp_path):
    data = _fixture(tmp_path)
    data[0]["schema"] = "fcp_u100_full_original_source_only_v1"
    with pytest.raises(ValueError, match="not a genuine complete"):
        _check(data)
    data = _fixture(tmp_path)
    data[0]["genuinely_measured_photo_rows"] = 199700
    with pytest.raises(ValueError, match="denominators invalid"):
        _check(data)


def test_old_frozen_wave4_50km_effect_cannot_change(tmp_path):
    receipt, joined, spatial, h1w, paths = _fixture(tmp_path)
    spatial["all_sources"]["wave4"]["50"] = _row(190, 0.025)
    with pytest.raises(ValueError, match="previously frozen spatial result changed"):
        m.verify_packet(receipt, joined, spatial, h1w, paths)


def test_reject_mismatched_local_pair_denominator_even_when_five_source_counts_match(tmp_path):
    data = _fixture(tmp_path)
    data[2]["primary_50km_evaluable_species"] += 3
    with pytest.raises(ValueError, match="pooled true local-pair eligible n"):
        _check(data)


def test_reject_tampered_null_or_radius(tmp_path):
    for change in ("radius", "p"):
        data = _fixture(tmp_path)
        row = data[2]["all_sources"]["new_U100"]["50"]
        if change == "radius":
            row["radius_km"] = 100
            with pytest.raises(ValueError, match="invalid 50-km test"):
                _check(data)
        else:
            row["mean_depletion_p_upper"] = 0.004
            with pytest.raises(ValueError, match="plus-one p"):
                _check(data)


def test_a_negative_final_effect_is_explicitly_reported_as_unsupported(tmp_path):
    data = _fixture(tmp_path)
    new = data[2]["all_sources"]["new_U100"]["50"]
    new.update(_row(400, -0.10, 1.0))
    rows = [c["50"] for c in data[2]["all_sources"].values()]
    new_effect = sum(x["n_species"] * x["mean_depletion"] for x in rows) / 1141
    data[2]["pooled_original_50km"] = _row(1141, new_effect, 1.0)
    packet = _check(data)
    assert packet["pooled_50km_positive_under_null"] is False
    assert "not supported" in m.manuscript_text(packet)
    assert "not independent-provider" not in m.manuscript_text(packet)


def test_missing_full_h1w_is_not_equivalent_to_a_successful_workflow(tmp_path):
    data = _fixture(tmp_path)
    data[3]["H2_W_original_fixed_white_axis"] = {}
    with pytest.raises(ValueError, match="H1/D/W did not complete"):
        _check(data)


def test_cli_writes_only_from_sha_authenticated_inputs(tmp_path, monkeypatch):
    receipt, joined, spatial, h1w, paths = _fixture(tmp_path)
    recpath = tmp_path / "receipt.json"
    recpath.write_text(json.dumps(receipt), encoding="utf-8")
    out = tmp_path / "packet"
    monkeypatch.setattr(sys, "argv", [
        "packet", "--receipt", str(recpath), "--combine", str(paths["combine"]),
        "--spatial", str(paths["spatial"]), "--h1w", str(paths["h1w"]),
        "--photos", str(paths["photos"]), "--out", str(out),
    ])
    m.main()
    assert (out / "fcp_all4730_verified_manuscript_results_packet.md").exists()
    result = json.loads((out / "fcp_all4730_verified_manuscript_results_packet.json").read_text())
    assert result["source_inference_complete"] is True
    assert result["source_sha256"]["spatial"] == m.digest(paths["spatial"])
