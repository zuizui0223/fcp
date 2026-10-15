"""Source-only tests for no-cap 128-worker photographic opportunity capacity planning."""
import csv
import hashlib
import json
import math

import pytest

from scripts.acquisition.plan_fcp_real_42111_partition_capacity import (
    FRAME_SHA, SOURCE_DATE, plan, partition_for,
)


def tiny_real_shaped_source(tmp_path):
    counts = tmp_path / "real_source_counts.csv"
    with counts.open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"])
        wr.writerows([
            [11, "Species alba", 65, SOURCE_DATE],
            [12, "Species rubra", 0, SOURCE_DATE],
            [13, "Species lutea", 25, SOURCE_DATE],
            [14, "Species flava", 105, SOURCE_DATE],
        ])
    digest = hashlib.sha256(counts.read_bytes()).hexdigest()
    report = {
        "schema": "fcp_42111_original_photo_opportunity_audit_v1",
        "status": "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
        "source_frame_sha256": FRAME_SHA,
        "snapshot_date": SOURCE_DATE,
        "count_ledger_sha256": digest,
        "original_source_taxon_IDs_and_names_independently_matched": True,
        "all_source_photo_id_manifests_frozen": False,
        "all_original_image_pixels_measured": False,
        "W_D_H2_or_climate_soil_results": False,
        "source_frame_species": 4,
        "all_source_photo_opportunity_ids": 195,
        "taxa_with_zero_photos": 1,
        "minimum_source_photo_depth_upper_bound_for_W": {"40": 2}
    }
    audit = tmp_path / "independent_audit.json"
    audit.write_text(json.dumps(report), encoding="utf-8")
    return counts, audit


def test_original_partition_hash_matches_frozen_fcp_hash_function():
    for tid, species in [(11, "Species alba"), (12, "Species rubra"),
                         (52821, "Achillea millefolium")]:
        h = hashlib.sha256(f"FCP_ALL_OBSERVABLE_20261014|{tid}|{species}".encode()).hexdigest()
        assert partition_for(tid, species, 128) == int(h[:16], 16) % 128


def test_real_shaped_counts_keep_all_photos_and_zero_source_taxa(tmp_path):
    counts, audit = tiny_real_shaped_source(tmp_path)
    rows, receipt = plan(counts, audit, _synthetic_test=True)
    assert len(rows) == 128
    assert sum(r["source_species"] for r in rows) == 4
    assert sum(r["source_photo_ids"] for r in rows) == 195
    assert sum(r["species_with_no_source_photos"] for r in rows) == 1
    assert sum(r["source_species_ge_40_photos"] for r in rows) == 2
    assert all(r["estimated_64_photo_work_units"] == math.ceil(r["source_photo_ids"] / 64) for r in rows)
    assert receipt["no_species_or_photo_inclusion_cap"] is True
    assert receipt["source_photo_id_manifests_complete"] is False
    assert receipt["actual_original_ROI_v4_photo_measurements_completed"] is False
    assert receipt["H2_W_999world_or_ecological_results"] is False
    assert receipt["photos"] == 195


def test_rejects_tampered_count_ledger(tmp_path):
    counts, audit = tiny_real_shaped_source(tmp_path)
    counts.write_text(counts.read_text().replace("Species rubra", "Another species"))
    with pytest.raises(ValueError, match="independent original-taxon audit"):
        plan(counts, audit, _synthetic_test=True)


def test_rejects_false_full_pixel_measurement_in_audit(tmp_path):
    counts, audit = tiny_real_shaped_source(tmp_path)
    doc = json.loads(audit.read_text())
    doc["all_original_image_pixels_measured"] = True
    audit.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="independent original-taxon audit"):
        plan(counts, audit, _synthetic_test=True)


def test_rejects_wrong_photo_count_totals_even_if_digest_matches(tmp_path):
    counts, audit = tiny_real_shaped_source(tmp_path)
    doc = json.loads(audit.read_text())
    doc["all_source_photo_opportunity_ids"] = 196
    audit.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="did not reconcile"):
        plan(counts, audit, _synthetic_test=True)
