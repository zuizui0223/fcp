"""Synthetic tests for the full 42,111-taxon ORIGINAL photo-source-count audit.

No colour images, photo measurements or provider API calls are used.
"""
import csv
import hashlib
import json

import pytest

from scripts.acquisition.summarize_fcp_real_photo_opportunities import (
    ARCHIVE_SHA, FRAME_SHA, TABLE_SHA, audit,
)


def fixture(tmp_path):
    counts = tmp_path / "photo_counts.csv"
    with counts.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"])
        w.writerows([
            [11, "Test species A", 2, "2026-09-27"],
            [12, "Test species B", 1, "2026-09-27"],
            [13, "Test species C", 0, "2026-09-27"],
            [14, "Test species D", 5, "2026-09-27"],
        ])
    digest = hashlib.sha256(counts.read_bytes()).hexdigest()
    source = {
        "schema": "fcp_full_source_photo_opportunity_counts_v1",
        "status": "COMPLETE_REAL_SOURCE_PHOTO_COUNTS_ONLY_NO_PHOTO_ID_FIREWALL",
        "source_species": 4,
        "source_species_frame_sha256": FRAME_SHA,
        "source_archive_attestation": {"source_archive_sha256": ARCHIVE_SHA},
        "source_input_sha256": TABLE_SHA,
        "snapshot_date": "2026-09-27",
        "counts_cover_every_target_taxon": True,
        "photo_cap": None,
        "counts_file_sha256": digest,
        "all_original_photo_ID_manifests_frozen": False,
        "full_source_census_complete": False,
        "full_study_original_images_measured": False,
        "historical_flower_colour_csv_used": False,
        "source_image_pixels_opened": False,
        "n_taxa_with_source_photos": 3,
        "n_taxa_without_source_photos": 1,
        "n_target_taxon_licensed_photos": 8,
    }
    execution = {
        "status": "COMPLETE_REAL_42111_PHOTO_OPPORTUNITY_COUNTS_ONLY_NOT_FULL_SOURCE",
        "target_frame_species": 4,
        "original_archive_sha256": ARCHIVE_SHA,
        "all_four_extracted_table_sha256": TABLE_SHA,
        "immutable_source_date": "2026-09-27",
        "source_photo_opportunity_counts_completed": True,
        "complete_metadata_census": False,
        "complete_original_images_measured": False,
        "new_W_or_environment_results_computed": False,
        "old_measured_flower_colours_used": False,
        "photo_id_manifest_universe_frozen": False,
        "photo_counts_csv_sha256": digest,
        "n_frozen_source_photo_opportunities": 8,
        "n_taxa_with_source_photos": 3,
        "n_taxa_without_source_photos": 1,
    }
    source_file = tmp_path / "source.json"
    execution_file = tmp_path / "execution.json"
    source_file.write_text(json.dumps(source), encoding="utf-8")
    execution_file.write_text(json.dumps(execution), encoding="utf-8")
    return counts, source_file, execution_file


def test_exact_counts_zero_taxa_and_no_fake_flower_morphs(tmp_path):
    counts, source, execution = fixture(tmp_path)
    result = audit(counts, source, execution, species=4)
    assert result["source_frame_species"] == 4
    assert result["all_source_photo_opportunity_ids"] == 8
    assert result["taxa_with_zero_photos"] == 1
    assert result["minimum_source_photo_depth_upper_bound_for_W"]["1"] == 3
    assert result["minimum_source_photo_depth_upper_bound_for_W"]["40"] == 0
    assert result["fraction_of_source_photos_in_top_one_percent_of_taxa"] == 5 / 8
    assert result["all_source_photo_id_manifests_frozen"] is False
    assert result["W_D_H2_or_climate_soil_results"] is False


def test_modified_ledger_fails_original_source_digest(tmp_path):
    counts, source, execution = fixture(tmp_path)
    counts.write_text(counts.read_text(encoding="utf-8").replace("Test species A", "Tampered"),
                      encoding="utf-8")
    with pytest.raises(ValueError, match="fingerprint"):
        audit(counts, source, execution, species=4)


def test_deceptive_claim_of_image_measurement_is_rejected(tmp_path):
    counts, source, execution = fixture(tmp_path)
    doc = json.loads(execution.read_text())
    doc["new_W_or_environment_results_computed"] = True
    execution.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(ValueError, match="biological gate"):
        audit(counts, source, execution, species=4)


def test_duplicate_source_taxon_id_cannot_be_inferred_as_two_species(tmp_path):
    counts, source, execution = fixture(tmp_path)
    content = counts.read_text(encoding="utf-8").replace(
        "14,Test species D,5", "12,Test species D,5")
    counts.write_text(content, encoding="utf-8")
    digest = hashlib.sha256(counts.read_bytes()).hexdigest()
    s = json.loads(source.read_text())
    e = json.loads(execution.read_text())
    s["counts_file_sha256"] = digest
    e["photo_counts_csv_sha256"] = digest
    source.write_text(json.dumps(s), encoding="utf-8")
    execution.write_text(json.dumps(e), encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        audit(counts, source, execution, species=4)
