"""Regression tests for the NEW literature-documented local-FCP-only study.

All evidence comes from the pre-existing, byte-frozen v2.2 positive-source audit;
no 42,111 source-photo opportunity taxon is assumed to have polymorphism.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.analysis.select_source_documented_fcp_cohort import (
    EXCLUDED_FIELDS, SELECTED_FIELDS, SOURCE, SOURCE_MANIFEST, STATUS,
    build, main, sha256, write_csv,
)

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_fixed_real_evidence_source_selects_only_natural_population_FCP():
    selected, excluded, manifest = build(ROOT / SOURCE, ROOT / SOURCE_MANIFEST)
    assert manifest["status"] == STATUS
    assert len(selected) == 26 and len(excluded) == 8
    assert manifest["local_coexistence_only"] == 11
    assert manifest["coexistence_and_segregation"] == 15
    assert all(row["C_local_coexistence_documented"] == "1" for row in selected)
    assert all(row["C_source_ids"] for row in selected)
    assert all(int(row["n_C_positive_sources"]) > 0 for row in selected)
    assert all(row["C_local_coexistence_documented"] == "0" for row in excluded)
    assert all(row["S_spatial_segregation_documented"] == "1" for row in excluded)
    assert all("human_duplicate_review" in row["evidence_review_status"] for row in selected)
    assert manifest["source_supported_not_independently_reviewer_verified"] is True
    assert manifest["all_42111_species_are_not_polymorphic_cases"] is True


def test_committed_26_taxa_and_8_quarantine_recreate_byte_for_byte(tmp_path):
    selection, exclusion, manifest = build(ROOT / SOURCE, ROOT / SOURCE_MANIFEST)
    keep = tmp_path / "strict.csv"
    drop = tmp_path / "excluded.csv"
    write_csv(keep, SELECTED_FIELDS, selection)
    write_csv(drop, EXCLUDED_FIELDS, exclusion)
    assert keep.read_bytes() == (ROOT / "data/frozen/strict_local_fcp_26_source_documented.csv").read_bytes()
    assert drop.read_bytes() == (
        ROOT / "data/frozen/geographic_only_8_excluded_from_strict_fcp.csv").read_bytes()
    assert not set(row["canonical_name"] for row in selection) & set(
        row["canonical_name"] for row in exclusion)
    assert "Silene littorea" in {row["canonical_name"] for row in selection}
    assert "Linanthus parryae" in {row["canonical_name"] for row in exclusion}


def test_no_image_census_or_biological_prevalence_derived_from_inclusion():
    _, _, result = build(ROOT / SOURCE, ROOT / SOURCE_MANIFEST)
    assert result["not_source_photo_validation"] is True
    assert result["previous_93209351_photo_ID_count_excluded_from_strict_FCP_measurement"] is True
    assert "cannot estimate prevalence" in result["warning"]


def test_rejects_changed_frozen_evidence_even_if_taxon_names_unchanged(tmp_path):
    source = tmp_path / "source.csv"
    source.write_bytes((ROOT / SOURCE).read_bytes().replace(b"Silene littorea", b"Silene armeria"))
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        build(source, ROOT / SOURCE_MANIFEST)


def test_rejects_changed_evidence_manifest(tmp_path):
    source_manifest = tmp_path / "manifest.json"
    original = json.loads((ROOT / SOURCE_MANIFEST).read_text(encoding="utf-8"))
    original["C_positive"] = 27
    source_manifest.write_text(json.dumps(original), encoding="utf-8")
    with pytest.raises(ValueError, match="evidence freeze changed"):
        build(ROOT / SOURCE, source_manifest)


def test_script_reproducibility_and_no_overwrite(tmp_path):
    output = tmp_path / "output"
    assert main(["--source", str(ROOT / SOURCE),
                 "--source-manifest", str(ROOT / SOURCE_MANIFEST),
                 "--output-dir", str(output)]) == 0
    result = json.loads((output / "cohort_receipt.json").read_text(encoding="utf-8"))
    assert result["selected_species"] == 26
    assert result["excluded_S_only"] == 8
    assert result["selected_csv_sha256"] == sha256(output / "strict_local_fcp_26_source_documented.csv")
    assert result["geographic_only_exclusion_csv_sha256"] == sha256(
        output / "geographic_only_8_excluded_from_strict_fcp.csv")
    with pytest.raises(ValueError, match="Never overwrite"):
        main(["--source", str(ROOT / SOURCE),
              "--source-manifest", str(ROOT / SOURCE_MANIFEST),
              "--output-dir", str(output)])
