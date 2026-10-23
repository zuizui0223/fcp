"""Real-photograph polymorphism screening contracts, using synthetic images only.

The tests do NOT claim a photograph from iNaturalist proves genetic local FCP.
Actual study inclusion depends on sealed freshly remeasured ROI-v4 outputs.
"""
from __future__ import annotations

import csv
import json

import pytest

from fcp_pipeline.photo_polymorphism_screen import (
    FRAME_SHA, MORPHS, PALETTE, POSITIVE, PROTOCOL, SNAPSHOT,
    observer_split, screen, sha256,
)


def fake_data(tmp_path):
    frame = tmp_path / "frame.csv"
    with frame.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["inat_taxon_id", "species"])
        w.writerows([[11, "Spec alba"], [12, "Spec rubra"],
                     [13, "Spec lutea"], [14, "Spec lilac"], [15, "Spec absent"]])
    measured = tmp_path / "new_sealed_original_photos.csv"
    fields = [
        "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
        "morph", "measurement_status", "global_classifiable",
        *(f"flower_fraction_{x}" for x in PALETTE),
    ]
    rows = []
    identity = {11: "Spec alba", 12: "Spec rubra",
                13: "Spec lutea", 14: "Spec lilac", 15: "Spec absent"}

    def photo(taxon, index, morph, observer="observer", obs=None):
        is_good = morph in MORPHS
        row = {
            "inat_taxon_id": taxon, "species": identity[taxon],
            "photo_id": f"{taxon}_{index}",
            "observation_id": obs or f"{taxon}_obs_{index}",
            "observer_id": observer,
            "morph": morph if is_good else "mixed_uncertain",
            "measurement_status": ("classified_four_state_morph"
                                   if is_good else "not_evaluable_roi_or_flip_gate"),
            "global_classifiable": str(is_good).lower(),
        }
        main_pal = {
            "white": "white", "yellow_orange": "yellow",
            "red_pink": "pink", "blue_purple": "blue",
        }.get(morph, "white")
        row.update({f"flower_fraction_{x}": "1" if x == main_pal else "0"
                    for x in PALETTE})
        rows.append(row)

    def observers_in_split(taxon, split, n):
        found = []
        k = 0
        while len(found) < n:
            label = f"camera_{taxon}_{k}"
            if observer_split(taxon, label) == split:
                found.append(label)
            k += 1
        return found

    # A: independent observer-disjoint replication of white/red, 60 photos.
    for fold in range(2):
        ids = observers_in_split(11, fold, 3)
        for i in range(20):
            photo(11, fold * 30 + i, "white", ids[i % 3])
        for i in range(10):
            photo(11, fold * 30 + 20 + i, "red_pink", ids[i % 3])
    # B: two independent observers but both fall into the same split.
    same_split_observers = observers_in_split(12, 0, 2)
    for i in range(50):
        photo(12, i, "white" if i < 35 else "blue_purple",
              same_split_observers[i % 2])
    # C: 60 photos all white, never claim biological monomorphism.
    for i in range(60):
        photo(13, i, "white", f"single_{i % 6}")
    # D: <40 classifiable including unclassifiable original flowers.
    for i in range(20):
        photo(14, i, "white" if i % 2 else "red_pink", f"d_{i % 4}")
    for i in range(20, 40):
        photo(14, i, "mixed_uncertain", f"d_{i % 4}")

    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, fields)
        writer.writeheader()
        writer.writerows(rows)

    counts = tmp_path / "snapshot_photo_counts.csv"
    with counts.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"])
        writer.writerows([[11, identity[11], 200, SNAPSHOT],
                         [12, identity[12], 250, SNAPSHOT],
                         [13, identity[13], 100, SNAPSHOT],
                         [14, identity[14], 100, SNAPSHOT],
                         [15, identity[15], 0, SNAPSHOT]])
    receipt = tmp_path / "sealed_photo_screen_receipt.json"
    receipt.write_text(json.dumps({
        "schema": "fcp_sealed_photo_screen_terminal_v1",
        "status": "SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL",
        "protocol": PROTOCOL,
        "source_date": SNAPSHOT,
        "source_frame_sha256": FRAME_SHA,
        "source_species": 5,
        "new_photo_table_sha256": sha256(measured),
        "original_ROI_v4_method_attested": True,
        "photo_ID_selection_pre_colour": True,
        "source_photo_IDs_sealed": True,
        "historical_measured_flower_colours_reused": False,
        "n_terminal_measurements": len(rows),
    }))
    return frame, counts, measured, receipt


def read_table(path):
    with path.open(newline="") as src:
        return list(csv.DictReader(src))


def test_real_cohort_inclusion_is_photo_based_not_literature(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    out = tmp_path / "output"
    summary = screen(measured, frame, counts, receipt, out, _synthetic=True)
    assert summary["source_species"] == 5
    assert summary["photo_supported_polymorphic_species"] == 1
    assert summary["complete_42111_photo_census"] is False
    assert summary["previous_34_literature_C_S_labels_used"] is False
    assert summary["historical_149900_measured_colour_rows_reused"] is False
    led = {int(r["inat_taxon_id"]): r for r in read_table(
        out / "all_42111_photo_polymorphism_screen.csv")}
    assert led[11]["photographic_colour_status"] == POSITIVE
    assert led[12]["photographic_colour_status"] == "photo_variation_candidate_holdout_unresolved"
    assert led[13]["photographic_colour_status"] == "single_photo_state_detected_not_monomorphic"
    assert led[14]["photographic_colour_status"] == "insufficient_classifiable_not_monomorphic"
    assert led[15]["photographic_colour_status"] == "not_screened_no_terminal_photo"
    assert all(x["genetic_polymorphism_proven"] == "False" for x in led.values())
    assert [r["species"] for r in read_table(
        out / "photo_supported_polymorphic_species.csv")] == ["Spec alba"]
    assert len(read_table(
        out / "not_proven_polymorphic_do_not_call_monomorphic.csv")) == 4


def test_observer_disjoint_split_prevents_duplicate_photo_pseudoreplication(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    # Convert B to one observer: many photos of both colours no longer count
    # as two independent natural individuals or independently replicated morphs.
    lines = read_table(measured)
    b = next(r["observer_id"] for r in lines if r["inat_taxon_id"] == "12")
    for row in lines:
        if row["inat_taxon_id"] == "12":
            row["observer_id"] = b
    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(lines[0]))
        writer.writeheader()
        writer.writerows(lines)
    doc = json.loads(receipt.read_text())
    doc["new_photo_table_sha256"] = sha256(measured)
    receipt.write_text(json.dumps(doc))
    out = tmp_path / "only_one_observer"
    screen(measured, frame, counts, receipt, out, _synthetic=True)
    records = {int(x["inat_taxon_id"]): x for x in read_table(
        out / "all_42111_photo_polymorphism_screen.csv")}
    assert records[12]["photographic_colour_status"] == "photo_colour_signal_weak_or_unreplicated"


def test_require_original_method_and_exact_measurement_digest(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    doc = json.loads(receipt.read_text())
    doc["historical_measured_flower_colours_reused"] = True
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="Incomplete or changed ORIGINAL"):
        screen(measured, frame, counts, receipt, tmp_path / "bad", _synthetic=True)


def test_broken_original_morph_status_and_palette_rejected(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    rows = read_table(measured)
    rows[0]["global_classifiable"] = "false"
    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    doc = json.loads(receipt.read_text())
    doc["new_photo_table_sha256"] = sha256(measured)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="Photo morph/classifiable/ROI status mismatch"):
        screen(measured, frame, counts, receipt, tmp_path / "bad", _synthetic=True)


def test_duplicate_photo_id_fails_closed(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    rows = read_table(measured)
    rows[1]["photo_id"] = rows[0]["photo_id"]
    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    doc = json.loads(receipt.read_text())
    doc["new_photo_table_sha256"] = sha256(measured)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="Duplicate/blank"):
        screen(measured, frame, counts, receipt, tmp_path / "bad", _synthetic=True)


def test_no_valid_sealed_receipt_no_picture_polymorphism_claim(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    receipt.write_text(json.dumps({"schema": "fcp_literature_only_26"}))
    with pytest.raises(ValueError, match="No genuine sealed"):
        screen(measured, frame, counts, receipt, tmp_path / "bad", _synthetic=True)


def test_reassembly_requires_exact_species_name_and_sorted_taxa(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    rows = read_table(measured)
    rows[0]["species"] = "Changed taxon"
    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    doc = json.loads(receipt.read_text())
    doc["new_photo_table_sha256"] = sha256(measured)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="Fresh photo taxon name drift"):
        screen(measured, frame, counts, receipt, tmp_path / "bad", _synthetic=True)


def test_no_40_classifiable_photos_cannot_be_proven_monomorphic(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    out = tmp_path / "output"
    screen(measured, frame, counts, receipt, out, _synthetic=True)
    status = {int(x["inat_taxon_id"]): x["photographic_colour_status"]
              for x in read_table(out / "all_42111_photo_polymorphism_screen.csv")}
    assert status[14].startswith("insufficient_")
    assert "monomorphic" not in status[15]


def test_double_colour_inside_one_original_flower_roi_is_qc_hold_not_fcp(tmp_path):
    frame, counts, measured, receipt = fake_data(tmp_path)
    rows = read_table(measured)
    for row in rows:
        if row["inat_taxon_id"] != "11":
            continue
        if row["morph"] == "red_pink":
            row["morph"] = "yellow_orange"
        # The ONE original photo ROI contains BOTH white and yellow:
        # high repeated label stability must not imply independent morphs.
        for colour in PALETTE:
            row[f"flower_fraction_{colour}"] = (
                ".52" if colour == "white" else
                ".48" if colour == "yellow" else "0")
    with measured.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    doc = json.loads(receipt.read_text(encoding="utf-8"))
    doc["new_photo_table_sha256"] = sha256(measured)
    receipt.write_text(json.dumps(doc), encoding="utf-8")
    out = tmp_path / "postpilot_structural_screen"
    result = screen(measured, frame, counts, receipt, out, _synthetic=True)
    species = {r["inat_taxon_id"]:r for r in read_table(
        out / "all_42111_photo_polymorphism_screen.csv")}
    source = species["11"]
    assert source["structural_flower_architecture_manual_review_required"] == "True"
    assert float(source["within_ROI_white_yellow_both_ge10_fraction"]) == 1.0
    assert source["photographic_colour_status"] == (
        "repeated_coarse_photo_labels_STRUCTURAL_QC_HOLD_not_FCP")
    assert result["n_repeated_photo_states_structural_ROI_review_holds"] == 1
    assert result["photo_supported_polymorphic_species"] == 0
    assert len(read_table(out / "photo_supported_polymorphic_species.csv")) == 0
    assert len(read_table(
        out / "not_proven_polymorphic_do_not_call_monomorphic.csv")) == 5
