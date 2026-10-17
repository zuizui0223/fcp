"""Synthetic audit for new-photo concatenation only, not empirical FCP inference."""
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.analysis import append_wave4_to_original_149900_20261017 as a


def build(tmp):
    targets = {"discovery": 1, "validation": 1, "prospective": 1, "wave4": 300}
    paths = {}
    contract = {}
    for label, n in targets.items():
        rows = []
        for i in range(n):
            taxon = len(rows) + 100000 + {"discovery": 0, "validation": 1000, "prospective": 2000, "wave4": 3000}[label] + i
            for k in range(2):
                rows.append({
                    "photo_id": ({"discovery": 10000, "validation": 20000, "prospective": 30000, "wave4": 50000}[label] + 2*i + k),
                    "inat_taxon_id": taxon, "species": f"{label} species_{i}",
                    "morph": "white" if k == 0 else "blue_purple",
                    "global_classifiable": True, "observer_id": f"observer_{i}_{k}",
                    "latitude": 0.0, "longitude": i * .001,
                    "flower_fraction_white": 1.0 if k == 0 else 0.0,
                })
        path = tmp / (label + ".csv")
        pd.DataFrame(rows).to_csv(path, index=False, lineterminator="\n")
        paths[label] = path
        contract[label] = {"sha256": a.sha256_file(path), "rows": 2*n, "species": n}
    receipt = {
        "schema": "fcp_same_149900_method_wave4_real_photo_terminals_v1",
        "original_results_changed": False,
        "measured_sha256": a.sha256_file(paths["wave4"]),
        "terminal_partitions": 256,
        "measurement_original_fixed_colour_worker": True,
        "wave4_100photo_species": 300,
        "wave4_new_photo_ids": 600,
        "wave4_classifiable_photo_rows": 600,
    }
    rp = tmp / "wave4_terminal.json"
    rp.write_text(json.dumps(receipt))
    return {k: v for k, v in paths.items() if k != "wave4"}, paths["wave4"], rp, {k: v for k, v in contract.items() if k != "wave4"}


def test_append_preserves_original_labels_and_adds_600_new_rows(tmp_path):
    old, new, receipt, old_contract = build(tmp_path)
    out, r = a.combine_old_and_wave4(old, new, receipt, old_frame_contract=old_contract, target=2)
    assert len(out) == 606
    assert r["original_source_rows"] == 6
    assert r["new_original_photo_rows"] == 600
    assert r["all_named_species"] == 303
    assert out["photo_id"].nunique() == len(out)
    assert r["append_only"] is True and r["four_state_method_relabelled"] is False
    assert r["existing_50km_null_reweighted"] is False
    assert out.loc[out.fcp_source_cohort == "discovery", "morph"].tolist() == ["white", "blue_purple"]


def test_changed_historical_source_fails_closed(tmp_path):
    old, new, receipt, contract = build(tmp_path)
    old["discovery"].write_text(old["discovery"].read_text().replace("white", "red_pink"))
    with pytest.raises(ValueError, match="immutable original source SHA256"):
        a.combine_old_and_wave4(old, new, receipt, old_frame_contract=contract, target=2)


def test_new_wave_original_photo_id_intersection_rejected(tmp_path):
    old, new, receipt, contract = build(tmp_path)
    df = pd.read_csv(new)
    df.loc[0, "photo_id"] = 10000  # a genuine old photo identity
    df.to_csv(new, index=False)
    rp = json.loads(receipt.read_text())
    rp["measured_sha256"] = a.sha256_file(new)
    receipt.write_text(json.dumps(rp))
    with pytest.raises(ValueError, match="identity reused"):
        a.combine_old_and_wave4(old, new, receipt, old_frame_contract=contract, target=2)


def test_missing_original_wave4_terminal_rejected(tmp_path):
    old, new, receipt, contract = build(tmp_path)
    record = json.loads(receipt.read_text())
    record["measurement_original_fixed_colour_worker"] = False
    receipt.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="completed unchanged-method"):
        a.combine_old_and_wave4(old, new, receipt, old_frame_contract=contract, target=2)


def test_species_disjoint_does_not_allow_temporary_taxonomy_identity_reuse(tmp_path):
    old, new, receipt, contract = build(tmp_path)
    df = pd.read_csv(new)
    df.loc[0:1, "species"] = "discovery species_0"
    df.to_csv(new, index=False)
    rp = json.loads(receipt.read_text())
    rp["measured_sha256"] = a.sha256_file(new)
    receipt.write_text(json.dumps(rp))
    with pytest.raises(ValueError, match="identity reused"):
        a.combine_old_and_wave4(old, new, receipt, old_frame_contract=contract, target=2)


def test_fewer_than_300_new_species_is_underidentified_not_zero_colour(tmp_path):
    old, new, receipt, contract = build(tmp_path)
    record = json.loads(receipt.read_text())
    record["wave4_100photo_species"] = 299
    record["wave4_new_photo_ids"] = 598
    receipt.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="sampling denominator"):
        a.combine_old_and_wave4(old, new, receipt, old_frame_contract=contract, target=2)
