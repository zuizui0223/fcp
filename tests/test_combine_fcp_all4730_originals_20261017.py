"""Synthetic protocol tests. Historical biological results are verified in real-data Actions."""
import json

import pandas as pd
import pytest

from scripts.analysis import combine_fcp_all4730_originals_20261017 as m


def fake_photo_table():
    frames=[]
    for i,cohort in enumerate(m.COHORTS):
        for k in range(40):
            frames.append({
                "photo_id":i*10000+k,
                "inat_taxon_id":i+500000,
                "species":f"Testus example{i}",
                "fcp_source_cohort":cohort,
                "global_classifiable":True,
                "morph":"white" if k<20 else "blue_purple",
                "latitude":35.0,
                "longitude":0.0 if k<20 else 10.0,
                "observer_id":f"test-{k}",
            })
    return pd.DataFrame(frames)


def test_exact_original_50km_replay_same_statistics_across_all_five_cohorts(monkeypatch):
    monkeypatch.setattr(m,"ORIGINAL_EXPECTED",{})
    species,rec=m.replay_same_old_50km(fake_photo_table())
    assert len(species)==5
    assert rec["all_evaluable_species"]==5
    assert rec["pooled_original_50km"]["n_species"]==5
    assert rec["old_three_original_149900_effects_exactly_reproduced"]
    assert rec["original_199_vertex_null_changed"] is False
    assert set(rec["all_sources"])==set(m.COHORTS)


def test_cannot_silently_change_original_discovery_50km_effect():
    with pytest.raises(RuntimeError,match="frozen statistic changed"):
        m.replay_same_old_50km(fake_photo_table())


def test_incomplete_wave4_actual_roi_receipt_cannot_be_appended(tmp_path):
    data=tmp_path/"some.csv"
    data.write_text("photo_id\n1\n")
    receipt=tmp_path/"terminal.json"
    receipt.write_text(json.dumps({
        "schema":"fcp_same_149900_method_wave4_real_photo_terminals_v1",
        "original_results_changed":False,
        "measurement_original_fixed_colour_worker":True,
        "terminal_partitions":255,
        "measured_sha256":m.sha(data),
        "wave4_100photo_species":498,
        "wave4_new_photo_ids":49800,
    }))
    with pytest.raises(RuntimeError,match="256-partition"):
        m.authenticate_terminal(receipt,data,wave4=True)


def test_incomplete_new_u100_roi_receipt_cannot_be_appended(tmp_path):
    data=tmp_path/"some.csv"
    data.write_text("photo_id\n1\n")
    receipt=tmp_path/"terminal.json"
    receipt.write_text(json.dumps({
        "schema":"fcp_all4730_new_original_method_pixel_terminal_v1",
        "original_roi_v4_and_efficientsam_unchanged":True,
        "original_frozen_149900_results_changed":False,
        "measurement_terminal_partitions":255,
        "measured_photo_file_sha256":m.sha(data),
        "new_complete_100_photo_species":1,
        "measurement_photo_rows":100,
    }))
    with pytest.raises(RuntimeError,match="SHA gate"):
        m.authenticate_terminal(receipt,data,wave4=False)


def test_no_new_species_imputed_for_missing_original_100_photo_gate():
    assert len(m.COHORTS)==5
    assert m.MAX_TAXA==4730
    assert 1499+498+2733==4730
