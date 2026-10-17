"""Synthetic safety tests; artificial flowers are never presented as results."""
import json
import numpy as np
import pandas as pd
import pytest

from scripts.analysis.replay_fcp_149900_plus_wave4_same_spatial_20261017 import (
    GROUPS, run_four_cohorts, require_reproduced_original_three,
    verify_existing_wave4_inference, import_exact_measured_rows,
)


def four_fake_cohorts():
    out=[]
    for ci, cohort in enumerate(GROUPS):
        for k in range(40):
            out.append({
                "fcp_source_cohort": cohort,
                "inat_taxon_id": ci + 500000,
                "species": f"Syntheticus {ci}",
                "photo_id": ci*1000+k+1,
                "morph": "white" if k < 20 else "blue_purple",
                "global_classifiable": True,
                "latitude": 10.0,
                "longitude": 0.0 if k < 20 else 10.0,
                "observer_id": f"observer{ci}_{k}",
            })
    return pd.DataFrame(out)


def test_same_original_50km_core_keeps_four_cohorts_and_original_vertex_null():
    stat, species = run_four_cohorts(four_fake_cohorts())
    assert len(species)==4
    assert stat["pooled"]["50"]["n_species"]==4
    assert stat["pooled_n_spatial_eligible_species"]==4
    observed = [stat["cohorts"][c]["radii"]["50"]["mean_depletion"] for c in GROUPS]
    assert np.isclose(stat["pooled"]["50"]["mean_depletion"],np.mean(observed))
    assert stat["pooled"]["50"]["mean_depletion_p_upper"] <= .05


def test_original_real_frozen_baseline_gate_does_not_accept_synthetic():
    stat, _ = run_four_cohorts(four_fake_cohorts())
    with pytest.raises(ValueError,match="FROZEN_ORIGINAL_50KM_REPLAY_MISMATCH"):
        require_reproduced_original_three(stat)


def test_wave4_link_checks_observed_effect_and_preserves_sampling_order_disclosure():
    stat,_ = run_four_cohorts(four_fake_cohorts())
    number=stat["cohorts"]["wave4"]["radii"]["50"]
    receipt={"input_sha256":{"wave4":"dummy"}}
    w={
        "schema":"fcp_same_149900_method_wave4_statistical_replication_v1",
        "original_frozen_results_unchanged":True,
        "wave4_measured_table_sha256":"dummy",
        "spatial_original_method":{"50":{
            "n_species":number["n_species"],
            "mean_depletion":number["mean_depletion"],
            "mean_depletion_p_upper":0.735,
        }},
        "h1_original_method":{"synthetic":True},
        "h2_original_W_method":{"synthetic":True},
    }
    res=verify_existing_wave4_inference(stat,receipt,w)
    assert res["new_wave4_original_50km_observed_effect_verified"]
    assert not res["wave4_50km_null_p_matches"]
    w["spatial_original_method"]["50"]["mean_depletion"]+=.01
    with pytest.raises(ValueError, match="changed during"):
        verify_existing_wave4_inference(stat,receipt,w)


def test_expanded_file_sha_gate_fail_closed(tmp_path):
    data=tmp_path/"fake.csv"
    four_fake_cohorts().to_csv(data,index=False)
    with pytest.raises(ValueError,match="checksum"):
        import_exact_measured_rows(data,{
            "schema":"fcp_original_149900_plus_same_method_wave4_merged_v1",
            "original_study_unchanged":True, "append_only":True,
            "four_state_method_relabelled":False,
            "existing_50km_null_reweighted":False,
            "merged_csv_sha256":"0"*64,
            "input_sha256":{},
        })
