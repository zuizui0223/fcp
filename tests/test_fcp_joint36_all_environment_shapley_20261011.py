"""Original FCP source and all 2**N environmental Shapley integrity tests."""
from __future__ import annotations
import math,os,sys,json
from pathlib import Path
import numpy as np
import pytest
SCRIPTS=Path(__file__).resolve().parents[1]/"scripts"/"analysis"
sys.path.insert(0,str(SCRIPTS))
from fcp_audit_36_preflight_20261011 import BLOCKS,FEATURES,CLASSES
from fcp_environment_joint_shapley_20261011 import SCHEMA
ROOT=Path(os.getenv("FCP_RESULT_ROOT","/tmp/fcp_joint_environment_results_20261011"))

def test_every_measured_environment_feature_in_exactly_one_group():
    assert len(BLOCKS)==7
    assert sorted(len(x) for x in BLOCKS.values())==[1,2,2,2,8,10,11]
    assert len(FEATURES)==36 and len(set(FEATURES))==36
    assert set(FEATURES)==set(z for seq in BLOCKS.values() for z in seq)
    assert tuple(CLASSES)==("white","yellow_orange","red_pink","blue_purple")

def test_shapley_weights_are_exactly_telescope_for_any_coalition_game():
    random=np.random.default_rng(31)
    for n in (6,7):
        values=random.normal(size=2**n)
        allocation=[]
        for j in range(n):
            allocation.append(sum(
                (values[S|1<<j]-values[S])/(n*math.comb(n-1,S.bit_count()))
                for S in range(2**n) if not S&(1<<j)))
        assert np.isclose(sum(allocation),values[-1]-values[0])

def test_environment_model_is_prediction_not_causal_fitness():
    import fcp_environment_joint_shapley_20261011 as model
    assert "Shapley" in (model.__doc__ or "")
    assert SCHEMA=="fcp_joint_abiotic_shapley_all_original_photo_species_geocell_v1"

@pytest.mark.parametrize("mode,region,n,k", [
    ("soil","world",13961,1167),
    ("soil","east",3307,346),
    ("soil","elsewhere",10307,914),
    ("climate","world",17183,1347),
    ("climate","east",4555,455),
    ("climate","elsewhere",12653,1071),
])
def test_exact_original_source_numerical_replay_if_case_has_run(mode,region,n,k):
    path=ROOT/f"{mode}_{region}"/"result.json"
    if not path.exists():pytest.skip("Use only on source replay results present")
    r=json.loads(path.read_text())
    assert r["source"]["photo_union"]==100543
    assert r["source"]["species_opportunity"]==42111
    assert r["source"]["photo_original_classified"]==45953
    assert r["source"]["original_unclassifiable_kept_missing"]==54590
    assert r["source"]["original_full_climate_soil_case_coverage"]["climate_complete"]==45836
    assert r["source"]["original_full_climate_soil_case_coverage"]["soil_complete"]==31925
    key=f"{mode}_{region}"
    a=r["cohorts"][key]
    assert a["schema"]==SCHEMA
    assert a["n_oof_species"]==n
    assert a["n_oof_genera"]==k
    blocks=a["grouped_shapley_of_joint_out_of_geo_cell_predictive_information"]
    assert len(blocks)==(7 if mode=="soil" else 6)
    assert a["n_subsets_scored"]==(128 if mode=="soil" else 64)
    assert a["shapley_telescopes_to_full_environment_gain"] is True
    tot=a["full_environment_increment_beyond_geo"]
    assert np.isclose(sum(z["genus_equal_shapley_Brier"] for z in blocks.values()),
                      tot["genus_equal_delta_brier"],atol=1e-10)
    assert np.isclose(sum(z["species_equal_shapley_Brier"] for z in blocks.values()),
                      tot["species_equal_delta_brier"],atol=1e-10)
    assert len(a["folds"])==5
    assert a["n_oof_species"]==sum(f["n_classifiable_eligible_evaluated"] for f in a["folds"])

@pytest.mark.parametrize("mode,n", [("soil",3399),("climate",4512)])
def test_strict_quality_controls_do_not_silently_change_photo_standard(mode,n):
    f=ROOT/f"{mode}_world_quality"/"result.json"
    if not f.exists():pytest.skip("Quality source replay pending")
    r=json.loads(f.read_text())
    assert r["quality_rule"]["median_flower_effective_pixels_threshold"]==23826.
    assert np.isclose(r["quality_rule"]["median_nuisance_pixel_fraction_max"],0.1375717742195135)
    assert r["cohorts"][f"{mode}_world_quality"]["n_oof_species"]==n
