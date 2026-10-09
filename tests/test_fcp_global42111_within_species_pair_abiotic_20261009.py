"""Source-exact pair matching and non-causal geographic/climate/soil model tests."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import analyze_fcp_global42111_within_species_pair_abiotic_20261009 as M


@pytest.fixture
def original_pair_fixture(monkeypatch):
    monkeypatch.setattr(M,"N_PAIRS",6)
    monkeypatch.setattr(M,"BOTH_CLASSIFIABLE",4)
    monkeypatch.setattr(M,"DISCORDANT",2)
    pair=[]
    env=[]
    for t in range(1,7):
        pair.append({
            "inat_taxon_id":t,"species":f"Genus{t%2} plant{t}",
            "photo_id_1":1000+t,"photo_id_2":2000+t,
            "observation_id_1":3000+t,"observation_id_2":4000+t,
            "observer_id_1":5000+t,"observer_id_2":6000+t,
            "cell_id_1":99,"cell_id_2":100,
            "pair_state":"unclassifiable" if t>4 else "discordant" if t%2==0 else "same",
            "both_endpoints_classifiable":t<=4,
        })
        for side in (1,2):
            row={
                "inat_taxon_id":t,"photo_id":side*1000+t,
                "observation_id":(side+2)*1000+t,
                "site_geo_status":M.SOURCE_GEO,"latitude":12.0,
                "longitude":10. if side==1 else 30.,
                "environment_climate_complete":True,
                "environment_soil_complete":True,
                "environment_all_complete":True,
                "wc_elevation_m":50. if side==1 else 75.,
            }
            for k in M.CLIMATE+M.SOIL:row[k]=t*2.+side
            env.append(row)
    return pd.DataFrame(pair),pd.DataFrame(env)


def test_immutable_sample_exact_identity_and_unclassifiable_preserved(original_pair_fixture):
    pairs,env=original_pair_fixture
    f,s=M.join_origins(pairs,env)
    assert s["n_original_species_unique_cross_cell_pairs"]==6
    assert s["n_original_pairs_both_four_state_classified"]==4
    assert s["n_original_pairs_discordant"]==2
    assert s["n_original_pairs_missing_classified_response"]==2
    assert len(f)==4
    assert s["n_classified_pairs_all_geo_climate_soil"]==4
    assert f["geodesic_distance_km"].gt(100).all()


def test_source_photo_id_and_observation_must_match(original_pair_fixture):
    pairs,env=original_pair_fixture
    env.loc[0,"observation_id"]=999999
    with pytest.raises(RuntimeError,match="endpoint"):
        M.join_origins(pairs,env)


def test_never_conflate_suppressed_soil_with_photo_colour_match(original_pair_fixture):
    pairs,env=original_pair_fixture
    env.loc[0,"soil_pH"]=np.nan
    env.loc[0,"environment_soil_complete"]=False
    env.loc[0,"environment_all_complete"]=False
    f,s=M.join_origins(pairs,env)
    assert len(f)==3
    assert s["n_original_pairs_both_four_state_classified"]==4


def test_cannot_swap_same_observer_or_same_cell(original_pair_fixture):
    pairs,env=original_pair_fixture
    pairs.loc[0,"observer_id_2"]=pairs.loc[0,"observer_id_1"]
    with pytest.raises(ValueError,match="Same cell/observer"):
        M.join_origins(pairs,env)


def test_antipodal_sites_have_undefined_midpoint():
    d,c=M.pair_distance_and_midpoint([0],[0],[0],[180])
    assert d[0]>20000
    assert c[0]==-1


def test_small_cohort_reports_hold_not_fake_effect(original_pair_fixture):
    p,e=original_pair_fixture
    d=M.analyze(p,e)
    assert d["status"]=="HOLD_SOURCE_PAIR_SOIL_COMPLETE_COVERAGE_OR_COLOUR_CLASS"
    assert d["original_pairs_unclassifiable_preserved_in_denominator"]
    assert d["models"]==[]


def test_out_of_genus_and_spatial_fold_same_fixed_pairs():
    rng=np.random.default_rng(102026)
    n=700
    ids=np.arange(n)
    y=(ids%3==0)
    d=pd.DataFrame({"pair_state":np.where(y,"discordant","same"),
                    "genus":[f"Genus{i%80}" for i in ids],
                    "pair_midpoint_cell_162":ids%35,
                    "log_geodesic_distance_km":rng.normal(6,1,n),
                    "delta_abs_latitude":rng.uniform(0,45,n),
                    "delta_elevation_m":rng.uniform(0,1000,n)})
    for feature in M.CLIM_FEATURES+M.SOIL_FEATURES:
        d[feature]=rng.normal(1,0.4,n)
    for by in ("genus","pair_midpoint_cell_162"):
        ans=M.oof_compare(d,by)
        assert ans["status"]=="FIXED_ORIGINAL_PAIRS_BLOCKED_OOF_DIAGNOSTIC"
        assert ans["n_pairs"]==n
        assert set(ans["models"])==set(M.FAMILIES)
        assert len(ans["group_bootstrap_95CI_logloss_gain_soil"])==2
        assert np.isfinite(ans["soil_increment_outofgroup_logloss"])
