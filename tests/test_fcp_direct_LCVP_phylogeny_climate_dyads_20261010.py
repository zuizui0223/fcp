"""Source-safe directly tipped LCVP photo dyad vs rainfall, no causal claim."""
from __future__ import annotations
from io import StringIO
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import assess_fcp_direct_LCVP_phylogeny_climate_dyads_20261010 as M


@pytest.fixture
def phylo_fixture():
    tree=Phylo.read(StringIO("(((G_a:1,G_b:1):1,G_c:2):1,((H_a:1,H_b:1):1,H_c:2):1);"),"newick")
    rows=[]
    for i,(genus,name,lat,lon,morph) in enumerate([
        ("G","G_a",0.,10.,"white"),
        ("G","G_b",0.,10.2,"white"),
        ("G","G_c",0.1,10.1,"red_pink"),
        ("H","H_a",4.,20.,"blue_purple"),
        ("H","H_b",4.,20.2,"yellow_orange"),
        ("H","H_c",4.1,20.1,"blue_purple"),
    ]):
        rows.append({"inat_taxon_id":i+1,"genus":genus,
           "genus_cell_id":genus+"|50","photo_cell_162":50,
           "original_direct_LCVP_tip":name,"morph":morph,
           "latitude":lat,"longitude":lon,"wc_elevation_m":100+i*3,
           "wc_bio1":12+i*.2,"wc_bio5":28+i*.1,
           "wc_bio12":500+i*20,"wc_bio15":40+i*.3})
    return pd.DataFrame(rows),tree


def test_exact_original_geo_phylo_pair_features_no_fake_species(phylo_fixture):
    data,tree=phylo_fixture
    pairs,cover=M.original_dyads(data,tree,50)
    assert len(pairs)==6
    assert cover["n_original_congeneric_photo_pair_dyads"]==6
    assert cover["n_original_species_in_local_direct_tip_comparisons"]==6
    assert cover["n_original_genera_in_local_tip_comparisons"]==2
    assert cover["n_photo_colour_mismatched_dyads"]==4
    assert pairs["log_LCVP_patristic"].nunique()==2
    assert pairs["log_km"].gt(0).all()
    assert cover["group_selection_did_not_use_photo_colour_or_rainfall"] is True


def test_local_groups_are_selected_blind_to_original_colour(phylo_fixture):
    source,tree=phylo_fixture
    original,stats=M.original_dyads(source,tree,50)
    changed=source.copy()
    changed["morph"]="white"
    alt,after=M.original_dyads(changed,tree,50)
    assert after["n_original_congeneric_photo_pair_dyads"]==stats["n_original_congeneric_photo_pair_dyads"]
    assert after["n_original_species_in_local_direct_tip_comparisons"]==stats["n_original_species_in_local_direct_tip_comparisons"]
    assert original[["taxon_1","taxon_2","log_km","log_LCVP_patristic"]].equals(
        alt[["taxon_1","taxon_2","log_km","log_LCVP_patristic"]])


def test_no_duplicate_original_species_or_invented_locality(phylo_fixture):
    x,tree=phylo_fixture
    altered=pd.concat([x,x.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicated"):
        M.original_dyads(altered,tree,50)
    x.loc[0,"latitude"]=np.nan
    with pytest.raises(ValueError,match="missing or invalid"):
        M.original_dyads(x,tree,50)


def synthetic_pairs():
    rng=np.random.default_rng(12026)
    rows=[]
    for g in range(40):
        for k in range(7):
            row={
                "genus":f"G{g}",
                "local_group":f"G{g}|original{g}",
                "photo_cell":g%30,
                "taxon_1":g*15+1+k,"taxon_2":g*15+9+k,
                "colour_photo_mismatch":int((g+k)%3!=0)
            }
            for feature in M.MODELS["GEO_THERMAL_PHYLO_RAIN"]:
                row[feature]=float(rng.uniform(0.001,6))
            rows.append(row)
    return pd.DataFrame(rows)


def test_genus_group_heldout_photo_pair_comparisons_same_dyads():
    d=synthetic_pairs()
    ans=M.predictive_check(d)
    assert ans["status"]=="EXPLORATORY_DIRECT_DATED_PHYLO_CLIMATE_PAIRWISE_PREDICTION"
    assert ans["n_original_photo_pair_dyads"]==280
    assert ans["n_actual_original_genera_held_out"]==40
    assert set(ans["models"])==set(M.MODELS)
    assert set(ans["incremental_prediction"])==set(M.COMPARE)
    assert {v["n_same_heldout_source_dyads"] for v in ans["models"].values()}=={280}
    assert all(len(v["original_genus_block_bootstrap_95CI"])==2
               for v in ans["incremental_prediction"].values())


def test_too_few_or_all_identical_photo_mismatch_is_hold():
    d=synthetic_pairs().iloc[:50]
    assert M.predictive_check(d)["status"]=="HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLO_DYAD_PREDICTION"
    z=synthetic_pairs()
    z["colour_photo_mismatch"]=0
    assert M.predictive_check(z)["status"]=="HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLO_DYAD_PREDICTION"


def test_original_full_plant_phylogeny_sufficiency_stays_hold(phylo_fixture):
    x,tree=phylo_fixture
    _,info=M.original_dyads(x,tree,100)
    assert info["n_original_species_in_local_direct_tip_comparisons"]==6
    assert M.MIN_SPECIES==100
    assert M.MIN_GENERA==20
