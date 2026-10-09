"""Colour-blind nearby source-species photo subset null: no outcome selection."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_microspatial_subset_scoring_20261009 as S
import audit_fcp_global42111_local_congener_moisture_null_20261009 as P
from audit_fcp_global42111_local_congeners_20261009 import fold_plan, CLASSES


@pytest.fixture
def original_local():
    rng=np.random.default_rng(20261009)
    z=[]
    for g in range(22):
        for j in range(15):
            lat=12+rng.uniform(-.01,.01)
            lon=-130+10*(g%9)+rng.uniform(-.01,.01)
            row={"inat_taxon_id":g*15+j+1,"genus":f"Genus{g}",
                 "genus_cell_id":f"Genus{g}|{g%9}","photo_cell_162":g%9,
                 "latitude":lat,"longitude":lon,
                 "morph":CLASSES[(j+g)%4],"abs_latitude":abs(lat),
                 "lon_sin":np.sin(np.deg2rad(lon)),"lon_cos":np.cos(np.deg2rad(lon)),
                 "wc_elevation_m":rng.uniform(20,300)}
            for k in ("wc_bio1","wc_bio5","wc_bio12","wc_bio15"):
                row[k]=rng.uniform(1,100)
            z.append(row)
    return pd.DataFrame(z)


def test_new_per_photo_gain_reproduces_exact_original_all_source_model(original_local):
    fold=fold_plan(original_local)
    v=S.original_fixed_fold_photo_gain(original_local,fold)
    prev=P.outofsource_photo_moisture_gain(original_local,fold)
    assert len(v)==len(original_local)
    assert float(v.mean())==pytest.approx(prev,abs=1e-12)


def test_eval_species_identity_chosen_without_original_flower_colour(original_local):
    a=S.evaluate_subset(original_local,50,nperm=1,seed=123)
    changed=original_local.copy()
    changed["morph"]="white"
    b=S.evaluate_subset(changed,50,nperm=1,seed=123)
    for field in ("n_original_evaluation_photo_ids_in_multispecies_microgroups",
                  "source_species_fraction_scored",
                  "n_original_multispecies_microgroups"):
        assert a[field]==b[field]
    assert a["test_species_selection_uses_only_original_photo_coordinates_and_species_ID"]
    assert a["n_fixed_null_model_refits"]==1
    assert .5<=a["permutation_p"]<=1


def test_micro_neighborhood_without_two_source_taxa_returns_hold(original_local):
    d=S.evaluate_subset(original_local,1e-5,nperm=3,seed=2)
    assert d["status"]=="HOLD_INSUFFICIENT_COLOUR_BLIND_MULTISPECIES_PHOTO_NEIGHBORHOODS"
    assert d["n_fixed_null_model_refits"]==0
    assert d["permutation_p"] is None


def test_source_multi_colour_microgroups_not_needed_to_define_eligible_photo_IDs(original_local):
    k=original_local.copy()
    k.loc[k.genus.eq("Genus0"),"morph"]="white"
    a=S.evaluate_subset(original_local,50,nperm=1,seed=1)
    b=S.evaluate_subset(k,50,nperm=1,seed=1)
    assert a["n_original_evaluation_photo_ids_in_multispecies_microgroups"]==330
    assert b["n_original_evaluation_photo_ids_in_multispecies_microgroups"]==330
    assert a["microgroup_four_colour_composition_preserved_per_null"]
