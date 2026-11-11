"""Protect nonphylogenetic genus-equal photo-white global sensitivity."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location("genus_fcp",SRC/"audit_fcp_global_white_balance_genus_20261008.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def sample_species():
    return pd.DataFrame([
        {"inat_taxon_id":i,"genus":"Largus","region":"low_0_30","n_classified":4,
         "fraction_white_when_classified":1.0,"n_photo_white":4,
         "white_and_chromatic_images_in_same_region":False} for i in range(1,21)
    ]+[{"inat_taxon_id":21,"genus":"Rarus","region":"low_0_30",
         "n_classified":4,"fraction_white_when_classified":0.0,"n_photo_white":0,
         "white_and_chromatic_images_in_same_region":False}])


def test_genus_equal_not_dominated_by_twenty_congeneric_white_species():
    p=sample_species()
    r,d=m.within_region_genus_stat(p,"low_0_30",43)
    assert r["n_species_classifiable"]==21
    assert r["n_genera_classifiable"]==2
    assert np.isclose(r["mean_species_equal_photo_white_fraction"],20/21)
    assert r["mean_genus_equal_photo_white_fraction"]==0.5
    assert r["max_species_per_genus"]==20
    assert len(r["genus_cluster_bootstrap_95CI"])==2
    assert len(d)==2


def test_within_species_paired_genus_effect_uses_independent_genera():
    p=pd.DataFrame([
        {"inat_taxon_id":1,"region":"low_0_30","n_classified":4,
         "fraction_white_when_classified":1,"white_and_chromatic_images_in_same_region":False,"genus":"A"},
        {"inat_taxon_id":1,"region":"middle_30_60","n_classified":4,
         "fraction_white_when_classified":0,"white_and_chromatic_images_in_same_region":False,"genus":"A"},
        {"inat_taxon_id":2,"region":"low_0_30","n_classified":4,
         "fraction_white_when_classified":0,"white_and_chromatic_images_in_same_region":False,"genus":"B"},
        {"inat_taxon_id":2,"region":"middle_30_60","n_classified":4,
         "fraction_white_when_classified":1,"white_and_chromatic_images_in_same_region":False,"genus":"B"},
    ])
    result=m.genus_paired_region_delta(p,"low_0_30","middle_30_60",3)
    assert result["n_matched_species"]==2
    assert result["n_matched_genera"]==2
    assert np.isclose(result["genus_equal_paired_photo_white_delta"],0)
    assert result["bootstrap_CI_95_genus_cluster"][0]<0
    assert result["bootstrap_CI_95_genus_cluster"][1]>0


def test_does_not_create_genus_from_favorable_colour():
    d=pd.DataFrame({"inat_taxon_id":range(1,42112),
                    "species":[f"Genus{i%10} herb{i}" for i in range(1,42112)]})
    tax=m.source_genus(d)
    assert len(tax)==42111
    assert tax.inat_taxon_id.nunique()==42111
    assert tax.genus.nunique()==10
    assert "white" not in tax.columns


def test_no_region_support_is_an_explicit_hold():
    x=m.within_region_genus_stat(sample_species(),"high_60_90",77)[0]
    assert x["status"]=="HOLD_NO_CLASSIFIABLE_SPECIES"


def test_missing_genus_or_duplicate_taxon_id_is_fatal():
    d=pd.DataFrame({"inat_taxon_id":range(1,42112),
                    "species":[f"Genus{i} herb" for i in range(1,42112)]})
    d.loc[1,"inat_taxon_id"]=d.loc[0,"inat_taxon_id"]
    with pytest.raises(RuntimeError,match="registry"):
        m.source_genus(d)


def test_original_source_SHA256_and_no_mechanism_fitness():
    s=(SRC/"audit_fcp_global_white_balance_genus_20261008.py").read_text()
    assert m.SOURCE_COMMIT=="2b5390ea74f8d196012499d2d35dd477f9795938"
    assert len(m.BREADTH_SHA)==64 and len(m.CELL_SHA)==64
    assert "image_gen" not in s
    assert m.N_BOOT==1999
