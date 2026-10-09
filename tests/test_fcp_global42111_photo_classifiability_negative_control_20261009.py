"""Synthetic classifiability negative control, retaining unclassified photos."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_photo_classifiability_negative_control_20261009 as M


@pytest.fixture
def original_photo_species():
    n=1200
    rng=np.random.default_rng(22026)
    i=np.arange(n)
    d=pd.DataFrame({
        "inat_taxon_id":i+1,
        "species":[f"Genus{j%100} plant{j}" for j in i],
        "morph":np.where(i%3==0,"white","UNKNOWN"),
        "measurement_status":np.where(i%3==0,"classified_four_state_morph","not_classified"),
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        "latitude":(i%18-9)*5.5,
        "longitude":(i%31-15)*11.,
        "environment_climate_complete":True,
        "environment_soil_complete":i%6!=0,
        "wc_elevation_m":rng.uniform(5,2000,n),
    })
    for col in M.TEMPERATURE+M.MOISTURE+M.SOIL:
        d[col]=rng.normal(10,2,n)
    d.loc[i%6==0,list(M.SOIL)]=np.nan
    return d


def test_all_source_photos_including_unclassifiable_remain_in_denominator(original_photo_species):
    pop,stat=M.source_sets(original_photo_species,strict=False)
    assert stat["original_source_taxa"]==1200
    assert stat["n_classifiable_original_photos"]==400
    assert stat["n_unclassifiable_original_photos"]==800
    assert len(pop["ALL_ORIGINAL_CLIMATE_PHOTO_OPPORTUNITY"])==1200
    assert len(pop["ALL_ORIGINAL_SOIL_PHOTO_OPPORTUNITY"])==1000


def test_training_only_genus_average_ignores_test_photo_classification():
    tr=pd.DataFrame({"inat_taxon_id":[1,2,3],
                     "genus":["A","A","A"],
                     "classifiable":[1,0,1],
                     "wc_bio12":[100.,150.,190.]})
    te=pd.DataFrame({"inat_taxon_id":[5],"genus":["A"],
                     "classifiable":[0],"wc_bio12":[125.]})
    base=M.predict_train_only_genus(tr,te,(),np.array([0]))
    moist=M.predict_train_only_genus(tr,te,("wc_bio12",),np.array([0]))
    assert base[0]==pytest.approx(2/3)
    assert 0<=moist[0]<=1
    te["classifiable"]=1
    new=M.predict_train_only_genus(tr,te,("wc_bio12",),np.array([0]))
    assert np.array_equal(moist,new)


def test_blocked_negative_control_identical_original_photo_denominators(original_photo_species):
    pop,_=M.source_sets(original_photo_species,strict=False)
    a=M.blocked_classifiability(pop["ALL_ORIGINAL_CLIMATE_PHOTO_OPPORTUNITY"],with_soil=False)
    b=M.blocked_classifiability(pop["ALL_ORIGINAL_SOIL_PHOTO_OPPORTUNITY"],with_soil=True)
    for z in (a,b):
        assert z["genus_means_estimated_on_training_photograph_opportunities_only"]
        assert z["same_original_source_photos_in_every_compared_model"]
        assert z["n_heldout_original_photos_genus_estimable"]>=100
        assert len(z["increments"]["unique_moisture_beyond_genus_geography_temperature"]["original_geographic_cell_bootstrap_95CI"])==2
        assert all(q["n_heldout_original_photos"]==z["n_heldout_original_photos_genus_estimable"] for q in z["models"].values())
    assert "GENUS_GEO_ALL_CLIMATE_SOIL" not in a["models"]
    assert "GENUS_GEO_ALL_CLIMATE_SOIL" in b["models"]


def test_failed_flower_roi_status_is_not_monomorphic_white(original_photo_species):
    original_photo_species.loc[0,"measurement_status"]="roi_unavailable"
    _,stat=M.source_sets(original_photo_species,strict=False)
    assert stat["n_classifiable_original_photos"]==399
    assert stat["n_unclassifiable_original_photos"]==801


def test_duplicated_source_species_identity_fails_closed(original_photo_species):
    original_photo_species.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="replicated"):
        M.source_sets(original_photo_species,strict=False)


def test_end_to_end_control_does_not_measure_genetic_colour(original_photo_species):
    d=M.run(original_photo_species,strict=False)
    assert d["schema"]=="fcp_global42111_colour_classifiability_moisture_negative_control_v1"
    assert d["original_source"]["n_classifiable_original_photos"]==400
    assert d["source_four_colour_labels_not_remeasured"]
    assert d["never_treat_unclassifiable_photographs_as_white_or_genetic_monomorphs"]
    assert d["climate_opportunity"]["outcome"].startswith("whether original photograph")
