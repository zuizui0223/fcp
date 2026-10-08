"""Synthetic phenotype ecology: no false imputations, species equality, soil increments."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import characterize_fcp_global42111_colour_climate_soil_20261009 as M


@pytest.fixture
def global_source():
    n=42111
    rng=np.random.default_rng(20261009)
    i=np.arange(n)
    x=pd.DataFrame({
        "inat_taxon_id":i+1,
        "species":[f"Genus{i%260} species{i}" for i in i],
        "latitude":rng.uniform(-60,60,n),
        "longitude":rng.uniform(-179,179,n),
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
    })
    soil=rng.standard_normal(n)
    cls=np.resize(np.array(M.COLOURS),n)
    x["morph"]=np.where(i<18457,cls,"UNKNOWN")
    x["measurement_status"]=np.where(i<18457,"classified_four_state_morph","unclassified")
    # 18457 classified photos but only first 1600 have all soil+climate data.
    x["environment_climate_complete"]=i<2000
    x["environment_soil_complete"]=i<1600
    for col in M.GEO[3:]+M.CLIM+M.SOIL:
        v=soil + rng.normal(0,1,n)
        x[col]=v
    missing=i>=1600
    x.loc[missing,M.SOIL]=np.nan
    x.loc[i>=2000,M.CLIM]=np.nan
    return x


def test_preserve_full_species_and_incomplete_photo_denominator(global_source):
    d,a=M.photo_equal_frame(global_source)
    assert len(global_source)==42111
    assert a["original_colour_classifiable_taxa"]==18457
    assert a["source_taxa"]==42111
    assert a["geo_climate_soil_eligible_classified"]==1600
    assert len(d)==1600 and d.inat_taxon_id.nunique()==1600
    assert d.site_cell_162.between(0,161).all()
    assert d.genus.nunique()>200


def test_genus_heldout_and_spatial_heldout_identical_source_population(global_source):
    d,_=M.photo_equal_frame(global_source)
    genus=M.split_loss(d,"genus")
    cell=M.split_loss(d,"site_cell_162")
    for output in (genus,cell):
        assert output["heldout_species"]==1600
        assert all(q["n_species_same_tested"]==1600 for q in output["feature_families"].values())
        assert all(np.isfinite(q["pooled_heldout_log_loss"]) for q in output["feature_families"].values())
        assert len(output["feature_families"]["GEO_CLIMATE_SOIL"]["heldout_fold_diagnostics"])==5


def test_full_results_are_noncausal_no_label_replacement(global_source):
    report,profiles=M.run(global_source)
    assert report["n_full_abiotic_colour_sample"]==1600
    assert report["global_original_species_denominator"]==42111
    assert report["photo_colour_outcomes_reclassified"] is False
    assert report["identical_complete_case_rows_for_every_model"] is True
    assert report["spatial_cell_centroid_as_plant_soil"] is False
    assert len(profiles)==4*(1+len(M.CLIM)+len(M.SOIL))


def test_unclassifiable_photos_never_imputed_as_white(global_source):
    global_source.loc[20000,"morph"]="white"
    # This photo has a white string but no classified status, thus never included.
    out,a=M.photo_equal_frame(global_source)
    assert len(out)==1600
    assert a["original_colour_classifiable_taxa"]==18457


def test_missing_site_geo_is_not_rescued_by_value_synthesis(global_source):
    global_source.loc[0,"site_geo_status"]="NO_UNOBSCURED_PUBLIC_COORDINATES"
    out,a=M.photo_equal_frame(global_source)
    assert out.inat_taxon_id.min()==2
    assert a["photo_classifiable_but_geolocated_missing"]==1


def test_original_all_42111_photo_opportunity_coverage_regions(global_source):
    x=M.region_coverage(global_source)
    assert int(x.source_species.sum())==42111
    assert int(x.classified_flower_photo.sum())==18457
    assert int(x.all_environment_complete_species.sum())==1600
    assert set(x.source_latitude_region).issubset({"0_30","30_60","60_90","NO_EXACT_PUBLIC_GEO"})
