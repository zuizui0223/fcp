"""Synthetic tests for full 42111 measured source and actual coordinates."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_global42111_abiotic_join_preflight_20261009 import (
    coordinate_contract, audit
)


def dummy_frozen():
    s=np.arange(42111,dtype=int)+1
    breadth=pd.DataFrame({
        "inat_taxon_id":s,
        "morph":np.where(s<=18457,"white","UNCLASSIFIABLE"),
        "measurement_status":np.where(s<=18457,"classified_four_state_morph","not_classified"),
    })
    k=np.arange(85337,dtype=int)
    taxon=(k%42111)+1
    geo=pd.DataFrame({
        "inat_taxon_id":taxon,
        "cell_id":k//42111,
        "morph":np.where(k<39075,"red_pink",""),
        "measurement_status":np.where(k<39075,"classified_four_state_morph","unclassifiable"),
    })
    pairs=pd.DataFrame({
        "inat_taxon_id":s[:13416],
        "cell_id_1":0,"cell_id_2":1,
        "pair_state":"unclassifiable",
    })
    return breadth,geo,pairs


def test_genuine_coordinate_coverage_is_not_faked_by_grid_cell():
    f=pd.DataFrame({"cell_id":[0,1],"inat_taxon_id":[1,2]})
    d=coordinate_contract(f)
    assert d["status"]=="NO_VERIFIED_OBSERVATION_COORDINATE_COLUMNS"
    assert d["n_missing"]==2
    assert d["site_level_soil_or_climate_sampling_permitted"] is False


def test_valid_coordinates_with_unrepresentative_ocean_positions_are_only_eligible():
    f=pd.DataFrame({"latitude":[10,90.1,float("nan"),-80],
                    "longitude":[130,0,10,-10],"cell_id":[0,0,0,0]})
    d=coordinate_contract(f)
    assert d["n_valid"]==2
    assert d["n_missing"]==2
    assert d["site_level_soil_or_climate_sampling_permitted"] is True


def test_conflicting_location_aliases_fail_closed_on_sampling():
    f=pd.DataFrame({"latitude":[10,20],"longitude":[100,110],
                    "observed_latitude":[10,21],"observed_longitude":[100,110]})
    d=coordinate_contract(f)
    assert d["n_conflicting_coord_alias_rows"]==1
    assert not d["site_level_soil_or_climate_sampling_permitted"]


def test_exact_full_source_denominators_and_missing_classification():
    d=audit(*dummy_frozen())
    assert d["n_original_taxa"]==42111
    assert d["n_original_taxon_cell_records"]==85337
    assert d["n_original_observer_disjoint_pairs"]==13416
    assert d["n_breadth_colour_classified"]==18457
    assert d["n_taxon_cell_colour_classified"]==39075
    assert not d["global_climate_soil_covariates_computed"]
    assert d["layers"][1]["coordinate_audit"]["site_level_soil_or_climate_sampling_permitted"] is False


def test_duplicate_cell_species_never_counted_twice():
    breadth,geo,pairs=dummy_frozen()
    geo.loc[0,"inat_taxon_id"]=int(geo.loc[42111,"inat_taxon_id"])
    geo.loc[0,"cell_id"]=int(geo.loc[42111,"cell_id"])
    with pytest.raises(ValueError,match="taxon-cell"):
        audit(breadth,geo,pairs)


def test_unclassified_not_treated_as_white_or_monomorphic():
    breadth,geo,pairs=dummy_frozen()
    geo.loc[:39074,"measurement_status"]="not_classified"
    with pytest.raises(ValueError,match="classification"):
        audit(breadth,geo,pairs)
