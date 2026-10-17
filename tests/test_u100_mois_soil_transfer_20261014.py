"""Synthetic-only contract tests for mois/FCP post hoc cross-source soil transfer."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SOURCE = Path(__file__).resolve().parents[1] / "scripts" / "analysis" / "u100_mois_soil_transfer_20261014.py"
# In staged local tree, script resides directly under scripts/.
if not SOURCE.exists():
    SOURCE=Path(__file__).resolve().parents[1]/"scripts"/"u100_mois_soil_transfer_20261014.py"
spec=importlib.util.spec_from_file_location("soil_transfer", SOURCE)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def source():
    d=pd.DataFrame([
        {"inat_taxon_id":1,"photo_id":"p1","latitude":0.0,"longitude":0.0,"positional_accuracy_m":20,"positional_accuracy_m":20,"soil_pH":5.0},
        {"inat_taxon_id":2,"photo_id":"p2","latitude":1.0,"longitude":1.0,"positional_accuracy_m":20,"soil_pH":np.nan},
        {"inat_taxon_id":3,"photo_id":"p3","latitude":0.0,"longitude":179.999,"positional_accuracy_m":20,"soil_pH":6.5},
        {"inat_taxon_id":4,"photo_id":"p4","latitude":0.0,"longitude":-179.999,"positional_accuracy_m":20,"soil_pH":7.0},
    ])
    for var in module.FEATURES:
        if var not in d:
            d[var]=12.5
    return d


def photos():
    return pd.DataFrame([
        {"inat_taxon_id":1,"photo_id":"p1","species":"A alba","observer_id":"o1","latitude":0.0,"longitude":0.0,"positional_accuracy_m":20,"white":True},
        {"inat_taxon_id":99,"photo_id":"new1","species":"B alba","observer_id":"o2","latitude":0.0005,"longitude":0.0,"positional_accuracy_m":20,"white":False},
        {"inat_taxon_id":99,"photo_id":"new2","species":"B alba","observer_id":"o3","latitude":1.0004,"longitude":1.0,"positional_accuracy_m":20,"white":True},
        {"inat_taxon_id":99,"photo_id":"new3","species":"B alba","observer_id":"o3","latitude":0.0,"longitude":-179.9991,"positional_accuracy_m":20,"white":False},
        {"inat_taxon_id":99,"photo_id":"new4","species":"B alba","observer_id":"o3","latitude":78.0,"longitude":20.0,"positional_accuracy_m":20,"white":True},
    ])


def test_mois_36_model_matches_exact_source_axes():
    assert len(module.FEATURES)==36
    assert len(module.SOIL)==10
    assert "wc_vapr_annual_mean" in module.FEATURES
    assert "soil_available_water_proxy" in module.FEATURES


def test_nearest_site_mapping_exact_and_dateline_and_not_filling_missing_soil():
    d=module.nearest_site_map(photos(),source())
    assert d.iloc[0].source_match_exact_photo
    assert d.iloc[0].distance_to_source_km==pytest.approx(0)
    assert d.iloc[1].distance_to_source_km < .1
    assert pd.isna(d.iloc[2].soil_pH)  # nearest source has missing soil: don't jump to another
    assert d.iloc[3].distance_to_source_km < .1  # antimeridian geometry
    assert not d.iloc[4].transfer_5km
    assert d.iloc[4].source_index >= 0


def test_exact_id_with_conflicting_taxon_fails():
    d=photos().copy();d.loc[0,"inat_taxon_id"]=999
    with pytest.raises(ValueError,match="conflicting taxon"):
        module.nearest_site_map(d,source())


def test_proximity_selection_uses_source_sites_not_outcome_or_soil_coverage():
    d=source()
    d.loc[d.photo_id.eq("p1"),"soil_pH"]=np.nan
    m=module.nearest_site_map(photos(),d)
    assert pd.isna(m.loc[m.photo_id.eq("new1"),"soil_pH"]).all()
    assert not m.loc[m.photo_id.eq("new1"),"soil_pH"].eq(7.0).any()


def test_holm_includes_full_36_family_when_not_estimable():
    p=np.array([.001, .02, np.nan]+[1.0]*33)
    a=module.holm_adjust(p)
    assert len(a)==36
    assert a[0]==pytest.approx(.036)
    assert a[1]>=.05
    assert a[2]==1


def test_species_effects_require_each_colour_and_same_observer_cell_support():
    r=[]
    for i in range(15):
        r.append({"inat_taxon_id":21,"species":"S alba","observer_id":f"o{i%3}",
                  "latitude":35.0,"longitude":135.0,"white":i<7,"soil_pH":float(i%4)})
    d=pd.DataFrame(r)
    a=module.per_species_effect(d,"soil_pH")
    assert len(a)==1
    assert a.iloc[0].paired_observers==3
    assert a.iloc[0].paired_1degree_cells==1
    assert np.isfinite(a.iloc[0].delta_white_minus_nonwhite_sd)
    d.loc[d.white,"soil_pH"]=np.nan
    assert module.per_species_effect(d,"soil_pH").empty


def test_empty_band_outputs_not_estimable_not_success(tmp_path):
    d=photos().iloc[:2].copy()
    receipt=module.run_analysis(d,source(),tmp_path)
    assert receipt["num_soil_variables"]==10
    assert receipt["num_total_variables"]==36
    assert receipt["holm_1km_supported"]==0
    tab=pd.read_csv(tmp_path/"mois_36_environment_transfer_screen.csv")
    assert tab.shape[0]==108
    assert receipt["holm_exact_supported"]==0
    assert tab.status.eq("not_estimable").all()


def test_target_accuracy_gate_does_not_approximate_uncertain_coordinates():
    d=photos().copy()
    d.loc[d.photo_id.eq("new1"),"positional_accuracy_m"]=3000
    mapped=module.nearest_site_map(d,source())
    assert not mapped.loc[mapped.photo_id.eq("new1"),"transfer_5km"].any()
    assert pd.isna(mapped.loc[mapped.photo_id.eq("new1"),"distance_to_source_km"]).all()


def test_vectorized_matches_reference_species_effect():
    rows=[]
    for i in range(20):
        rows.append({"inat_taxon_id":55,"species":"S test","observer_id":f"o{i%4}",
                     "latitude":35.0+i*.005,"longitude":135.0,"white":i<8,
                     "soil_pH":float(i%7),"soil_SOC":float(i//3)})
    d=pd.DataFrame(rows)
    a=module.per_species_effect(d,"soil_pH").iloc[0]
    b=module.all_species_effects(d,("soil_pH","soil_SOC"))
    v=b.set_index("variable").loc["soil_pH"]
    for field in ("delta_white_minus_nonwhite_sd","paired_observers",
                  "observer_paired_delta","paired_1degree_cells","same_cell_delta"):
        assert v[field]==pytest.approx(a[field])


def test_missing_observer_is_not_an_observer_pair():
    rows=[]
    for i in range(12):
        rows.append({"inat_taxon_id":7,"species":"S good", "observer_id":None,
                     "latitude":36.,"longitude":135.,"white":i<6,
                     "soil_pH":float(i)})
    r=module.all_species_effects(pd.DataFrame(rows),("soil_pH",))
    assert len(r)==1
    assert r.iloc[0].paired_observers==0
    assert pd.isna(r.iloc[0].observer_paired_delta)
