"""Synthetic true-site 100/250/500-km congeneric photo geometry guards."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_local_congener_diameters_20261009 as D
import audit_fcp_global42111_local_congeners_20261009 as M


@pytest.fixture
def source():
    rng=np.random.default_rng(920)
    row=[]
    for g in range(24):
        centre=-170+20*(g%10)
        for i in range(24):
            jitter=(0.08 if g<16 else 4.0)
            row.append({
                "inat_taxon_id":g*24+i+1,
                "species":f"Genus{g} sp{g*24+i}",
                "latitude":12+rng.uniform(-jitter,jitter),
                "longitude":centre+rng.uniform(-jitter,jitter),
                "morph":M.CLASSES[(i+g)%4],
                "measurement_status":"classified_four_state_morph",
                "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
                "environment_climate_complete":True,
                "environment_soil_complete":i%7!=0,
                "wc_elevation_m":rng.uniform(0,600),
            })
    z=pd.DataFrame(row)
    for v in M.TEMP+M.RAIN+M.SOIL:
        z[v]=rng.uniform(1,10,len(z))
    z.loc[z.environment_soil_complete.eq(False),list(M.SOIL)]=np.nan
    return z


def test_distance_is_great_circle_km_not_cell_centroid():
    assert D.maximum_great_circle_km([0,0],[0,0])==pytest.approx(0,abs=1e-3)
    assert D.maximum_great_circle_km([0,0],[0,1])==pytest.approx(111.195,abs=0.02)
    with pytest.raises(ValueError,match="Unverified"):
        D.maximum_great_circle_km([12,float("nan")],[10,20])


def test_photos_further_apart_than_500km_are_excluded_but_local_groups_kept(source):
    pops,_=M.photo_populations(source,strict=False)
    geo=D.group_geometry(pops["CLIMATE_ALL"])
    assert len(geo)==24
    assert geo["n_distinct_original_species"].eq(24).all()
    assert int(geo.max_source_photo_distance_km.le(100).sum())==16
    assert int(geo.max_source_photo_distance_km.gt(500).sum())==8


def test_source_group_selection_cannot_depend_on_colour_labels(source):
    pops,_=M.photo_populations(source,strict=False)
    before=D.group_geometry(pops["CLIMATE_ALL"])
    other=pops["CLIMATE_ALL"].copy()
    other["morph"]="white"
    after=D.group_geometry(other)
    pd.testing.assert_frame_equal(before,after)


def test_all_predeclared_100_250_500_band_support_heldout_safe(source,monkeypatch):
    pops,_=M.photo_populations(source,strict=False)
    calls=[]
    def fake_test_local(chosen,*,include_soil):
        calls.append((len(chosen),include_soil))
        support=M.source_group_coverage(chosen)
        return {"status":"PRECOMMITTED_SOURCE_GEOMETRY_PREFLIGHT",
                "original_photo_support":support}
    monkeypatch.setattr(D,"test_local_congeners",fake_test_local)
    out,groups=D.assess(pops["CLIMATE_ALL"],include_soil=False)
    assert out["distance_cutoff_km_prespecified_all_reported"]==[100,250,500]
    assert set(out["distance_bounded_studies"])=={"100","250","500"}
    assert out["n_original_complete_case_photo_species"]==576
    assert out["n_original_genus_cell_groups_min3"]==24
    assert out["distance_bounded_studies"]["100"]["n_original_source_species"]==384
    assert out["distance_bounded_studies"]["500"]["n_original_source_species"]==384
    assert calls==[(384,False)]*3
    assert len(groups)==24


def test_insufficient_under_stricter_threshold_retains_hold_not_rescued(source):
    pops,_=M.photo_populations(source,strict=False)
    old=D.THRESHOLDS_KM
    try:
        D.THRESHOLDS_KM=(0.01,)
        out,_=D.assess(pops["CLIMATE_ALL"],include_soil=False)
        z=out["distance_bounded_studies"]["0"]
        assert z["n_original_source_species"]==0
        assert z["model"]["status"]=="HOLD_INSUFFICIENT_SAME_GENUS_SAME_CELL_SPECIES"
    finally:
        D.THRESHOLDS_KM=old
