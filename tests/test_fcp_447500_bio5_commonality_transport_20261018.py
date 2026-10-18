"""Within-species BIO5 replication cannot overwrite frozen climate, W or FCP claims."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.analysis import fcp_447500_bio5_commonality_transport_20261018 as f


def test_signed_summary_species_equal_not_photo_pseudoreplicates():
    x=f.summarise_deltas([.2,.8,-.1,0])
    assert x["n_species"]==4
    assert x["median_white_minus_chromatic_BIO5_within_species_SD"]==pytest.approx(.1)
    assert x["mean_white_minus_chromatic_BIO5_within_species_SD"]==pytest.approx(.225)
    assert x["fraction_positive_among_all_evaluable_species"]==pytest.approx(.5)
    assert x["n_exact_zero"]==1
    assert 0<=x["p_wilcoxon"]<=1 and 0<=x["p_sign_two_sided"]<=1
    assert f.summarise_deltas([])["estimable"] is False
    assert f.summarise_deltas([0,0,0])["p_wilcoxon"]==1.


def test_observer_and_month_groups_need_both_morph_states():
    d=pd.DataFrame([
        {"inat_taxon_id":1,"observer_id":5,"calendar_yearmonth":"2026-04","white":1,"bio5_z":1.},
        {"inat_taxon_id":1,"observer_id":5,"calendar_yearmonth":"2026-04","white":0,"bio5_z":-1.},
        {"inat_taxon_id":1,"observer_id":6,"calendar_yearmonth":"2026-05","white":1,"bio5_z":.8},
        {"inat_taxon_id":1,"observer_id":6,"calendar_yearmonth":"2026-05","white":0,"bio5_z":-.2},
        {"inat_taxon_id":2,"observer_id":7,"calendar_yearmonth":"2026-04","white":1,"bio5_z":.4},
    ])
    o,n=f.paired_species_strata(d,("inat_taxon_id","observer_id"))
    assert n==2 and len(o)==1
    assert o.loc[1]==pytest.approx(1.5)
    m,mm=f.paired_species_strata(d,("inat_taxon_id","calendar_yearmonth"))
    assert mm==2 and m.loc[1]==pytest.approx(1.5)


def test_equal_genus_sampled_species_independence():
    x=pd.Series([.4,.2,-.1],index=[1,2,3])
    group=pd.Series(["Silene","Silene","Ipomoea"],index=[1,2,3])
    s=f.equal_genus(x,group)
    assert s["n_species"]==2
    assert s["mean_white_minus_chromatic_BIO5_within_species_SD"]==pytest.approx(.1)


def test_background_white_fraction_can_be_missing_without_imputation():
    row={"background_palette_count_white":3}
    for col in f.SOURCE_COLUMNS:
        if col.startswith("background_palette_count_") and col not in row:
            row[col]=0
    a=pd.DataFrame([row,{col:0 for col in row}])
    x=f.background_white_fraction(a)
    assert x.loc[0,"background_white_fraction"]==pytest.approx(1.)
    assert np.isnan(x.loc[1,"background_white_fraction"])


def test_frozen_negative_replication_and_original_source_provenance_pins():
    assert f.ORIGINAL_PHOTO_SHA256=="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
    assert f.BIO5_RASTER_SHA256=="32e02d85868734c32547da70cf9c620f4739a3cf5dbc4f2b751e16ebb2338ae1"
    assert f.BIO5_ARCHIVE_SHA256=="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
    assert f.SPECIES_PHOTOS==447500 and f.CLASSIFIABLE_PHOTOS==227107
    assert "not" in f.__doc__.lower()
    assert "original 149900-photo prospective white-environment" in f.__doc__.lower()
