"""Synthetic geometric and null checks for original-source site-thinning audit."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

p=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_polymorphism_original_site_thinning_20261010.py"
spec=importlib.util.spec_from_file_location("site_thinning",p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def photo_table(n=60):
    return pd.DataFrame({
        "inat_taxon_id":[82]*n,"species":["X test"]*n,
        "photo_id":np.arange(n).astype(str),
        "latitude":[0.]*n,
        "longitude":np.arange(n)*0.02,
        "morph":["white"]*(n//2)+["red_pink"]*(n-n//2),
        "observer_id":np.arange(n)//2,
    })

def test_min_spacing_does_not_reuse_nearby_original_photo_sites():
    g=photo_table()
    dist=m.orig.pairwise_geo_km(g.latitude,g.longitude)
    for spacing in (1.,5.,10.):
        order=np.arange(len(g))
        idx=m.greedy_min_separation(dist,spacing,order)
        assert len(idx)>=2 and len(np.unique(idx))==len(idx)
        u,v=np.triu_indices(len(idx),1)
        assert np.all(dist[idx[u],idx[v]]>spacing)

def test_matching_selection_is_colour_blind():
    g=photo_table()
    d=m.orig.pairwise_geo_km(g.latitude,g.longitude)
    first=m.greedy_min_separation(d,5.,np.random.default_rng(m.seed("discovery",82,5.,0,"photo_selection")).permutation(len(g)))
    g.loc[:,"morph"]=g.morph.iloc[::-1].to_numpy()
    second=m.greedy_min_separation(d,5.,np.random.default_rng(m.seed("discovery",82,5.,0,"photo_selection")).permutation(len(g)))
    assert np.array_equal(first,second)

def test_thinning_makes_exactly_original_four_label_null_without_new_species():
    g=photo_table()
    d=m.orig.pairwise_geo_km(g.latitude,g.longitude)
    labels=pd.Categorical(g.morph,categories=m.orig.MORPHS).codes
    r,n=m.study_one_realization(g,d,labels,"discovery",82,1.,0)
    assert n==len(g)
    for policy in m.OBSERVER_POLICIES:
        x=r[policy]
        assert x is not None
        assert len(x["null_depletion"])==m.N_PERM
        assert abs(x["null_depletion"].mean())<.12

def test_ten_km_source_thinning_can_restrict_opportunity_without_creating_false_zero():
    g=photo_table()
    d=m.orig.pairwise_geo_km(g.latitude,g.longitude)
    y=pd.Categorical(g.morph,categories=m.orig.MORPHS).codes
    r,n=m.study_one_realization(g,d,y,"third",82,10.,0)
    assert n<m.MIN_RETAINED_PHOTOS
    assert all(r[k] is None for k in m.OBSERVER_POLICIES)

def test_species_summary_holds_and_does_not_count_thinning_repeats_as_species():
    records=[]
    nulls=[]
    for i in range(34):
        row={"cohort":"discovery","inat_taxon_id":i}
        n={}
        for dist in m.SPACINGS_KM:
            for policy in m.OBSERVER_POLICIES:
                key=f"{int(dist)}km_{policy}"
                row[f"{key}_depletion"]=.07 if i<31 else np.nan
                row[f"{key}_n_retained_median_all_repeats"]=37
                row[f"{key}_n_photo_pairs_median_valid"]=20
                n[(int(dist),policy)]=np.zeros(m.N_PERM) if i<31 else None
        records.append(row);nulls.append(n)
    out=m.summarize_source(pd.DataFrame(records),nulls)
    assert out["5km_all_photos"]["n_species"]==31
    assert out["5km_all_photos"]["positive_bounded_evidence"] is True
    assert np.allclose(out["5km_all_photos"]["species_bootstrap_ci95"],[.07,.07])
    assert out["10km_different_observer"]["n_species"]==31
