"""Original-photo whole-genus matched geography study: causal-inference guard tests."""
import importlib.util,sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

p=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_fcp_same_genus_east_elsewhere_rain_tradeoff_20261011.py"
sys.path.insert(0,str(p.parent))
spec=importlib.util.spec_from_file_location("fcp_shared_genus_source",p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def make_source():
    rows=[]
    for genus in range(42):
        for east in [False, True]:
            for number in range(3):
                row={"genus":f"G{genus}","inat_taxon_id":(genus*1000+100*int(east)+number),
                     "species":f"G{genus} species{number}",
                     "longitude":120. if east else -50.,"latitude":30.,
                     "photo_id":len(rows)}
                rows.append(row)
    return pd.DataFrame(rows)


def test_same_genus_source_filter_requires_independent_species_each_side():
    x=make_source()
    x=x[~((x.genus=="G1")&(x.longitude<0)&(x.inat_taxon_id%100==1))]
    selected,ng=m.shared_genus_subset(x,3)
    assert ng==41
    assert "G1" not in selected.genus.unique()
    assert selected.groupby(["genus","east"]).inat_taxon_id.nunique().min()==3


def test_distinct_taxa_single_side_cannot_identify_interaction():
    x=make_source()
    x=x.loc[x.longitude>0]
    with pytest.raises(ValueError,match="missing one geographic side"):
        m.shared_genus_subset(x,3)


def test_paired_same_genus_resampling_recovers_induced_white_slope_difference():
    rng=np.random.default_rng(13)
    a={};b={}
    for i in range(45):
        Z=rng.normal(size=(10,8))
        Z-=Z.mean(axis=0)
        mat=Z.T@Z
        xa=np.zeros((8,4));xa[0,:]=[.1,-.1,0.,0.]
        xb=np.zeros((8,4))
        a[f"G{i}"]=(mat,mat@xa,10)
        b[f"G{i}"]=(mat,mat@xb,10)
    v=m.paired_genus_bootstrap(a,b,repetitions=89,seed=12)
    assert v["status"]=="OK"
    assert v["n_shared_genera"]==45
    assert abs(v["response_four_original_photo_hues"]["white"]["east_minus_elsewhere_rainfall_slope"]-.1)<1e-11
    assert v["response_four_original_photo_hues"]["white"]["paired_same_genus_difference_CI95"][0]>.099


def test_no_generalizable_claim_for_fewer_than_30_shared_genera():
    A={str(i):(np.eye(8),np.ones((8,4)),3) for i in range(20)}
    out=m.paired_genus_bootstrap(A,A,repetitions=10,seed=1)
    assert out["status"]=="HOLD_FEWER_THAN_30_SHARED_GENERA"


def test_species_paired_east_west_requires_100_original_named_taxa():
    x=make_source()
    x["white"]=0.
    for col in m.COLS:x[col]=1.
    z=m.paired_same_species(x,np.zeros(8),np.ones(8),1)
    assert z["status"]=="HOLD_INSUFFICIENT_MATCHED_SOURCE_SPECIES"


def test_four_photo_hue_coefficients_sum_to_zero_with_fixed_source_probabilities():
    rng=np.random.default_rng(14)
    a={};b={}
    for i in range(42):
        Z=rng.normal(size=(14,8));Z-=Z.mean(axis=0)
        mat=Z.T@Z
        beta=np.zeros((8,4))
        beta[0,:]=[.04,-.02,-.01,-.01]
        a[f"G{i}"]=(mat,mat@beta,14)
        b[f"G{i}"]=(mat,np.zeros((8,4)),14)
    z=m.paired_genus_bootstrap(a,b,repetitions=19,seed=99)
    s=sum(z["response_four_original_photo_hues"][h]["east_rainfall_slope"] for h in m.HUES)
    assert abs(s)<1e-12


def test_study_has_explicit_noncausal_source_referent():
    assert "photo" in (m.__doc__ or "").lower()
    assert len(m.COLS)==8
    assert m.MIN_REGION_SPECIES == (3,5,10)
