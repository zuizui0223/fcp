"""Exact dated source-phylogeny covariance and 50/250km spatial kernel safety."""
from __future__ import annotations
import sys
from io import StringIO
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_real_spatial_phylo_covariance_20261010 as M


def test_true_Brownian_shared_branch_not_genus_dummy():
    tree=Phylo.read(StringIO("((A_a:1,A_b:1):1,B_c:2);"),"newick")
    K=M.BM_shared_tree_cov(tree,["A_a","A_b","B_c"])
    assert K.shape==(3,3)
    assert K[0,0]==pytest.approx(1)
    assert K[0,1]==pytest.approx(.5)
    assert K[0,2]==pytest.approx(0)
    assert np.linalg.eigvalsh(K).min()>=-1e-8
    with pytest.raises(ValueError,match="exactly"):
        M.BM_shared_tree_cov(tree,["A_a","A_b","Other_sp"])


def test_true_geodesic_kernel_decays_with_photo_distance():
    points=pd.DataFrame({"latitude":[0.,0.,0.],"longitude":[0.,1.,45.]})
    a=M.geographic_kernel(points,50)
    b=M.geographic_kernel(points,250)
    assert np.allclose(np.diag(a),1)
    assert a[0,2]<a[0,1]<1
    assert b[0,1]>a[0,1]
    assert np.linalg.eigvalsh(a).min()>-1e-7


@pytest.fixture
def sampled(monkeypatch):
    monkeypatch.setattr(M,"MIN_SOURCE",100)
    n=160
    rng=np.random.default_rng(20261010)
    i=np.arange(n)
    data=pd.DataFrame({
        "inat_taxon_id":i+1,
        "species":[f"G{i%40} sp{i}" for i in i],
        "genus":[f"G{j%40}" for j in i],
        "morph":np.array(M.CLASSES)[i%4],
        "latitude":rng.uniform(-60,60,n),
        "longitude":rng.uniform(-178,178,n),
    })
    for f in (v for block in M.BLOCKS.values() for v in block):
        data[f]=rng.normal(size=n)
    data=M.add_source_geo(data)
    Kphy=np.eye(n)*.7+np.ones((n,n))*.3
    return data,Kphy


def test_abiotic_increment_only_after_both_phylo_and_spatial_kernel(sampled):
    data,Kphy=sampled
    Ksp=M.geographic_kernel(data,250)
    for factor in ("genus","true_photo_cell"):
        out=M.fixed_krr_scores(data,Ksp,Kphy,factor,soil=False,nboot=19)
        assert out["status"]=="EXPLORATORY_EXACT_PHYLO_BROWNIAN_KERNEL_PLUS_SPATIAL_KERNEL_ENVIRONMENT"
        assert out["same_original_species_and_fold_assignments"]
        assert out["modelled_true_phylogeny_BM_shared_branch_covariance"]
        assert out["modelled_original_photograph_spatial_exponential_covariance"]
        assert "all_abiotic_beyond_spatial_plus_real_phylogeny" in out["conditional_feature_gains"]
        assert "unique_precipitation_beyond_space_real_phylogeny_and_other_abiotic" in out["conditional_feature_gains"]
        assert all(len(v["original_group_cluster_95CI"])==2 for v in out["conditional_feature_gains"].values())
        assert {v["n_same_species"] for v in out["source_model_scores"].values()}=={160}


def test_same_photo_source_complete_cases_not_reclassified(sampled):
    d,_=sampled
    d.loc[0,"soil_pH"]=np.nan
    climate,c0=M.source_features(d,soil=False)
    soil,c1=M.source_features(d,soil=True)
    assert len(climate)==160
    assert len(soil)==159
    assert c1["n_direct_original_photo_taxa_excluded_by_environment_mask"]==1


def test_early_hold_on_too_few_real_direct_tips(sampled,monkeypatch):
    d,Kphy=sampled
    monkeypatch.setattr(M,"MIN_SOURCE",300)
    z=M.fixed_krr_scores(d,M.geographic_kernel(d,50),Kphy,"genus",soil=False,nboot=2)
    assert z["status"]=="HOLD_DIRECT_LCVP_MODEL_COVERAGE"
