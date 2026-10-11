"""Synthetic original photo hue-niche null validity and source population guards."""
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pytest

p=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_fcp_photo_white_colour_climate_niche_segregation_20261011.py"
sys.path.insert(0,str(p.parent))
spec=importlib.util.spec_from_file_location("fcp_hue_niche",p)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_source_sha_mismatch_fails_closed(tmp_path):
    f=tmp_path/"bad_source.zip";f.write_bytes(b"wrong")
    with pytest.raises(RuntimeError,match="SHA256"):mod.load_original(f)

def test_standardized_white_nonwhite_environment_centroids():
    x=np.array([[-1.,0.],[0.,1.],[1.,0.],[2.,2.]])
    assert np.allclose(mod.contrast(x,np.array([1,1,0,0],np.int8))[0],[-2.,-.5])

def test_species_photo_colour_count_unchanged_under_unrestricted_shuffle():
    label=np.array([1,1,0,0,1,0,0,0],np.int8)
    perm=mod.shuffle(label,None,np.random.default_rng(4),100)
    assert np.all(perm.sum(axis=1)==3)

def test_local_shuffle_preserves_exact_original_photo_colour_counts_per_cell():
    label=np.array([1,1,0,0,1,0,0,0],np.int8)
    blocks=np.array([0,0,0,1,1,1,1,1])
    perm=mod.shuffle(label,blocks,np.random.default_rng(4),100)
    for block in (0,1):
        idx=np.flatnonzero(blocks==block)
        assert np.all(perm[:,idx].sum(axis=1)==label[idx].sum())

def test_no_mix_block_gives_zero_informative_spatial_support():
    labels=np.array([1,1,0,0],np.int8)
    fake=[{"taxon":7,"genus":"Example","labels":labels,
           "cell10":np.array([1,1,2,2]),"cell20":np.array([1,1,2,2])}]
    z=mod.source_spatial_exchange(fake)
    assert z["10degree"]["n_species_with_both_original_photo_colours_in_one_geocell"]==0
    assert z["20degree"]["informative_source_photo_count"]==0

def test_original_10degree_photo_geography_within_bounds():
    loc=mod.geographic_cell(np.array([-90.,0.,90.]),np.array([-180.,0.,180.]),10)
    assert np.array_equal(loc,[0,9*37+18,18*37+36])

def test_no_difference_when_species_colour_groups_share_identical_environment():
    lab=np.array([1,1,0,0],np.int8)
    sp=[{"taxon":2,"genus":"Example","n":4,"n_white":2,"labels":lab,
         "z":np.zeros((4,4)),"cell10":np.array([0,0,1,1]),
         "cell20":np.array([0,0,1,1])}]
    z=mod.evaluate(sp,nrep=20)
    assert "results" not in z  # output is a thresholds dictionary, no eligible fake gate
