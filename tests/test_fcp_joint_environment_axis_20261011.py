from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'/'analysis'))
from fcp_joint_environment_axis_20261011 import components, paired_summary


def test_exact_orthogonal_loss_decomposition():
    rng=np.random.default_rng(14)
    y=rng.dirichlet(np.ones(4),size=40)
    p=rng.dirichlet(np.ones(4),size=40)
    out=components(y,p)
    assert np.allclose(out[:,0]+out[:,1],out[:,2],atol=1e-14)
    assert np.allclose(out[:,0],4/3*(p[:,0]-y[:,0])**2)
    q=np.array([3,-1,-1,-1])/np.sqrt(12)
    assert np.allclose(out[:,0],((p-y)@q)**2)
    assert np.all(out>=-1e-13)


def test_colour_only_change_has_no_white_axis():
    y=np.array([[.5,.2,.2,.1]])
    p=np.array([[.5,.3,.1,.1]])
    a=components(y,p)
    assert a[0,0]==0
    assert a[0,1]>0


def test_white_only_shift_has_no_remaining_hue():
    y=np.tile([.5,1/6,1/6,1/6],(30,1))
    p=np.tile([.2,.8/3,.8/3,.8/3],(30,1))
    a=components(y,p)
    assert np.allclose(a[:,1],0,atol=1e-13)


def test_normalization_removes_proportional_improvement_illusion():
    rng=np.random.default_rng(35)
    y=rng.dirichlet([8,1,1,1],size=80)
    p=np.tile([.2,.4,.2,.2],(80,1))
    p1=y+np.sqrt(.9)*(p-y)
    a=components(y,p);b=components(y,p1)
    s=paired_summary(a,b,np.repeat(np.arange(20),4),199,8)
    for weight in ('equal_genus','equal_species'):
        r=s[weight]
        assert np.isclose(r['white_axis']['relative_error_reduction'],.1)
        assert np.isclose(r['remaining_hue']['relative_error_reduction'],.1)
        assert np.isclose(r['axis_minus_hue_relative_skill']['estimate'],0,atol=1e-12)


def test_identical_models_cannot_create_environmental_gain():
    rng=np.random.default_rng(24)
    y=rng.dirichlet(np.ones(4),size=36);p=rng.dirichlet(np.ones(4),size=36)
    a=components(y,p)
    z=paired_summary(a,a,np.repeat(np.arange(12),3),49,1)
    assert z['equal_genus']['total']['gain']==0
    assert z['equal_genus']['white_axis']['gain_ci95']==[0.,0.]
    assert z['equal_genus']['white_share_of_total_gain'] is None


def test_zero_baseline_component_is_not_imputed_as_skill():
    y=np.tile([.5,.2,.2,.1],(40,1));p=np.tile([.5,.3,.1,.1],(40,1))
    z=paired_summary(components(y,p),components(y,y),np.repeat(np.arange(10),4),19,1)
    assert z['equal_genus']['white_axis']['relative_error_reduction'] is None
    assert z['equal_genus']['axis_minus_hue_relative_skill']['estimate'] is None


def test_invalid_simplex_rejected():
    with pytest.raises(ValueError):components(np.array([[1,0,0,0]]),np.array([[.5,.5,.5,.5]]))
    with pytest.raises(ValueError):components(np.ones((2,3))/3,np.ones((2,3))/3)


def test_species_weighting_and_genus_weighting_are_different_targets():
    y=np.tile([.5,.2,.2,.1],(12,1))
    p=np.tile([.2,.3,.3,.2],(12,1))
    p1=p.copy();p1[:2]=y[:2]
    out=paired_summary(components(y,p),components(y,p1),np.array(['a']*2+['b']*10),49,2)
    assert out['equal_genus']['total']['gain']>out['equal_species']['total']['gain']
