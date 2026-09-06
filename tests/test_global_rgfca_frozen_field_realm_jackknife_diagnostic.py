import importlib.util
from pathlib import Path
import numpy as np
import pytest
P=Path(__file__).resolve().parents[1]/'scripts/analysis/run_global_rgfca_frozen_field_realm_jackknife_diagnostic.py'
s=importlib.util.spec_from_file_location('realm_diagnostic',P)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def test_weighted_variance_manual():
    assert np.isclose(m.weighted_variance(np.array([1.,2.,4.]),np.array([1.,2.,1.])),1.1875)

def test_variance_recentered_after_deletion():
    y=np.arange(6.)
    null=np.tile(np.arange(6.)/2,(999,1))
    labels=np.array(['Nearctic']*3+['Palearctic']*3)
    rows,b,n=m.deletion_rows(y,null,np.ones(6),np.ones(6,dtype=bool),labels)
    a=next(x for x in rows if x['deleted_group']=='Nearctic')
    assert np.isclose(a['observed_concentration'],2/3)
    assert a['observed_to_null_mean_ratio']==pytest.approx(4.)
    assert a['removed_opportunity_fraction']==.5
    assert b['reconstructed_existing_primary_p_upper']==.001

def test_empty_group_is_not_a_pass():
    rows,_,_=m.deletion_rows(np.arange(6.),np.tile(np.arange(6.),(999,1)),np.ones(6),np.ones(6,dtype=bool),['Nearctic']*6)
    a=next(x for x in rows if x['deleted_group']=='Oceania')
    assert a['status']=='not_evaluable_no_evaluable_cells'
    assert 'positive_excess' not in a
    b=next(x for x in rows if x['deleted_group']=='Nearctic')
    assert b['status']=='not_evaluable_fewer_than_three_retained_cells'

def test_nonfinite_and_nonpositive_weights_rejected():
    with pytest.raises(ValueError):m.weighted_variance([1.,np.nan,2.],[1.,1.,1.])
    with pytest.raises(ValueError):m.weighted_variance([1.,2.,3.],[1.,0.,1.])

def test_null_count_must_remain_frozen():
    with pytest.raises(ValueError):m.deletion_rows(np.arange(6.),np.ones((99,6)),np.ones(6),np.ones(6,dtype=bool),['Nearctic']*6)
