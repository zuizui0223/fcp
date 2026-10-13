import math
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from within_line import (load_source, line_score, B_and_line_null, informative_exposure_line,
                         distance_residual, site_conditioned, plusone, holm_two)


def test_source_is_exact_512_field_rows():
    d=load_source(Path(__file__).resolve().parents[1]/'source'/'frequency.txt')
    assert len(d)==512
    assert d['v'].sum()==5081 and d['w'].sum()==8723
    assert d['l'].sum()==4463 and d['s'].sum()==29
    assert ((d['n']==d['v']+d['w'])&(d['n']>=0)).all()
    assert set(d['year'])=={'Y11','Y13'}


def test_line_score_is_colourwise_congener_share_difference():
    n=np.array([2,2]);v=np.array([2,0]); a=np.array([.5,0.])
    assert math.isclose(line_score(n,v,a),.5,abs_tol=1e-12)
    assert math.isclose(line_score(n,np.array([1,1]),a),0,abs_tol=1e-12)


def test_line_randomisation_expectation_and_partition():
    n=np.array([2,2]);v=np.array([2,0]);line=['A','A']
    b, between, enull=B_and_line_null(n,v,line)
    assert math.isclose(b,1,abs_tol=1e-12)
    assert math.isclose(between,0,abs_tol=1e-12)
    assert math.isclose(enull,1/3,abs_tol=1e-12)


def test_between_line_only_difference_is_not_within_line_signal():
    n=np.array([2,2,2,2]);v=np.array([2,2,0,0]);line=['A','A','B','B']
    b, between, expected=B_and_line_null(n,v,line)
    assert math.isclose(b,1,abs_tol=1e-12)
    assert math.isclose(between,1,abs_tol=1e-12)
    assert math.isclose(expected,1,abs_tol=1e-12)


def test_informative_line_requires_two_occupied_plots_color_overlap_and_exposure_variation():
    n=np.array([5,5]);v=np.array([1,4]);a=np.array([.4,.2])
    assert informative_exposure_line(n,v,a)
    assert not informative_exposure_line(np.array([10]),np.array([5]),np.array([.3]))
    assert not informative_exposure_line(n,np.array([0,0]),a)
    assert not informative_exposure_line(n,v,np.array([.2,.2]))


def test_distance_residual_eliminates_exact_linear_trend():
    n=np.array([2,3,5]);d=np.array([1,2,4]);a=np.array([.1,.2,.4])
    r=distance_residual(n,d,a)
    assert np.max(np.abs(r))<1e-12


def test_exact_null_matches_multivariate_hypergeom():
    n=np.array([2,2]);v=np.array([2,0]);a=np.array([.5,0.]);line=['A','A']
    out=site_conditioned(n,v,a,line,R=3999,seed=2026101309,min_line_n=2)
    assert out['n_replicates']==3999
    assert abs(np.mean(out['T_null']))<.03
    assert abs(np.mean(out['B_null'])-(1/3))<.03
    assert math.isclose(out['T_observed'],.5,abs_tol=1e-12)


def test_plusone_holm_guard():
    assert plusone(-1,np.arange(20),upper=True)==1.0
    assert plusone(5,np.zeros(9999),upper=True)==1/10000
    assert holm_two(.012,.034)==(.024,.034)


def test_line_synthetic_exact_exposure_by_original_other_species_counts():
    # The 4 white L-species individuals are NOT white M-species individuals.
    n=np.array([5,5]);v=np.array([4,1]);L=np.array([4,0]);S=np.array([0,0]);a=(L+S)/(L+S+n)
    assert line_score(n,v,a)>0
    assert sum(n)==10 and sum(L+S+n)==14


def test_line_conditional_baseline_respects_line_inventory_even_if_extreme():
    n=np.array([5,5,10,10]);v=np.array([5,0,10,0]);line=['a','a','b','b']
    b,bb,expected=B_and_line_null(n,v,line)
    assert 0<=bb<=expected<b<=1+1e-12
    assert math.isclose(b,1,abs_tol=1e-12)


def test_absent_target_plot_does_not_enter_target_colour_denominator():
    n=np.array([2,2,0]);v=np.array([2,0,0]);line=['A','A','B']
    b,between,expected=B_and_line_null(n,v,line)
    assert math.isclose(b,1,abs_tol=1e-12)
    assert math.isclose(expected,1/3,abs_tol=1e-12)


def test_exact_source_total_stays_separate_from_2011_2013_siteyear_variation():
    d=load_source(Path(__file__).resolve().parents[1]/'source'/'frequency.txt')
    assert d.groupby(['year','population']).ngroups==16
    assert d['population'].nunique()==13
    assert (d['n']==0).sum()==126
    assert d[['v','w']].sum().sum()==13804
    assert d[['l','s']].sum().sum()==4492
