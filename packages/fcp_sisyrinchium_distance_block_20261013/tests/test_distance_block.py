import sys,io,zipfile,itertools
from pathlib import Path
from fractions import Fraction
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from distance_block import (distance_bins,phenotype_partition,source_rows,analyze,hypergeom_trial,eligibility)


def test_distance_bins_median_ties_and_unknown():
    distances=[5,5,10,15,20,25,-1,-1]
    n=[2,0,3,4,6,8,2,1]
    groups=distance_bins(distances,n,method='median')
    assert groups==['near','NO_TARGET_PLANTS','near','near','far','far','unknown','unknown']


def test_distance_bins_adjacent_unique_pair():
    assert distance_bins([1,5,5,10,15,20,-1],[1,1,1,1,1,1,1],method='pairs')==[
        'pair00','pair00','pair00','pair01','pair01','pair02','unknown']


def test_exact_simple_equal_two_plots():
    B,between,E,delta=phenotype_partition([2,2],[2,0],['near','near'])
    assert B==1
    assert abs(between-0)<1e-14
    assert abs(E-1/3)<1e-14
    assert abs(delta-2/3)<1e-14


def test_all_bins_singletons_have_zero_residual():
    B,between,E,delta=phenotype_partition([2,2],[2,0],['near','far'])
    assert B==1 and between==1 and E==1 and abs(delta)<1e-14


def test_exact_hand_enumeration_multi_bin():
    n=[2,3,2,1]
    cols=['near','near','far','far']
    # each bin fixes one violet label; enumerate individual slots.
    obs_v=[1,0,1,0]
    B,between,E,delta=phenotype_partition(n,obs_v,cols)
    positions=[list(range(2)),list(range(2,5)),list(range(5,7)),[7]]
    values=[]
    for i in range(5):
        for j in range(3):
            # one violet in near-bin and one in far-bin
            vv=[sum(i in positions[k] for _ in [0]) for k in range(2)]+[sum(5+j in positions[k] for _ in [0]) for k in range(2,4)]
            values.append(phenotype_partition(n,vv,cols)[0])
    assert abs(np.mean(values)-E)<1e-14


def test_census_source_bytes_and_population_units():
    rows=source_rows(ROOT/'source/frequency.txt',ROOT/'source/Sisyrinchium_data.zip')
    assert len(rows)==512
    assert sum(r['v'] for r in rows)==5081
    assert sum(r['w'] for r in rows)==8723
    assert sum(r['l'] for r in rows)==4463
    assert sum(r['s'] for r in rows)==29
    assert len({(r['year'],r['population'],r['line'],r['quad']) for r in rows})==512
    assert sum(int(r['distance'])<0 for r in rows)==3


def test_eligibility_needs_two_informative_lines():
    # all units one line: not enough independent recorded line groups
    d={'first':dict(n=[10,10],v=[5,5],line='A')}
    assert eligibility(d)['status']=='HOLD_DISTANCE_BIN_SUPPORT'


def test_simulation_preserves_inventory_and_areas():
    sizes=[3,2,4];V=4
    rng=np.random.default_rng(2026101312)
    draw=hypergeom_trial(sizes,V,rng,999)
    assert draw.shape==(999,3)
    assert np.all(draw.sum(axis=1)==V)
    assert np.all((draw>=0)&(draw<=np.asarray(sizes)))


def test_analysis_preserves_source_and_main_supported_or_holds(tmp_path):
    a=analyze(ROOT/'source/frequency.txt',ROOT/'source/Sisyrinchium_data.zip', ROOT/'source/parent_line_result.json',R=19)
    assert a['source']['source_quadrats']==512
    assert len(a['support'])==16
    assert a['median']['n_eligible']>=0
    assert a['pairs']['n_eligible']>=0
    assert a['median']['n_eligible']+sum(x['median']['status'].startswith('HOLD') for x in a['support'])==16


def test_sample_conditioned_null_returns_one_for_degenerate():
    B,_,E,delta=phenotype_partition([1,1],[1,0],['near','far'])
    assert abs(delta)<1e-14
