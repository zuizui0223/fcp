"""Pre-result gates for NHM independent botanical colour vs photo-site rainfall."""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1]/'scripts/analysis/audit_independent_nhm_colour_rain_20261013.py'
spec=importlib.util.spec_from_file_location('fcp_nhm_independent',P)
m=importlib.util.module_from_spec(spec)
assert spec and spec.loader
sys.modules[spec.name]=m
spec.loader.exec_module(m)


def test_single_organ_is_not_white_when_mixed():
    assert m.classify_botanical_description({'petal colour':'white'})[0]=='white'
    assert m.classify_botanical_description({'corolla colour':'blue or violet'})[0]=='blue_purple'
    assert m.classify_botanical_description({'corolla colour':'white with purple veins'})[0]=='HOLD'
    assert m.classify_botanical_description({'corolla colour':'pink to white'})[0]=='HOLD'
    assert m.classify_botanical_description({'flower colour':'cream'})[0]=='OTHER'


def test_conflicting_organ_and_missing_fields_remain_hold():
    assert m.classify_botanical_description({'corolla colour':'white','petal colour':'violet'})[0]=='HOLD'
    assert m.classify_botanical_description({'corolla colour':'pale'})[0]=='HOLD'
    assert m.classify_botanical_description({'petal colour':'blue','flower colour':'violet'})[0]=='blue_purple'


def test_nhm_duplicate_conflicts_and_trinomials_removed_without_fuzzy_join():
    d=pd.DataFrame({'taxon':['Aster albus','Aster albus','Aster bluus','Aster bluus var. minor','Aster ruber'],
       'petal colour':['white','blue','purple','white',None],
       'corolla colour':[None,None,None,None,'red'],'flower colour':[None]*5})
    out,s=m.extract_botanical(d)
    assert s['ambiguous_conflicted_duplicate_names']==1
    assert out.loc[out.taxon.eq('Aster albus'),'label'].item()=='HOLD'
    assert not out.taxon.str.contains('var.').any()
    assert out.loc[out.taxon.eq('Aster bluus'),'label'].item()=='blue_purple'


def test_species_means_accept_unknown_photo_classification_and_use_external_label():
    values=[]
    for species,label in [('Aster albus','white'),('Aster bluus','blue_purple')]:
        for i in range(3):
            d=dict(species=species,genus='Aster',label=label,photo_id=species+str(i))
            for k in m.COLS:d[k]=float(i+1)
            d['measurement_status']='not_classifiable_photo'
            values.append(d)
    t=m.species_means(pd.DataFrame(values),np.zeros(8),np.ones(8))
    assert len(t)==2 and t.n_photo.tolist()==[3,3]
    assert sorted(t.white)==[0,1]


def make_ready_source(seed=8,ng=35):
    rng=np.random.default_rng(seed)
    z=[]
    for g in range(ng):
        for sector in (0,3):
            for s in range(4):
                name=f'Genus{g} species{s}{sector}'
                y='white' if s%2 else 'blue_purple'
                row=dict(species=name,genus=f'Genus{g}',label=y,photo_id=name,sector=sector)
                for key in m.COLS:row[key]=float(rng.normal(0,1))
                z.append(row)
    return pd.DataFrame(z)


def test_ready_shared_genus_rainfall_binary_coefficients_sum_zero():
    f=make_ready_source(); f['east']=f.sector.eq(3)
    x=m.shared_genera(f,3)
    assert x.genus.nunique()==35
    left=m.species_means(x.loc[x.east],np.zeros(8),np.ones(8))
    right=m.species_means(x.loc[~x.east],np.zeros(8),np.ones(8))
    a,b=m.regression_matrices(left,True),m.regression_matrices(right,True)
    q=m.paired_bootstrap(a,b,nreps=59,seed=1)
    assert q['status']=='OK'
    assert np.isfinite(q['point'])
    assert q['CI95'][0]<=q['CI95'][1]
    assert abs(float(m.fit_sum(a,sorted(a)).sum()))<1e-10


def test_entire_four_region_family_not_evaluated_on_one_sector():
    f=make_ready_source(ng=29)
    a=m.single_sector(f,np.zeros(8),np.ones(8),3,3,True,39)
    assert a['status']=='HOLD_INDEPENDENT_BOTANICAL_COLOUR_SUPPORT'
    assert a['support']['n_genus']==29
