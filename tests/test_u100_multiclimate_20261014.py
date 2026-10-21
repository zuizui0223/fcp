"""Synthetic QA only: never treat these as observed ecological measurements."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

MODULE=Path(__file__).resolve().parents[1]/'scripts'/'analysis'/'u100_multiclimate_20261014.py'
spec=importlib.util.spec_from_file_location('multiclimate',MODULE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_frozen_predictors_and_multiplicity():
    assert len(m.PREDICTORS)==10
    assert m.PREDICTORS[:5]==('bio1','bio4','bio5','bio6','bio7')
    assert m.PREDICTORS[-1]=='srad_mean'
    p=np.array([.001,.02,.03,np.nan]+[1.]*6)
    v=m.holm_adjust(p)
    assert v[0]==pytest.approx(.01)
    assert v[1]>=.05
    assert v[3]==1.
    assert (v>=0).all() and (v<=1).all()


def test_equal_species_and_observer_local_cells():
    rows=[]
    for sp in ['spA','spB']:
        for ob,k in [('o1',0),('o2',1),('o3',2)]:
            for j in range(12):
                w=j<6
                rows.append({'species':sp,'observer_id':ob,'latitude':21.1+k,'longitude':111.1,
                     'white':w,'z_bio5':(1 if w else -1)+(0.1*k),
                     'near_clip_fraction':0.2+(0.04 if w else 0.)})
    p=pd.DataFrame(rows)
    # All nine predictors can be equal in the synthetic counterfactual.
    for x in m.PREDICTORS:
        p['z_'+x]=p.z_bio5
    tab,summary,contract=m.contrasts(p)
    assert len(tab)==2*10
    assert len(summary)==10
    assert tab.n_paired_observers.min()==3
    assert tab.n_mixed_grid_cells.min()==3
    assert np.allclose(tab.white_minus_nonwhite_sd,2)
    assert np.allclose(tab.observer_balanced_delta,2)
    assert np.allclose(tab.paired_observer_delta,2)
    assert np.allclose(tab.site_cell_delta,2)
    assert not contract['confirmatory_evidence']


def test_observer_or_site_no_overlap_does_not_impute_zero():
    p=[]
    for j in range(12):
        for w in (True,False):
            p.append({'species':'sp','observer_id':str(w),'latitude':0 if w else 24,
                      'longitude':100,'white':w,'near_clip_fraction':0.01,
                      **{'z_'+v:float(w) for v in m.PREDICTORS}})
    tab,_,_=m.contrasts(pd.DataFrame(p))
    assert pd.isna(tab.paired_observer_delta).all()
    assert pd.isna(tab.site_cell_delta).all()
    assert (tab.n_paired_observers==0).all()


def test_frozen_source_check_cannot_be_silently_changed():
    arr=pd.DataFrame([{'species':f'S{i}','variable':var,'white_minus_nonwhite_sd':i/100}
                      for var in ('bio5','bio14','srad_mean') for i in range(281)])
    import json
    ref={'results':[{'variable':var,'n_species_delta':281,'mean_delta_white_minus_nonwhite_SD':float(arr.query('variable==@var').white_minus_nonwhite_sd.mean()),
                  'median_delta_white_minus_nonwhite_SD':float(arr.query('variable==@var').white_minus_nonwhite_sd.median())}
                    for var in ('bio5','bio14','srad_mean')]}
    panel=pd.DataFrame({'species':[f'S{i%281}' for i in range(11000)]})
    assert len(m.frozen_gate(arr,ref,panel))==3
    arr.loc[0,'white_minus_nonwhite_sd']=22
    with pytest.raises(ValueError,match='frozen'):
        m.frozen_gate(arr,ref,panel)


def test_raster_and_frozen_quality_join(tmp_path):
    root=tmp_path;bd=root/'bio';sd=root/'srad';bd.mkdir();sd.mkdir()
    transform=from_origin(-180,90,1,1)
    a=np.tile(np.linspace(1,20,360,dtype='float32'),(180,1))
    for index in m.BIO_KEYS:
        with rasterio.open(bd/f'wc2.1_10m_bio_{index}.tif','w',driver='GTiff',width=360,height=180,
            count=1,dtype='float32',crs='EPSG:4326',transform=transform,nodata=-9999) as dst:
            dst.write(a+index,1)
    for index in range(1,13):
        with rasterio.open(sd/f'wc2.1_10m_srad_{index:02}.tif','w',driver='GTiff',width=360,height=180,
            count=1,dtype='float32',crs='EPSG:4326',transform=transform,nodata=-9999) as dst:
            dst.write(a+index*2,1)
    photos=[];technical=[];join=[]
    for sp in ('spA','spB'):
        for n in range(40):
            pid=f'{sp}_{n}'
            photos.append({'photo_id':pid,'species':sp,'observer_id':f'ob{n%4}',
                'morph':'white' if n<20 else 'red_pink','global_classifiable':True,
                'latitude':20.3,'longitude':100.0+(n%10)})
            technical.append({'measurement_id':pid,'near_clip_fraction':0.02})
            join.append({'measurement_id':pid,'photo_id':pid})
    measured=root/'measured.csv';pd.DataFrame(photos).to_csv(measured,index=False)
    tech=root/'technical.csv.gz';pd.DataFrame(technical).to_csv(tech,index=False,compression='gzip')
    jo=root/'join.csv';pd.DataFrame(join).to_csv(jo,index=False)
    high=root/'high.csv';pd.DataFrame({'measurement_id':['spA_0']}).to_csv(high,index=False)
    panel,meta=m.make_panel(measured,tech,jo,high,bd,sd,verify_sha=False)
    assert meta['eligible_species']==2
    assert len(panel)==79
    assert 'z_bio17' in panel
    assert np.isfinite(panel['z_srad_mean']).all()
    p,summary,meta2=m.contrasts(panel)
    assert len(summary)==10
    assert len(p)==20
    assert int(summary.n.max())==2
