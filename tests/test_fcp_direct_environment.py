"""Synthetic raster checks; not biological evidence."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from fcp_pipeline.fcp_direct_environment import (
    MOIS_36, MOIS_SOIL, read_manifest, sample_rows, stream_all_photos,
)


def raster(path:Path,values):
    with rasterio.open(path,'w',driver='GTiff',width=2,height=2,count=1,dtype='float32',
        crs='EPSG:4326',transform=from_origin(0,2,1,1),nodata=-9999) as ds:
        ds.write(np.asarray(values,dtype='float32'),1)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(tmp_path):
    f=tmp_path/'soil.tif'
    sha=raster(f,[[10,-9999],[20,30]])
    obj={'schema':'fcp_mois36_direct_photo_raster_manifest_v1',
         'features':{'soil_pH':{'unit':'source_pH_scale','source':'synthetic test raster',
                                'reducer':'identity','rasters':[{'path':'soil.tif','sha256':sha}]}}}
    path=tmp_path/'m.json';path.write_text(json.dumps(obj))
    return path


def test_frozen_mois36_names_and_direct_raster_sampling(tmp_path):
    assert len(MOIS_36)==36 and len(MOIS_SOIL)==10
    name=manifest(tmp_path);meta=read_manifest(name)
    dat=pd.DataFrame({'inat_taxon_id':[11]*4,'photo_id':['1','2','3','4'],
                      'latitude':[1.5,.5,1.5,44.], 'longitude':[.5,.5,1.5,139.],
                      'positional_accuracy_m':[100,1500,100,25]})
    out=sample_rows(dat,meta)
    assert out.soil_pH.iloc[0]==10
    assert out.soil_pH.iloc[1]==20
    assert np.isnan(out.soil_pH.iloc[2])  # nodata, not soil = 0
    assert np.isnan(out.soil_pH.iloc[3])  # out of bounds, not soil = 0
    assert out['location_source'].eq('direct_source_photo_coordinate').all()
    assert out.accuracy_le1000m.tolist()==[True,False,True,True]
    with pytest.raises(ValueError,match='All 36'):
        read_manifest(name,require_all36=True)


def test_full_metadata_photo_identity_and_sha_guard(tmp_path):
    name=manifest(tmp_path)
    measured=tmp_path/'measure.csv'
    pd.DataFrame({'inat_taxon_id':[11,11],'photo_id':['3','4'],'latitude':[0.5,1.5],
                  'longitude':[0.5,0.5], 'positional_accuracy_m':[20,20]}).to_csv(measured,index=False)
    dest=tmp_path/'out.csv'
    receipt=stream_all_photos([measured],name,dest)
    assert receipt['n_photo_ids']==2
    assert receipt['per_feature_available_photo_count']['soil_pH']==2
    assert receipt['is_full_mois36'] is False
    assert pd.read_csv(dest).soil_pH.tolist()==[20,10]
    with pytest.raises(ValueError,match='overwrite'):
        stream_all_photos([measured],name,dest)
    obj=json.loads(name.read_text())
    obj['features']['soil_pH']['rasters'][0]['sha256']='0'*64
    name.write_text(json.dumps(obj))
    with pytest.raises(ValueError,match='SHA256 mismatch'):
        read_manifest(name)
