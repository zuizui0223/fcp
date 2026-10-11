from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/analysis'))
import fcp_joint_environment_taxongeo_20261011 as m


def sample_raw():
    # Taxon2 has a usable-looking western photo AND an unclassified eastern photo.
    return pd.DataFrame({'inat_taxon_id':[1,2,2,3,4,5,6],
       'photo_id':list(range(7)), 'latitude':[0.,0.,0.,0.,0.,np.nan,0.],
       'longitude':[100.,-30.,100.5,89.,-50.,-60.,85.],
       'morph':['white','white','unknown','blue_purple','white','white','white']})

def test_longitude_sectors_cover_boundaries():
    x=m.sector_ids(np.array([-180.,-90.01,-90.,-.01,0.,89.99,90.,180.]))
    assert np.array_equal(x,[0,0,1,1,2,2,3,3])
    with pytest.raises(ValueError):m.sector_ids(np.array([181.]))

def test_full_species_footprint_including_unknown_colours_is_purged():
    raw=sample_raw(); masks,stats=m.spatial_exclusions(raw,buffer_km=500)
    assert {1,2,5}.issubset(masks[3]['sector_purged'])
    assert 4 not in masks[3]['sector_purged']
    assert masks[3]['sector_purged']<=masks[3]['buffer500']

def test_buffer_distance_is_geodesic_and_not_longitude_difference():
    raw=pd.DataFrame({'inat_taxon_id':[1,2,3,4], 'photo_id':range(4),
       'latitude':[0.,0.,0.,0.], 'longitude':[90.1,89.9,80.,-179.9]})
    masks,stats=m.spatial_exclusions(raw,buffer_km=500)
    assert 2 not in masks[3]['sector_purged']
    assert 2 in masks[3]['buffer500']
    assert 3 not in masks[3]['buffer500']

def test_antimeridian_buffer_keeps_close_cross_boundary_photos_apart():
    raw=pd.DataFrame({'inat_taxon_id':[1,2,3], 'photo_id':range(3),
       'latitude':[0.,0.,0.], 'longitude':[179.9,-179.9,-90.]})
    masks,_=m.spatial_exclusions(raw,buffer_km=500)
    assert 2 in masks[3]['buffer500']
    assert 3 not in masks[3]['buffer500']

def test_source_location_values_not_colour_determine_exclusions():
    d=sample_raw();first,_=m.spatial_exclusions(d)
    d['morph']='wrong';d['photo_id']=np.arange(len(d))+100
    second,_=m.spatial_exclusions(d)
    assert first==second

def test_training_selector_excludes_whole_genus_fold_and_forbidden_species():
    sp=pd.DataFrame({'inat_taxon_id':[1,2,3,4],'genus':['A','A','B','C'],'fold':[0,0,1,2]})
    v=m.select_train(sp,0,{4})
    assert v.inat_taxon_id.tolist()==[3]

def test_repeated_sector_species_average_not_photographic_abundance():
    d=pd.DataFrame({'inat_taxon_id':[1,1,2],'genus':['A','A','B'],'sector':[0,1,0],
        'base_total':[1.,3.,10.],'full_total':[0.,0.,9.],
        'base_white_axis':[.5,1.5,5.],'full_white_axis':[0.,0.,4.5],
        'base_remaining_hue':[.5,1.5,5.],'full_remaining_hue':[0.,0.,4.5]})
    d=pd.concat([d.assign(inat_taxon_id=d.inat_taxon_id+10*i,genus=d.genus+str(i)) for i in range(5)],ignore_index=True)
    result=m.summarize_rows(d,reps=99)
    assert result['n_species']==10
    assert np.isclose(result['species_equal']['total']['gain'],1.5)
    assert np.isclose(result['genus_equal']['total']['gain'],1.5)

def test_feature_support_identical_training_and_far_point():
    rng=np.random.default_rng(3);train=rng.normal(size=(500,4))
    supported,dist,threshold=m.environmental_support(train,np.vstack((train.mean(0),[100.]*4)))
    assert supported.tolist()==[True,False]
    assert dist[0]<threshold<dist[1]

def test_feature_support_fits_only_training_distribution():
    rng=np.random.default_rng(4);x=rng.normal(size=(100,4));one=np.ones((1,4))
    _,a,t1=m.environmental_support(x,one)
    _,b,t2=m.environmental_support(x,np.vstack((one,[1e6]*4)))
    assert t1==t2 and np.isclose(a[0],b[0])

def test_model_does_not_read_target_hue_labels():
    rng=np.random.default_rng(6);tr=pd.DataFrame(rng.normal(size=(100,2)),columns=['a','b'])
    tr['genus']=['A']*100
    for i,h in enumerate(m.CLASSES):tr['hue_'+h]=.25
    te=pd.DataFrame({'a':[0.,1.],'b':[0.,1.],'genus':['B','C']})
    one=m.predict(tr,te,('a','b'),'absolute')
    for h in m.CLASSES:te['hue_'+h]=1e9
    two=m.predict(tr,te,('a','b'),'absolute')
    assert np.array_equal(one,two)

def test_full_environment_and_climate_feature_counts():
    assert len(m.FULL)==36 and len(m.CLIMATE)==26
    assert set(m.CLIMATE)<set(m.FULL)
