import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from fcp_pipeline.fcp_environment_profiles import paired_profiles
from fcp_pipeline.fcp_shared_structure import baseline,join_metric,with_environment_heterogeneity,inspect_commonality


def fixtures():
    photo=[]
    env=[]
    for j in range(24):
        for i in range(40):
            taxon=100+j
            pid=f'{taxon}-{i}'
            white=(i%3==0)
            lat=35 + (i%10)*.12
            lon=139 + (i//10)*.12
            photo.append(dict(photo_id=pid,inat_taxon_id=taxon,species=f'Sp{taxon}',
                              morph='white' if white else 'red_pink',
                              global_classifiable=True,observer_id=f'u{i%5}'))
            env.append(dict(photo_id=pid,inat_taxon_id=taxon,latitude=lat,longitude=lon,
                            soil_pH=4.5+(i%10)*.02*(j+1),
                            location_source='direct_source_photo_coordinate'))
    return pd.DataFrame(photo),pd.DataFrame(env)


def test_all_photo_interspecies_W_environment_spread_join(tmp_path):
    photo,env=fixtures()
    spread,sort=paired_profiles(photo,env,grid_degrees=.1)
    assert len(spread)==24 and len(sort)==24
    assert set(spread.variable)=={'soil_pH'}
    assert spread.n_environment_photos.eq(40).all()
    assert sort.n_white.eq(14).all()
    assert sort.n_nonwhite.eq(26).all()
    source=tmp_path/'source.csv'
    pd.DataFrame({'inat_taxon_id':range(100,124),'species':[f'Sp{i}' for i in range(100,124)],
        'n_classifiable':[40]*24,'n_classifiable_observations':[40]*24,
        'n_classifiable_observers':[5]*24,'D_raw':[.4]*24,'D_observer_equal':[.4]*24,
        'measured_status':['multiple_photo_states_detected']*24}).to_csv(source,index=False)
    w=tmp_path/'w.csv'
    pd.DataFrame({'inat_taxon_id':range(100,124),'w_alignment_sq':np.linspace(.3,.95,24),
                  'method':['frozen_study_specific_H2']*24}).to_csv(w,index=False)
    spread_file=tmp_path/'spread.csv';spread.to_csv(spread_file,index=False)
    s=join_metric(baseline(source),w,'W')
    s,het=with_environment_heterogeneity(s,spread_file)
    report,_=inspect_commonality(s,environmental_spread=het)
    results=report['W_vs_withinspecies_environmental_heterogeneity']
    assert len(results)==1 and results[0]['n_species']==24
    assert results[0]['outcome']=='W_white_chromatic_axis'
    assert report['W_environment_frozen_construction_null_completed'] is False


def test_no_nearest_mois_photo_proxy(tmp_path):
    photo,env=fixtures()
    env['location_source']='mois_nearest_photo_5km'
    with pytest.raises(ValueError,match='Nearest mois photo'):
        paired_profiles(photo,env)
    _,env=fixtures()
    z=env.copy()
    z['photo_id']=z.photo_id.apply(lambda x:'changed-'+x)
    with pytest.raises(ValueError,match='does not cover every supplied photo'):
        paired_profiles(photo,z)
