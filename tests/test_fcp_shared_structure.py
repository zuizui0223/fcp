"""Synthetic-only comparative ecological contrasts."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from fcp_pipeline.fcp_shared_structure import baseline,join_metric,with_environment,inspect_commonality


def prep(tmp_path):
    n=50
    path=tmp_path/'photo_FCP.csv'
    pd.DataFrame({'inat_taxon_id':range(100,100+n),'species':[f'Species {i}' for i in range(n)],
        'n_classifiable':[60]*40+[10]*10,'n_classifiable_observations':[55]*40+[10]*10,
        'n_classifiable_observers':[7]*40+[2]*10,
        'D_raw':list(np.linspace(.1,.5,40))+[np.nan]*10,
        'D_observer_equal':list(np.linspace(.08,.44,40))+[np.nan]*10,
        'measured_status':['multiple_photo_states_detected']*40+['insufficient_classifiable']*10}).to_csv(path,index=False)
    return path


def test_structural_comparison_is_species_equal_and_no_inferred_null(tmp_path):
    src=baseline(prep(tmp_path))
    w=tmp_path/'w.csv';sp=tmp_path/'sp.csv'
    pd.DataFrame({'inat_taxon_id':range(100,135),'w_alignment_sq':np.linspace(.2,.6,35),'method':['frozen_study_specific_H2']*35}).to_csv(w,index=False)
    pd.DataFrame({'inat_taxon_id':range(100,140),'local_pairs_50km':[32]*40,
                  'D_pair':np.linspace(.1,.4,40),'depletion_50km':np.linspace(-.02,.08,40)}).to_csv(sp,index=False)
    src=join_metric(src,w,'W');src=join_metric(src,sp,'spatial')
    receipt,_=inspect_commonality(src)
    assert receipt['n_D_evaluable_species']==40
    assert receipt['n_source_species']==50
    assert receipt['species_equal_descriptive_summaries']['within_species_W_direction']['n_species']==35
    assert receipt['frozen_primary_H2_W_null_recomputed'] is False
    assert receipt['frozen_primary_50km_null_recomputed'] is False
    assert receipt['cross_species_structural_relations']['D_raw_vs_depletion_50km']['n_species']==40
    assert receipt['no_colour_prevalence_or_local_genetic_FCP_claim']


def test_photo_coordinate_required_for_environment_main(tmp_path):
    src=baseline(prep(tmp_path))
    env=tmp_path/'env.csv'
    df=pd.DataFrame({'inat_taxon_id':range(100,140),'variable':['soil_pH']*40,
           'white_minus_nonwhite_sd':np.linspace(-1,1,40),'n_white':[6]*40,'n_nonwhite':[6]*40,
           'location_source':['nearest_mois_photo_5km']*40})
    df.to_csv(env,index=False)
    with pytest.raises(ValueError,match='Proxy or nearest-photo'):
        with_environment(src,env)
    df['location_source']='direct_source_photo_coordinate'
    df.to_csv(env,index=False)
    source,effects=with_environment(src,env)
    receipt,_=inspect_commonality(source,effects)
    assert receipt['ecological_environment_associations'][0]['n_species']==40
    assert not receipt['environment_confirmatory_support']


def test_fake_h2_method_and_spatial_missing_pair_guard(tmp_path):
    src=baseline(prep(tmp_path))
    w=tmp_path/'w.csv'
    pd.DataFrame({'inat_taxon_id':[100],'w_alignment_sq':[.6],'method':['generic_disttrait_refit']}).to_csv(w,index=False)
    with pytest.raises(ValueError,match='identical frozen study'):
        join_metric(src,w,'W')
    s=tmp_path/'s.csv'
    pd.DataFrame({'inat_taxon_id':[100],'local_pairs_50km':[2], 'depletion_50km':[.1], 'D_pair':[.2]}).to_csv(s,index=False)
    with pytest.raises(ValueError,match='below frozen'):
        join_metric(src,s,'spatial')
