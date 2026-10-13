import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from robbed_morph import group_plant_units, site_contrasts, permutation_main, date_overlap_sensitivity, verify_source


def fixture():
    return pd.DataFrame([
        dict(COUNTRY='Sweden',Year=2024,**{'GPS-code':'G_1'},Plant_ID='1',Flower_morph='violet',Robbed_flower='yes',Length_flower='40',Date='1-Jul-24'),
        dict(COUNTRY='Sweden',Year=2024,**{'GPS-code':'G_1'},Plant_ID='1',Flower_morph='violet',Robbed_flower='no',Length_flower='42',Date='1-Jul-24'),
        dict(COUNTRY='Sweden',Year=2024,**{'GPS-code':'G_1'},Plant_ID='2',Flower_morph='pink',Robbed_flower='no',Length_flower='44',Date='1-Jul-24'),
        dict(COUNTRY='Sweden',Year=2024,**{'GPS-code':'G_1'},Plant_ID='3',Flower_morph='pink',Robbed_flower='eaten',Length_flower='43',Date='1-Jul-24'),
        dict(COUNTRY='Sweden',Year=2024,**{'GPS-code':'G_1'},Plant_ID='4',Flower_morph='violet',Robbed_flower='no',Length_flower='39',Date='2-Jul-24'),
    ])


def test_repeated_flowers_collapsed_one_plant_with_correct_fraction():
    d=group_plant_units(fixture())
    assert len(d)==4
    x=d[d.Plant_ID=='1'].iloc[0]
    assert x['n_flower_rows']==2 and x['n_yes']==1 and x['n_no']==1
    assert x['robbed_fraction']==0.5 and x['length_plant_mean']==41


def test_unknown_robbery_never_coded_as_no():
    d=group_plant_units(fixture())
    x=d[d.Plant_ID=='3'].iloc[0]
    assert x['n_yes']==0 and x['n_no']==0
    assert np.isnan(x['robbed_fraction'])


def test_conflicting_morph_holds_complete_plant():
    f=fixture()
    row=f.iloc[0].copy();row['Flower_morph']='pink';f=pd.concat([f,pd.DataFrame([row])],ignore_index=True)
    x=group_plant_units(f); assert x[x.Plant_ID=='1'].iloc[0]['exclusion_reason']=='CONFLICTING_MORPH'


def test_no_empty_site_comparison_masquerades_as_zero():
    x=group_plant_units(fixture())
    with pytest.raises(ValueError, match='HOLD_INSUFFICIENT'):
        site_contrasts(x,['G_1'],min_each=5)


def test_permutation_preserves_units_and_is_deterministic():
    x=group_plant_units(fixture())
    a=permutation_main(x,['G_1'],reps=19,seed=1,min_each=1)
    b=permutation_main(x,['G_1'],reps=19,seed=1,min_each=1)
    assert a==b
    assert 1/20 <= a['p_two_sided']<=1


def test_mixed_dates_excluded_from_date_strata():
    f=fixture();extra=f.iloc[0].copy();extra['Date']='3-Jul-24'
    f=pd.concat([f,pd.DataFrame([extra])],ignore_index=True)
    x=group_plant_units(f)
    assert x[x.Plant_ID=='1'].iloc[0]['n_dates']==2
    a=date_overlap_sensitivity(x,['G_1'],min_sites=1,min_plants=1,reps=19,seed=2)
    assert a['status']=='HOLD_INSUFFICIENT_DATE_OVERLAP'


def test_missing_plant_id_rejected_not_accidentally_joined():
    f=fixture();f.loc[0,'Plant_ID']=np.nan
    with pytest.raises(ValueError, match='missing plant key'):
        group_plant_units(f)


def test_source_bytes_match_publisher_exact_digest(tmp_path):
    f=tmp_path/'x';f.write_bytes(b'csv bytes')
    with pytest.raises(ValueError, match='source MD5'):
        verify_source(f)


def source_path():
    return ROOT/'source'/'digitalis_field_individuals.csv'


def test_real_file_is_authenticated_1323_flower_rows_and_site_exclusion():
    from robbed_morph import load_source, PRIMARY_SITES
    r,s=load_source(source_path())
    assert len(r)==1323 and s['encoding']=='cp1252'
    assert len(r[(r.Year==2024)&r['GPS-code'].isin(PRIMARY_SITES)])==977
    assert len(r[(r.Year==2024)&r['GPS-code'].eq('H_2')])==153


def test_source_2024_sample_has_226_cluster_units_and_records_conflicts():
    from robbed_morph import load_source, PRIMARY_SITES
    source,_=load_source(source_path())
    plants=group_plant_units(source)
    q=plants[(plants.Year==2024)&(plants['GPS-code'].isin(PRIMARY_SITES))]
    assert len(q)==226
    assert (q.exclusion_reason=='CONFLICTING_MORPH').sum()==2
    assert (q.exclusion_reason=='MISSING_MORPH').sum()==1
    assert (q.exclusion_reason=='NO_ASSESSED_FLOWERS').sum()==9


def test_primary_site_equal_contrast_excludes_selected_H2_and_2023():
    from robbed_morph import load_source, PRIMARY_SITES
    source,_=load_source(source_path())
    plants=group_plant_units(source)
    rows=site_contrasts(plants,PRIMARY_SITES,min_each=5)
    assert len(rows)==5
    assert sum(x['n_pink_groups']+x['n_violet_groups'] for x in rows)==181
    assert all(x['site']!='H_2' for x in rows)
    assert any(x['site']=='G_1' and x['violet_mean']==x['pink_mean']==0 for x in rows)


def test_main_real_robbery_inference_not_significant_with_frozen_seed():
    from robbed_morph import load_source, PRIMARY_SITES
    source,_=load_source(source_path())
    plants=group_plant_units(source)
    result=permutation_main(plants,PRIMARY_SITES,reps=9999,seed=2026101301,min_each=5)
    assert abs(result['site_equal_estimate']-0.04705668169298781)<1e-14
    assert result['p_two_sided']==0.2647


def test_descriptive_variance_partition_reassembles_and_no_causal_claim():
    from robbed_morph import descriptive_variance_partition
    f=pd.DataFrame([
        {'Year':2024,'GPS-code':'X','morph':'violet','exclusion_reason':'','robbed_fraction':0.0},
        {'Year':2024,'GPS-code':'X','morph':'pink','exclusion_reason':'','robbed_fraction':0.0},
        {'Year':2024,'GPS-code':'Y','morph':'violet','exclusion_reason':'','robbed_fraction':1.0},
        {'Year':2024,'GPS-code':'Y','morph':'pink','exclusion_reason':'','robbed_fraction':0.5},
    ])
    got=descriptive_variance_partition(f,sites=['X','Y'])
    assert got['n_plant_patch_units']==4
    assert got['site_component_fraction']>got['morph_given_site_component_fraction']
    assert sum(got[k] for k in ('site_component_fraction','morph_given_site_component_fraction','within_site_morph_residual_fraction'))==pytest.approx(1)
