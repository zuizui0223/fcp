"""TDD guards; deliberately initially red before any production implementation."""
import csv, math, sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from field_spacetime import load_rows,eligible_years,compute_variance,cv_predictions,run


def artificial(vals):
    return [{'year':y,'transect':tr,'side':side,'blue':n-w,'white':w,'n':n,'p':w/n,'section':tr+'_'+side}
           for (y,tr,side,n,w) in vals]


def test_missing_and_zero_counts_are_preserved_not_white():
    d=load_rows(ROOT/'source/linanthus_counts.csv')
    assert len(d)==28
    assert sum(x['n'] for x in d)==100500
    assert all(x['p'] is None for x in d if x['n']==0)
    assert sum(x['year']==1997 for x in d)==4


def test_support_year_gate_only_counts_and_two_level_threshold():
    d=load_rows(ROOT/'source/linanthus_counts.csv')
    assert eligible_years(d,20)==[1992,1993,1995,1996,1998]
    assert eligible_years(d,100)==[1992,1993,1995,1998]
    changed=[{**x,'white':0,'blue':x['n'],'p':(0. if x['n'] else None)} for x in d]
    assert eligible_years(changed,20)==eligible_years(d,20)


def test_variance_exact_identity_on_synthetic_balanced_design():
    vals=[(y,tr,side,100,w) for y in (1,2,3) for tr in ('north','south') for side in ('east','west') for w in [80 if side=='east' else 20]]
    d=artificial(vals); v=compute_variance(d,[1,2,3])
    assert v['spatial_fraction'] == pytest.approx(1)
    assert v['temporal_variance'] == pytest.approx(0)
    assert v['interaction_variance'] == pytest.approx(0)
    assert v['total_variance'] == pytest.approx(v['spatial_variance']+v['temporal_variance']+v['interaction_variance'])


def test_loyo_never_uses_heldout_year_as_training():
    d=load_rows(ROOT/'source/linanthus_counts.csv');years=eligible_years(d,20)
    pred=cv_predictions(d,years)
    assert len(pred)==4*len(years)*3
    original={(x['year'],x['section'],x['model']):x['prediction'] for x in pred}
    changed=[]
    for x in d:
        if x['year']==1995 and x['n']:
            changed.append({**x,'white':0,'blue':x['n'],'p':0.})
        else: changed.append(x.copy())
    again=cv_predictions(changed,years)
    for x in again:
        if x['year']==1995:
            assert x['prediction']==pytest.approx(original[x['year'],x['section'],x['model']])


def test_longevity_1994_1997_not_scored():
    z=run(ROOT/'source/linanthus_counts.csv')
    assert [y['year'] for y in z['excluded_years']]==[1994,1997]
    assert z['excluded_years'][1]['n_plant_records']==0


def test_loyo_correctly_predicts_constant_spatial_cline():
    d=artificial([(y,tr,side,100,80 if side=='east' else 20) for y in (1,2,3) for tr in ('north','south') for side in ('east','west')]);p=cv_predictions(d,[1,2,3])
    assert all(x['error_sq']<1e-15 for x in p if x['model'] in ('SIDE','SECTION'))
    assert all(x['error_sq']>0 for x in p if x['model']=='GLOBAL')


def test_whole_observed_dataset_has_invariant_raw_source_ledger():
    z=run(ROOT/'source/linanthus_counts.csv')
    assert z['n_counted_plant_records_all_years']==100500
    assert z['source_sha256']=='fa408c7ad6e521a8410252d903b024a057b1ee4f7e27406927e77fcc22bf4ac7'
    assert z['n_years_main']==5 and z['n_section_years_main']==20


def test_row_order_does_not_change_results():
    d=load_rows(ROOT/'source/linanthus_counts.csv')
    p=cv_predictions(d,eligible_years(d,20))
    q=cv_predictions(list(reversed(d)),eligible_years(d,20))
    assert p==q


def test_source_percent_mismatch_is_flagged():
    z=run(ROOT/'source/linanthus_counts.csv')
    mismatch=[d for d in z['printed_percent_mismatches'] if d['year']==1992 and d['section']=='south_west']
    assert len(mismatch)==1
    assert mismatch[0]['actual_white_pct']==pytest.approx(100*26/1383)


def test_year_prediction_ledger_has_three_models_and_4_sections():
    z=run(ROOT/'source/linanthus_counts.csv')
    assert len(z['primary']['cv_predictions'])==5*4*3
    assert all(sum(1 for x in z['primary']['cv_predictions'] if x['year']==y)==12 for y in [1992,1993,1995,1996,1998])
    assert z['primary']['out_of_year_global_mse']>z['primary']['out_of_year_side_mse']
