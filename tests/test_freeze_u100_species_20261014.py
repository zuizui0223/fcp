"""Synthetic identity and outcome-blind input contract tests."""
import importlib.util
from pathlib import Path
import pandas as pd
import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/analysis/freeze_fcp_u100_intraspecific_roster_20261014.py'
spec=importlib.util.spec_from_file_location('u100_freeze',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def fixtures():
    taxa=[(i,f'Species {i}',130) for i in range(1,10)]
    frame=pd.DataFrame(taxa,columns=mod.FIELDS)
    p100=frame[frame.inat_taxon_id>2].copy()
    p500=frame[frame.inat_taxon_id.isin([3,4])].copy()
    selected=frame[frame.inat_taxon_id.isin([5,6])][mod.FIELDS[:2]].copy()
    return frame,p100,p500,selected


def test_full_roster_yields_disjoint_original_cohorts():
    x=fixtures()
    roster,counts=mod.make_roster(*x)
    assert len(roster)==9
    assert [counts[k] for k in mod.STATUS]==[2,2,2,3]
    assert set(roster.loc[roster.fresh_outcome_unexposed_candidate,'inat_taxon_id'])=={7,8,9}
    assert roster.inat_taxon_id.is_unique


def test_identifiers_only_even_if_outcomes_accidentally_supplied():
    f,p,q,t=fixtures()
    f['photo_hue']='white'; p['W']=1234; t['local_environment']='fake'
    r,_=mod.make_roster(f,p,q,t)
    assert 'photo_hue' not in r.columns and 'W' not in r.columns
    assert 'local_environment' not in r.columns


def test_colour_and_environment_columns_cannot_define_source_capacity():
    f,p,q,t=fixtures()
    f=f.drop(columns=['after_observer_cap']); f['colour_count']=1000
    with pytest.raises(ValueError,match='after_observer_cap'):
        mod.make_roster(f,p,q,t)


def test_names_must_agree_for_same_taxon_id():
    f,p,q,t=fixtures();t.loc[t.inat_taxon_id.eq(5),'species']='Other incorrect ID'
    with pytest.raises(ValueError,match='identity'):
        mod.make_roster(f,p,q,t)


def test_p500_and_third_cohort_nonoverlap():
    f,p,q,t=fixtures();t.loc[t.inat_taxon_id.eq(6),['inat_taxon_id','species']]=[4,'Species 4']
    with pytest.raises(ValueError,match='overlap'):
        mod.make_roster(f,p,q,t)


def test_old_p100_capacity_drift_fails_closed():
    f,p,q,t=fixtures();p.loc[p.inat_taxon_id.eq(8),'after_observer_cap']=999
    with pytest.raises(ValueError,match='capacity'):
        mod.make_roster(f,p,q,t)


def test_unexpected_cohort_size_fails_closed():
    f,p,q,t=fixtures()
    with pytest.raises(ValueError,match='attrition'):
        mod.make_roster(f,p,q,t,expect_counts=mod.STATUS)


def test_low_capacity_is_not_u100():
    f,p,q,t=fixtures();f.loc[f.inat_taxon_id.eq(1),'after_observer_cap']=99
    r,counts=mod.make_roster(f,p,q,t)
    assert len(r)==8 and counts['legacy_discovery_or_validation']==1
