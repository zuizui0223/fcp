import importlib.util
from pathlib import Path
import pandas as pd
import pytest
f=Path(__file__).resolve().parents[1]/'scripts/analysis/freeze_fcp_u100_intraspecific_roster_20261014.py'
spec=importlib.util.spec_from_file_location('u100',f);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def fixture():
 frame=pd.DataFrame([(i,f'Species {i}') for i in range(1,10)],columns=m.ID)
 pool=frame[frame.inat_taxon_id.gt(2)].copy();pool['after_observer_cap']=130
 p500=pool[pool.inat_taxon_id.isin([3,4])].copy()
 third=pool[pool.inat_taxon_id.isin([5,6])][m.ID].copy()
 disc=frame[frame.inat_taxon_id.eq(1)].copy();val=frame[frame.inat_taxon_id.eq(2)].copy()
 return frame,pool,p500,third,disc,val

def test_identity_partition():
 r,c=m.roster_from_identities(*fixture());assert len(r)==9
 assert list(c.values())==[2,2,2,3]
 assert set(r.loc[r.fresh_outcome_unexposed_candidate,'inat_taxon_id'])=={7,8,9}
 assert r.loc[r.historical_membership.eq('legacy_discovery_or_validation'),'after_observer_cap'].isna().all()

def test_invariant_outcome_column_ignored():
 fs=list(fixture());fs[0]['colour']='pink';fs[1]['W']=0.999
 r,_=m.roster_from_identities(*fs)
 assert 'colour' not in r.columns and 'W' not in r.columns

def test_cap_only_historical_parent_no_invented_legacy():
 r,_=m.roster_from_identities(*fixture())
 assert r.loc[r.inat_taxon_id.lt(3),'capacity_provenance'].str.contains('attested_parent').all()
 assert r.loc[r.inat_taxon_id.ge(3),'after_observer_cap'].eq(130).all()

def test_wrong_mapped_taxon_name_holds():
 f,p,q,t,a,b=fixture();t.loc[t.inat_taxon_id.eq(5),'species']='Incorrect new botanical name'
 with pytest.raises(ValueError,match='identity'):m.roster_from_identities(f,p,q,t,a,b)

def test_prior_P500_and_third_overlap_holds():
 f,p,q,t,a,b=fixture();t.loc[t.inat_taxon_id.eq(6),m.ID]=[4,'Species 4']
 with pytest.raises(ValueError,match='overlap'):m.roster_from_identities(f,p,q,t,a,b)

def test_prior_discovery_and_validation_overlap_holds():
 f,p,q,t,a,b=fixture();b.loc[b.inat_taxon_id.eq(2),m.ID]=[1,'Species 1']
 with pytest.raises(ValueError,match='overlap'):m.roster_from_identities(f,p,q,t,a,b)

def test_cap_threshold_cannot_be_relaxed():
 f,p,q,t,a,b=fixture();p.loc[p.inat_taxon_id.eq(8),'after_observer_cap']=99
 with pytest.raises(ValueError,match='observer-cap'):m.roster_from_identities(f,p,q,t,a,b)

def test_parent_observer_cap_conflict():
 f,p,q,t,a,b=fixture();q.loc[q.inat_taxon_id.eq(3),'after_observer_cap']=100
 with pytest.raises(ValueError,match='capacity diverged'):m.roster_from_identities(f,p,q,t,a,b)

def test_historical_count_no_rescue():
 with pytest.raises(ValueError,match='attrition'):m.roster_from_identities(*fixture(),expected=m.EXPECTED)
