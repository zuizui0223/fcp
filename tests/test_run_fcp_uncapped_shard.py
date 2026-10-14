"""No network calls. Provenance/cover-all checks use synthetic provider pages."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import pandas as pd
import pytest

PATH=Path(__file__).resolve().parents[1]/'scripts/acquisition/run_fcp_uncapped_shard.py'
spec=importlib.util.spec_from_file_location('fcp_uncapped_shard',PATH)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def frame(tmp_path):
    p=tmp_path/'source.csv'
    pd.DataFrame({'inat_taxon_id':[11,12,13,14],
                  'species':['Species alba','Species rubra','Species caerulea','Species lutea'],
                  'after_observer_cap':[0,10,100,200]}).to_csv(p,index=False)
    return p


def observation(n):
    return {'id':n,'taxon':{'id':11,'ancestor_ids':[11,12,13,14]},'user':{'id':42},
            'geojson':{'coordinates':[135.2,35.6]},'observed_on':'2023-04-23',
            'photos':[{'id':1000+n*3+i,'url':f'https://example.org/{n}-{i}/medium.jpg',
                       'license_code':'cc-by','attribution':'synthetic'} for i in range(3)]}


def test_all_taxa_are_preserved_and_resumed_past_page_budget(tmp_path):
    p=frame(tmp_path);calls=[]
    def request(params):
        calls.append((params['taxon_id'],params.get('id_above',0)))
        start=params.get('id_above',0)
        if start==0:
            return {'results':[observation(i) for i in range(1,201)]}
        if start==200:
            return {'results':[observation(i) for i in range(201,203)]}
        return {'results':[]}
    out=tmp_path/'out'
    a=m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-14',
                  pages_per_species=1,species_budget=2,request=request,expected_frame_species=4)
    assert a['status']=='in_progress' and a['partial_species']==2
    assert a['not_started_species']==2 and a['completed_species']==0
    assert not a['all_available_photo_ids_confirmed']
    b=m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-14',
                  pages_per_species=1,request=request,expected_frame_species=4)
    assert b['status']=='in_progress' and b['completed_species']==2
    assert b['partial_species']==2
    assert not b['all_available_photo_ids_confirmed']
    c=m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-14',
                  pages_per_species=1,request=request,expected_frame_species=4)
    assert c['status']=='complete' and c['new_species_worked_this_invocation']==2
    assert c['photo_manifest_rows_including_retry_duplicates']==4*202*3
    assert c['all_available_photo_ids_confirmed']
    d=m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-14',
                  pages_per_species=1,request=request,expected_frame_species=4)
    assert d['status']=='complete' and d['new_species_worked_this_invocation']==0
    assert (out/'shard_000_of_001'/'source_contract.json').exists()
    assert calls.count((11,0))==1
    assert calls.count((11,200))==1
    with pytest.raises(ValueError,match='source/date/geometry'):
        m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-13',
                    pages_per_species=1,request=request,expected_frame_species=4)


def test_provider_error_stops_and_receipt_remains_partial(tmp_path):
    p=frame(tmp_path);out=tmp_path/'out'
    def request(_):raise RuntimeError('429 rate limit')
    with pytest.raises(RuntimeError,match='checkpoint'):
        m.run_shard(p,out,shard=0,nshards=1,created_before='2026-10-14',
                    pages_per_species=1,request=request,expected_frame_species=4)
    report=json.loads((out/'shard_000_of_001'/'shard_progress.json').read_text())
    assert report['status']=='in_progress'
    assert report['last_source_error']['type']=='RuntimeError'
    assert not report['all_available_photo_ids_confirmed']


def test_partition_is_stable_and_never_sampling_cap(tmp_path):
    p=frame(tmp_path);out=tmp_path/'out'
    def request(_):return {'results':[]}
    rs=[]
    for shard in range(3):
        r=m.run_shard(p,out,shard=shard,nshards=3,created_before='2026-10-14',
                      request=request,expected_frame_species=4)
        rs.append(r)
    assert sum(x['n_source_species_in_shard'] for x in rs)==4
    assert sum(x['completed_species'] for x in rs)==4
    assert all(x['photo_cap'] is None for x in rs)
    assert all(x['status']=='complete' for x in rs)