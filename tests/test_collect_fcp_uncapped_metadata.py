import importlib.util
import json
from pathlib import Path
import pytest

SCRIPT = Path(__file__).resolve().parents[1]/'scripts/acquisition/collect_fcp_uncapped_metadata.py'
spec=importlib.util.spec_from_file_location('uncapped_metadata_collector',SCRIPT)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def obj(index, obs_taxon=44, photo_count=1):
    return {'id':index,'taxon':{'id':obs_taxon,'ancestor_ids':[44]},
            'user':{'id':index%7+1},'geojson':{'coordinates':[135.0,35.0]},
            'observed_on':'2025-04-03','photos':[
              {'id':index*10+j,'url':f'https://example.org/photos/{index}/{j}/medium.jpg',
               'license_code':'cc-by','attribution':'creator'} for j in range(photo_count)]}


def test_cursor_iterates_without_photo_limit(tmp_path):
    calls=[]
    def request(params):
        calls.append(params.copy())
        cur=params.get('id_above',0)
        if cur==0:return {'results':[obj(i,photo_count=2) for i in range(1,201)]}
        if cur==200:return {'results':[obj(i,photo_count=2) for i in range(201,204)]}
        return {'results':[]}
    done=module.collect(44,'Species alba','2026-10-14',tmp_path,request)
    assert done['status']=='complete'
    assert done['photo_cap'] is None
    assert done['cursor']==203
    assert calls[1]['id_above']==200
    rows=[json.loads(x) for x in (tmp_path/'photo_manifest.ndjson').read_text().splitlines()]
    assert len(rows)==406  # all two licensed photos per iNat observation
    assert rows[0]['photo_url_large'].endswith('/large.jpg')
    assert len({x['photo_id'] for x in rows})==406


def test_incomplete_checkpoint_is_not_full_census(tmp_path):
    def request(params):
        return {'results':[obj(i) for i in range(params.get('id_above',0)+1,params.get('id_above',0)+201)]}
    a=module.collect(44,'Species alba','2026-10-14',tmp_path,request,max_batches=1)
    assert a['status']=='in_progress'
    assert a['photo_cap'] is None
    assert a['cursor']==200
    b=module.collect(44,'Species alba','2026-10-14',tmp_path,request,max_batches=1)
    assert b['status']=='in_progress' and b['cursor']==400
    with pytest.raises(ValueError,match='current frozen species/date'):
        module.collect(45,'Other','2026-10-14',tmp_path,request,max_batches=1)


def test_bad_cursor_and_rights_filtered(tmp_path):
    def bad(params):return {'results':[obj(5),obj(4)]}
    with pytest.raises(ValueError,match='not strictly monotonic'):
        module.collect(44,'Species alba','2026-10-14',tmp_path,bad)
    q=obj(100,photo_count=2);q['photos'][1]['license_code']=None
    assert len(module.observation_rows(q,species='S',taxon_id=44,freeze='2026-10-14',query_label='test'))==1
    q['taxon']={'id':99,'ancestor_ids':[]}
    assert not module.observation_rows(q,species='S',taxon_id=44,freeze='2026-10-14',query_label='test')
