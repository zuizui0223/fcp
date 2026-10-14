"""Synthetic only: official Open Data snapshot-to-fresh-FCP source contract tests."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition.import_inat_open_data_snapshot import (
    ingest, lineage_source_taxon, normalize_license, SNAPSHOT_SOURCE,
)
from fcp_pipeline.fcp_allphotos import opportunity_schedule, opportunity_frame, sha256_file
from fcp_pipeline.fcp_fresh_study import required_source_dirs


def tsv_gz(path:Path, cols:list[str], rows:list[list[str]]) -> Path:
    with gzip.open(path,"wt",encoding="utf-8",newline="") as f:
        f.write("\t".join(cols)+"\n")
        for row in rows:
            f.write("\t".join(str(v) for v in row)+"\n")
    return path


def fixture(tmp_path:Path):
    frame=tmp_path/'frame.csv'
    pd.DataFrame({'inat_taxon_id':[11,12,13,14],
                  'species':['Species alba','Species rubra','Species caerulea','Species lutea'],
                  'after_observer_cap':[0,1,2,100]}).to_csv(frame,index=False)
    taxa=tsv_gz(tmp_path/'taxa.csv.gz',['taxon_id','ancestry','rank'],[
        [11,'1/5','species'],[12,'1/5','species'],[13,'1/5','species'],
        [14,'1/5','species'],[111,'1/5/11','subspecies'],[99,'1/5','species']])
    obs=tsv_gz(tmp_path/'observations.csv.gz',
        ['observation_uuid','taxon_id','observer_id','latitude','longitude','positional_accuracy','quality_grade','observed_on'],[
        ['o-1',111,1,35.6,135.2,15,'Research Grade','2023-08-01'],
        ['o-2',12,2,36.0,136.0,200,'research','2023-07-10'],
        ['o-3',13,3,36.2,136.2,10,'Needs ID','2023-06-10'],
        ['o-4',99,4,35,135,25,'Research Grade','2023-03-11'],
        ['o-5',12,5,102,20,20,'research','2023-04-01']])
    photos=tsv_gz(tmp_path/'photos.csv.gz',
        ['photo_id','observation_uuid','observer_id','extension','license','width','height'],[
        [1001,'o-1',1,'jpeg','CC-BY',1024,1024],
        [1002,'o-1',1,'jpg','CC-BY-NC',1024,768],
        [1003,'o-2',2,'png','cc0',640,480],
        [1004,'o-2',2,'jpeg','CC-BY-ND',1024,768],
        [1005,'o-3',3,'jpg','cc-by',1024,768],
        [1006,'o-4',4,'jpg','cc-by',1024,768],
        [1007,'o-5',5,'jpg','cc-by',1024,768]])
    pinned={k:sha256_file(v) for k,v in [('taxa',taxa),('observations',obs),('photos',photos)]}
    return frame,taxa,obs,photos,pinned


def test_source_snapshot_retains_all_photos_and_all_zero_photo_taxa(tmp_path):
    frame,taxa,obs,photos,pinned=fixture(tmp_path)
    out=tmp_path/'source'
    r=ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,expected_species=4,input_sha256=pinned)
    assert r['status']=='COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA'
    assert r['n_source_observations_scanned']==5
    assert r['n_target_taxon_research_grade_observations']==2
    assert r['n_source_photos_scanned']==7
    assert r['n_target_taxon_licensed_photos']==3
    assert r['n_taxa_with_source_photos']==2
    assert not r['image_pixels_opened'] and not r['old_photo_measurements_read']
    dirs,audit=required_source_dirs(frame,out,shards=2,expected=4)
    assert len(dirs)==4 and audit['status']=='SOURCE_CENSUS_COMPLETE'
    plan=opportunity_schedule(opportunity_frame(frame,4),2).set_index('inat_taxon_id')
    p=out/f"shard_{int(plan.loc[11,'partition']):03}_of_002"/'species'/'11'/'photo_manifest.ndjson'
    values=[json.loads(line) for line in p.read_text().splitlines()]
    assert len(values)==2 and {x['photo_id'] for x in values}=={'1001','1002'}
    assert all(x['source_snapshot_date']=='2026-10-14' for x in values)
    assert values[0]['observation_id']=='o-1'
    assert values[0]['photo_url_large'].endswith('/1001/large.jpeg')
    assert all(x['source_query']==SNAPSHOT_SOURCE for x in values)
    species13=out/f"shard_{int(plan.loc[13,'partition']):03}_of_002"/'species'/'13'/'photo_manifest.ndjson'
    assert species13.exists() and species13.stat().st_size==0
    source_contract=json.loads((p.parents[2]/'source_contract.json').read_text())
    assert source_contract['created_d2_field_semantics']=='immutable_source_snapshot_date_not_API_created_d2_query'


def test_original_frame_and_all_input_byte_hashes_are_mandatory(tmp_path):
    frame,taxa,obs,photos,pinned=fixture(tmp_path)
    with pytest.raises(ValueError,match='pre-pinned SHA256'):
        ingest(frame,taxa,obs,photos,tmp_path/'bad','2026-10-14',shards=2,
               expected_species=4,input_sha256={**pinned,'photos':'0'*64})
    changed=pd.read_csv(frame);changed.loc[0,'species']='Different species';changed.to_csv(frame,index=False)
    # The census may only start from the exact externally frozen source frame.
    assert sha256_file(frame)!='5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc'


def test_lineage_and_licenses_are_metadata_only():
    source={11,12,13}
    assert lineage_source_taxon(111,'1/5/11',source)==11
    assert lineage_source_taxon(12,'1\\5',source)==12
    assert lineage_source_taxon(333,'1/5/444',source) is None
    assert normalize_license('CC-BY-NC')=='cc-by-nc'
    assert normalize_license('CC-BY-ND') is None


def test_bulk_snapshot_cannot_mix_with_api_or_altered_photo_manifest(tmp_path):
    frame,taxa,obs,photos,pinned=fixture(tmp_path)
    out=tmp_path/'bulk'
    ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,expected_species=4,input_sha256=pinned)
    contract=next(out.glob('shard_000_*/source_contract.json'))
    original=contract.read_text()
    changed=json.loads(original)
    changed['source_kind']='inat_api_observations_cursor'
    contract.write_text(json.dumps(changed))
    with pytest.raises(ValueError,match='Mixed API and monthly Open Data'):
        required_source_dirs(frame,out,shards=2,expected=4)
    contract.write_text(original)
    plan=opportunity_schedule(opportunity_frame(frame,4),2).set_index('inat_taxon_id')
    m=out/f"shard_{int(plan.loc[11,'partition']):03}_of_002"/'species'/'11'/'photo_manifest.ndjson'
    with m.open('a') as f:f.write('{"unexpected":true}\n')
    with pytest.raises(ValueError,match='snapshot manifest identity/hash mismatch'):
        required_source_dirs(frame,out,shards=2,expected=4)
