"""Synthetic only: official Open Data snapshot-to-fresh-FCP source contract tests."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition.import_inat_open_data_snapshot import (
    ingest, lineage_source_taxon, normalize_license, SNAPSHOT_SOURCE, verify_source_archive_pin,
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
    observers=tsv_gz(tmp_path/'observers.csv.gz',
        ['observer_id','login','name'],[
        [1,'alice_i','Alice Flower'],[2,'bob_i','Bob Botany'],
        [3,'charlie_i','Charlie'],[4,'dora_i','Dora'],[5,'erin_i','Erin']])
    pinned={k:sha256_file(v) for k,v in [('taxa',taxa),('observations',obs),
                                           ('observers',observers),('photos',photos)]}
    return frame,taxa,obs,photos,observers,pinned


def test_source_snapshot_retains_all_photos_and_all_zero_photo_taxa(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    out=tmp_path/'source'
    r=ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,expected_species=4,input_sha256=pinned,observers=observers)
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
    assert 'Alice Flower' in values[0]['attribution']
    assert 'www.inaturalist.org/users/1' in values[0]['attribution']
    assert 'CC-BY' in values[0]['attribution']
    assert r['n_target_photos_with_named_creator']==3
    assert r['n_target_photos_creator_unresolved']==0
    assert r['n_source_observers_scanned']==5
    assert all(x['source_query']==SNAPSHOT_SOURCE for x in values)
    species13=out/f"shard_{int(plan.loc[13,'partition']):03}_of_002"/'species'/'13'/'photo_manifest.ndjson'
    assert species13.exists() and species13.stat().st_size==0
    source_contract=json.loads((p.parents[2]/'source_contract.json').read_text())
    assert source_contract['created_d2_field_semantics']=='immutable_source_snapshot_date_not_API_created_d2_query'


def test_original_frame_and_all_input_byte_hashes_are_mandatory(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    with pytest.raises(ValueError,match='pre-pinned SHA256'):
        ingest(frame,taxa,obs,photos,tmp_path/'bad','2026-10-14',shards=2,
               expected_species=4,input_sha256={**pinned,'photos':'0'*64},observers=observers)
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

def test_official_observer_metadata_required_before_any_source_manifest(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    with pytest.raises(ValueError,match='observer table required'):
        ingest(frame,taxa,obs,photos,tmp_path/'without_authors',
               '2026-10-14',shards=2,expected_species=4,input_sha256=pinned)
    assert not (tmp_path/'without_authors').exists()



def test_bulk_snapshot_cannot_mix_with_api_or_altered_photo_manifest(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    out=tmp_path/'bulk'
    ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,expected_species=4,input_sha256=pinned,observers=observers)
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


def test_full_snapshot_source_pin_requires_matching_original_dated_archive(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    receipt={
        'schema':'fcp_dated_inaturalist_metadata_bundle_pin_v1',
        'status':'FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED',
        'snapshot_date':'2026-09-27',
        's3_object_key':'metadata/inaturalist-open-data-20260927.tar.gz',
        'metadata_files_sha256':dict(pinned),
        'full_bundle_sha256':'a'*64,
        'source_images_opened':False,
        'historical_photo_colours_read':False,
    }
    pin=tmp_path/'source_archive_pin.json'
    pin.write_text(json.dumps(receipt))
    result=verify_source_archive_pin(pin,'2026-09-27',pinned)
    assert result['source_archive_sha256']=='a'*64
    assert result['source_archive_key'].endswith('20260927.tar.gz')
    receipt['s3_object_key']='metadata/inaturalist-open-data-latest.tar.gz'
    pin.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='Dated immutable snapshot archive key'):
        verify_source_archive_pin(pin,'2026-09-27',pinned)
    receipt['s3_object_key']='metadata/inaturalist-open-data-20260927.tar.gz'
    receipt['metadata_files_sha256']['observers']='b'*64
    pin.write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='attestation missing or disagrees'):
        verify_source_archive_pin(pin,'2026-09-27',pinned)


def test_case_sensitive_s3_extension_is_not_lowercased(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    with gzip.open(photos,'rt',encoding='utf-8') as h:
        original=h.read()
    assert '\t1003\t' not in original  # Original rows are one photo ID per first column.
    assert '1003\to-2\t2\tpng' in original
    with gzip.open(photos,'wt',encoding='utf-8') as h:
        h.write(original.replace('1003\to-2\t2\tpng','1003\to-2\t2\tJPG'))
    pinned['photos']=sha256_file(photos)
    out=tmp_path/'newcase'
    got=ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,
               expected_species=4,input_sha256=pinned,observers=observers)
    plan=opportunity_schedule(opportunity_frame(frame,4),2).set_index('inat_taxon_id')
    loc=out/f"shard_{int(plan.loc[12,'partition']):03}_of_002"/'species'/'12'/'photo_manifest.ndjson'
    row=json.loads(loc.read_text().strip())
    assert row['photo_id']=='1003' and row['photo_url_large'].endswith('/1003/large.JPG')
    assert got['original_S3_extension_case_preserved']


def test_cross_observation_duplicate_photo_id_is_quarantined_without_downsampling(tmp_path):
    frame,taxa,obs,photos,observers,pinned=fixture(tmp_path)
    # Photo #1003 was first credited to observation o-2, later re-used in o-1.
    with gzip.open(photos,'rt',encoding='utf-8') as h:
        original=h.read()
    with gzip.open(photos,'wt',encoding='utf-8') as h:
        h.write(original+'1003\to-1\t1\tJPG\tCC-BY\t1024\t768\n')
    pinned['photos']=sha256_file(photos)
    out=tmp_path/'deduplicated'
    got=ingest(frame,taxa,obs,photos,out,'2026-10-14',shards=2,
               expected_species=4,input_sha256=pinned,observers=observers)
    assert got['n_source_photos_scanned']==8
    assert got['n_target_taxon_licensed_photo_metadata_rows']==4
    assert got['n_target_taxon_licensed_photos']==2
    assert got['n_conflicting_photo_ids_quarantined']==1
    assert got['n_taxa_with_source_photos']==1
    dirs,audit=required_source_dirs(frame,out,shards=2,expected=4)
    assert audit['status']=='SOURCE_CENSUS_COMPLETE'
    manifest_rows=[json.loads(line) for d in dirs for line in
                   (d/'photo_manifest.ndjson').read_text().splitlines()]
    assert {r['photo_id'] for r in manifest_rows}=={'1001','1002'}
    assert len(manifest_rows)==2
