"""Synthetic image outcomes ONLY; never count these as biological photographs."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from fcp_pipeline.fcp_allphotos import opportunity_frame, opportunity_schedule, sha256_file
from fcp_pipeline.fcp_fresh_study import (
    SOURCE_MEASUREMENT_BLOB, SOURCE_METHOD_COMMIT, SOURCE_ENCODER_SHA256,
    SOURCE_DECODER_SHA256, SOURCE_DETECTOR_SHA256, BIO, freeze_fresh,
    reassemble_fresh, required_source_dirs,
)


def setup_source(tmp:Path, *, complete:bool=True):
    p=tmp/'frozen_taxa.csv'
    taxa=pd.DataFrame({'inat_taxon_id':[11,12,13,14],
                       'species':['S alba','S lutea','S caerulea','S rosea'],
                       'after_observer_cap':[0,1,3,150]})
    taxa.to_csv(p,index=False)
    root=tmp/'metadata'
    plan=opportunity_schedule(opportunity_frame(p,4),2)
    by_tax={int(r.inat_taxon_id):int(r.partition) for r in plan.itertuples()}
    for part in range(2):
        folder=root/f'shard_{part:03}_of_002'
        folder.mkdir(parents=True)
        n=sum(i==part for i in by_tax.values())
        (folder/'source_contract.json').write_text(json.dumps({
            'source_species_frame_sha256':sha256_file(p),
            'expected_source_species':4,'n_source_species_in_shard':n,
            'shard':part,'nshards':2,'created_d2':'2026-10-14',
            'photo_cap':None,'outcome_blind':True}))
    for tax in [11,12,13,14]:
        folder=root/f'shard_{by_tax[tax]:03}_of_002'/'species'/str(tax)
        folder.mkdir(parents=True)
        (folder/'checkpoint.json').write_text(json.dumps({
            'taxon_id':tax,'species':taxa.loc[taxa.inat_taxon_id.eq(tax),'species'].iloc[0],
            'created_d2':'2026-10-14','photo_cap':None,'colour_pixels_opened':False,
            'status':'complete' if complete or tax!=14 else 'in_progress'}))
        with (folder/'photo_manifest.ndjson').open('w') as f:
            # Photo 12345 might have been used in an earlier project, but the
            # FRESH pipeline DOES NOT import or reuse its previous colour result.
            for i in range(2):
                v={'inat_taxon_id':tax,'species':taxa.loc[taxa.inat_taxon_id.eq(tax),'species'].iloc[0],
                   'observation_id':str(tax*10+i),'photo_id':str(12345 if tax==11 and i==0 else tax*100+i),
                   'photo_url_large':f'https://example.org/{tax}/{i}/large.jpg',
                   'photo_license':'cc-by','attribution':'synthetic',
                   'observer_id':'100','latitude':35,'longitude':135,'positional_accuracy_m':35,
                   'observed_on':'2023-07-01','source_created_d2':'2026-10-14',
                   'source_query':'synthetic'}
                f.write(json.dumps(v)+'\n')
    return p,root


def attestation(tmp:Path)->Path:
    p=tmp/'method.json'
    p.write_text(json.dumps({'source_commit':SOURCE_METHOD_COMMIT,
             'measurement_git_blob':SOURCE_MEASUREMENT_BLOB,
             'detector_sha256':SOURCE_DETECTOR_SHA256,
             'efficient_sam_encoder_sha256':SOURCE_ENCODER_SHA256,
             'efficient_sam_decoder_sha256':SOURCE_DECODER_SHA256,
             'location_blind':True,'no_old_measurement_csv_inputs':True,
             'persisted_image_pixels':False}))
    return p


def terminal(tmp:Path, fw:Path,*,missing:bool=False)->Path:
    target=tmp/'sealed';target.mkdir()
    from fcp_pipeline.fcp_fresh_study import SOURCE_METHOD_COMMIT, SOURCE_MEASUREMENT_BLOB, SOURCE_DETECTOR_SHA256, SOURCE_ENCODER_SHA256, SOURCE_DECODER_SHA256
    from fcp_pipeline.fcp_allphotos import sha256_file
    method={'source_commit':SOURCE_METHOD_COMMIT,'measurement_git_blob':SOURCE_MEASUREMENT_BLOB,
            'detector_sha256':SOURCE_DETECTOR_SHA256,'efficient_sam_encoder_sha256':SOURCE_ENCODER_SHA256,
            'efficient_sam_decoder_sha256':SOURCE_DECODER_SHA256,
            'location_blind':True,'no_old_measurement_csv_inputs':True,'persisted_image_pixels':False}
    all_parts=sorted(fw.glob('worker_partition_*.csv'))
    nonempty=[j for j,f in enumerate(all_parts) if len(list(csv.DictReader(f.open())))>0]
    last=nonempty[-1]
    for part,worker_file in enumerate(all_parts):
        workers=list(csv.DictReader(worker_file.open()))
        if not workers:continue
        filename=target/f'terminal_p{part:03}_b0000000.csv'
        with filename.open('w',newline='') as f:
            fields=['measurement_id','morph','measurement_status',*[f'flower_fraction_{c}' for c in BIO]]
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for j,row in enumerate(workers[:-1] if missing and part==last else workers):
                is_good=(j%2==0)
                values={f'flower_fraction_{c}':float((c=='white') if is_good else 0) for c in BIO}
                w.writerow({'measurement_id':row['measurement_id'],
                            'morph':'white' if is_good else 'mixed_uncertain',
                            'measurement_status':'classified_four_state_morph' if is_good else 'image_acquisition_failed',
                            **values})
        written=max(0,len(workers)-1) if missing and part==last else len(workers)
        receipt={'status':'COMPLETE_WORK_UNIT','partition':part,'batch':0,
                 'compute_batch_size_not_photo_cap':100,'number_batches_for_this_partition':1,
                 'n_source_photo_ids':written,'n_terminal_photo_ids':written,
                 'terminal_csv_sha256':sha256_file(filename),'image_bytes_persisted':False,
                 'reused_historical_measurement_csv':False,'full_study_complete':False,
                 'method':method}
        filename.with_suffix('.json').write_text(json.dumps(receipt))
    return target


def test_fresh_freeze_includes_every_taxon_and_remeasures_even_old_id(tmp_path):
    p,meta=setup_source(tmp_path)
    fw=tmp_path/'fresh_worker_only'
    got=freeze_fresh(p,meta,fw,shards=2,workers=3,expected=4)
    assert got['status']=='FRESH_BLIND_METADATA_FROZEN'
    assert got['all_available_licensed_photo_ids']==8
    assert got['historical_measured_CSV_read'] is False
    assert got['reused_old_colour_labels'] is False
    assert (fw/'metadata_join_key.csv').read_text().count('12345')==1
    assert not any('photo_url_large' in (fw/f'worker_partition_{i:03}.csv').read_text().splitlines()[0] for i in range(3))
    seal=terminal(tmp_path,fw)
    result=reassemble_fresh(fw,seal,attestation(tmp_path),tmp_path/'fresh_output')
    assert result['status']=='COMPLETE_FRESH_MEASURED_PHOTO_CENSUS'
    assert result['n_terminal_measurements']==8
    assert result['n_classifiable']==4
    assert result['historical_measurements_reused'] is False
    measured=pd.read_csv(tmp_path/'fresh_output'/'newly_remeasured_photos.csv')
    assert measured.photo_id.nunique()==8
    assert measured.global_classifiable.sum()==4
    assert measured[~measured.global_classifiable].morph.eq('mixed_uncertain').all()
    assert 12345 in measured.photo_id.to_numpy()
    assert measured['attribution'].eq('synthetic').all()
    assert 'attribution' in (fw/'metadata_join_key.csv').read_text().splitlines()[0]
    assert not any('attribution' in (fw/f'worker_partition_{i:03}.csv').read_text().splitlines()[0] for i in range(3))


def test_incomplete_metadata_cannot_freeze_or_be_rebranded(tmp_path):
    p,meta=setup_source(tmp_path,complete=False)
    dirs,status=required_source_dirs(p,meta,shards=2,expected=4)
    assert status['status']=='INCOMPLETE_METADATA_NOT_AUTHORIZED'
    assert status['incomplete_species']==1
    assert len(dirs)==3
    with pytest.raises(ValueError,match='Cannot open any pixels'):
        freeze_fresh(p,meta,tmp_path/'no',shards=2,expected=4)
    assert not (tmp_path/'no'/'worker_partition_000.csv').exists()


def test_old_method_drifts_and_missing_outcomes_fail_closed(tmp_path):
    p,meta=setup_source(tmp_path)
    fw=tmp_path/'firewall'
    freeze_fresh(p,meta,fw,shards=2,expected=4)
    sealed=terminal(tmp_path,fw,missing=True)
    with pytest.raises(ValueError,match='Incomplete fresh ROI-v4'):
        reassemble_fresh(fw,sealed,attestation(tmp_path),tmp_path/'partial')
    assert not (tmp_path/'partial'/'newly_remeasured_photos.csv').exists()
    assert not (tmp_path/'partial'/'fresh_complete_receipt.json').exists()
    with next(sealed.glob('terminal_*.csv')).open('w',newline='') as f:
        f.write('measurement_id,morph,measurement_status\n')
    with pytest.raises(ValueError,match='photo checksum'):
        reassemble_fresh(fw,sealed,attestation(tmp_path),tmp_path/'wrong_checksum')
    corrupted=next(sealed.glob('terminal_*.csv'))
    corrupted_receipt=corrupted.with_suffix('.json')
    meta=json.loads(corrupted_receipt.read_text())
    from fcp_pipeline.fcp_allphotos import sha256_file
    meta['terminal_csv_sha256']=sha256_file(corrupted)
    corrupted_receipt.write_text(json.dumps(meta))
    with pytest.raises(ValueError,match='nine palette'):
        reassemble_fresh(fw,sealed,attestation(tmp_path),tmp_path/'wrong_columns')
    method=attestation(tmp_path)
    a=json.loads(method.read_text());a['detector_sha256']='someone_else';method.write_text(json.dumps(a))
    with pytest.raises(ValueError,match='detector_sha256'):
        reassemble_fresh(fw,sealed,method,tmp_path/'wrong_model')


def test_no_colour_based_species_eligibility_and_refuse_source_swap(tmp_path):
    p,meta=setup_source(tmp_path)
    # Frame changes cause absolute contract failure, NOT optimistic partial run.
    altered=pd.read_csv(p)
    altered.loc[0,'species']='S changed'
    altered.to_csv(p,index=False)
    with pytest.raises(ValueError,match='Source'):
        freeze_fresh(p,meta,tmp_path/'different',shards=2,expected=4)