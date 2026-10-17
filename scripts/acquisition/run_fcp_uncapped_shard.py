#!/usr/bin/env python3
"""Checkpointed, outcome-blind iNaturalist full-photo census for one taxon shard.

Run again with identical frame/date/shard to resume. A per-run page budget is a
compute/rate-limit checkpoint, NEVER a total photo or per-species sampling cap.
A shard is COMPLETE only if every species checkpoint has status=complete.

This script does NOT download images, run segmentation, or compute FCP outcomes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from typing import Callable

from fcp_pipeline.fcp_allphotos import opportunity_frame, opportunity_schedule, sha256_file


def collector_module():
    path=Path(__file__).with_name('collect_fcp_uncapped_metadata.py')
    spec=importlib.util.spec_from_file_location('fcp_source_metadata_collector',path)
    if spec is None or spec.loader is None:
        raise ValueError('Original uncapped metadata collector unavailable')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _atomic_json(path:Path, obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.part')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    os.replace(tmp,path)


def _species_state(directory:Path,taxon_id:int,species:str,created_before:str)->dict:
    checkpoint=directory/'checkpoint.json'
    if not checkpoint.exists():
        return {'taxon_id':taxon_id,'species':species,'status':'not_started',
                'photo_manifest_rows':0,'batches':0,'cursor':0}
    meta=json.loads(checkpoint.read_text(encoding='utf-8'))
    for k,v in [('taxon_id',taxon_id),('species',species),('created_d2',created_before)]:
        if meta.get(k)!=v:
            raise ValueError(f'Metadata checkpoint mismatch for {taxon_id}: {k}')
    if meta.get('photo_cap') is not None or meta.get('colour_pixels_opened') is not False:
        raise ValueError(f'Metadata checkpoint violates no-cap/blind contract: {taxon_id}')
    if meta.get('status') not in ('in_progress','complete'):
        raise ValueError(f'Metadata checkpoint has invalid state for {taxon_id}')
    manifest=directory/'photo_manifest.ndjson'
    if not manifest.exists():
        raise ValueError(f'Missing checkpoint photo manifest for {taxon_id}')
    # Photo-row count is provenance, not an inclusion or stopping criterion.
    with manifest.open('rb') as fh:
        rows=sum(bool(b.strip()) for b in fh)
    return {'taxon_id':taxon_id,'species':species,'status':meta['status'],
            'photo_manifest_rows':rows,'batches':int(meta.get('batches',0)),
            'cursor':int(meta.get('cursor',0)),
            'photo_manifest_sha256':sha256_file(manifest)}


def run_shard(frame_path:Path, output_root:Path, *, shard:int, nshards:int,
              created_before:str, pages_per_species:int=10,
              species_budget:int|None=None,
              request:Callable[[dict],dict]|None=None,
              expected_frame_species:int=42111)->dict:
    """Advance selected shard without dropping unvisited or partial taxa.

    species_budget controls *invocations*, not biological species selection:
    all source taxa are included in the receipt even when not yet visited.
    """
    from datetime import date
    try:
        parsed=date.fromisoformat(created_before)
        if parsed.isoformat()!=created_before:
            raise ValueError('bad date normalization')
    except ValueError as exc:
        raise ValueError('created_before must be frozen YYYY-MM-DD') from exc
    if nshards<=0 or not 0<=shard<nshards:
        raise ValueError('Invalid shard or nshards')
    if pages_per_species<1 or (species_budget is not None and species_budget<1):
        raise ValueError('A work-unit must have a positive page/species budget')
    frame=opportunity_frame(frame_path,expected_frame_species)
    plan=opportunity_schedule(frame,nshards)
    selected=plan[plan['partition'].eq(shard)]
    # Stable selection sorted by canonical source taxon ID, not apparent phenotype.
    identities=selected[['inat_taxon_id','species']].itertuples(index=False,name=None)
    ids=[(int(t),str(s)) for t,s in identities]
    root=output_root/f'shard_{shard:03}_of_{nshards:03}'
    root.mkdir(parents=True,exist_ok=True)
    contract=root/'source_contract.json'
    source={'schema':'fcp_uncapped_shard_source_contract_v1',
            'source_species_frame_sha256':sha256_file(frame_path),
            'expected_source_species':expected_frame_species,
            'n_source_species_in_shard':len(ids),
            'shard':shard,'nshards':nshards,'created_d2':created_before,
            'photo_cap':None,'outcome_blind':True}
    if contract.exists():
        if json.loads(contract.read_text(encoding='utf-8'))!=source:
            raise ValueError('Shard source/date/geometry changed after freeze')
    else:
        _atomic_json(contract,source)
    collector=collector_module()
    fetch=request
    if fetch is None:
        fetch=lambda params:collector.network_fetch(params,pause_seconds=1.2)
    progressed=0
    fail=None
    for taxon_id,species in ids:
        directory=root/'species'/str(taxon_id)
        state=_species_state(directory,taxon_id,species,created_before)
        if state['status']=='complete':
            continue
        if species_budget is not None and progressed>=species_budget:
            break
        progressed+=1
        try:
            collector.collect(taxon_id,species,created_before,directory,fetch,
                              max_batches=pages_per_species)
        except Exception as exc:
            # Source service errors do not turn a partial ledger into a complete one.
            fail={'taxon_id':taxon_id,'type':type(exc).__name__,
                  'message':str(exc)[:400]}
            break  # fail-closed so repeated HTTP 429 does not create a request storm
    states=[_species_state(root/'species'/str(t),t,s,created_before) for t,s in ids]
    ncomplete=sum(x['status']=='complete' for x in states)
    npartial=sum(x['status']=='in_progress' for x in states)
    receipt={**source,'status':'complete' if ncomplete==len(ids) else 'in_progress',
             'completed_species':ncomplete,'partial_species':npartial,
             'not_started_species':len(ids)-ncomplete-npartial,
             'new_species_worked_this_invocation':progressed,
             'photo_manifest_rows_including_retry_duplicates':sum(x['photo_manifest_rows'] for x in states),
             'species':states,'last_source_error':fail,
             'pixels_opened':False,'FCP_outcomes_measured':False,
             'all_available_photo_ids_confirmed':ncomplete==len(ids) and fail is None}
    _atomic_json(root/'shard_progress.json',receipt)
    if fail is not None:
        raise RuntimeError(f'Source acquisition stopped without losing checkpoint: {fail}')
    return receipt


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frame',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--shard',type=int,required=True)
    p.add_argument('--shards',type=int,default=128)
    p.add_argument('--created-before',type=str,required=True)
    p.add_argument('--pages-per-species-per-run',type=int,default=10)
    p.add_argument('--species-per-run',type=int,default=None)
    p.add_argument('--expected-source-species',type=int,default=42111)
    a=p.parse_args(argv)
    result=run_shard(a.frame,a.output_root,shard=a.shard,nshards=a.shards,
                     created_before=a.created_before,
                     pages_per_species=a.pages_per_species_per_run,
                     species_budget=a.species_per_run,
                     expected_frame_species=a.expected_source_species)
    print(json.dumps({k:result[k] for k in ('status','completed_species','partial_species',
                         'not_started_species','photo_manifest_rows_including_retry_duplicates')},
                     ensure_ascii=False))
    return 0


if __name__=='__main__':
    raise SystemExit(main())