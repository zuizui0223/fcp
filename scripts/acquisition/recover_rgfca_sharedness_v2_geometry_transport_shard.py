#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from typing import Mapping, Sequence
import pandas as pd
from fcp_pipeline.random_photo_h9_pool import h9_query_for_species, observer_cap, parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient
from scripts.acquisition.freeze_rgfca_sharedness_v2_geometry_shard import selected_species, load_exclusions, row_hash
ROOT=Path(__file__).resolve().parents[2]
C=ROOT/'docs/supporting/rgfca_sharedness_v2_geometry_transport_recovery_contract_v1.json'
OLD=ROOT/'docs/supporting/rgfca_sharedness_v2_geometry_freeze_contract_v1.json'

def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--shard-index',type=int,required=True); ap.add_argument('--shard-count',type=int,required=True); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args()
 c=json.loads(C.read_text()); old=json.loads(OLD.read_text()); q=c['query']
 if c['status']!='frozen_after_random-batch-geometry-failure_before_recovery_requests' or not c['scientific_design_unchanged']: raise RuntimeError('recovery contract drift')
 if c['selected_species']!=old['activation']['selected_species'] or c['photos_per_species']!=old['fixed_photo_frame']['photos_per_species']: raise RuntimeError('scientific frame changed')
 s=selected_species(); shard=s.loc[(s.species_order%a.shard_count)==a.shard_index].copy().reset_index(drop=True)
 oldq=old['metadata_sampling']; allowed=frozenset(x.casefold() for x in q['allowed_photo_licenses']); exobs,exphoto=load_exclusions([ROOT/p for p in oldq['prior_experiment_exclusion_sources']])
 client=InaturalistObservationClient(request_interval_seconds=q['request_interval_seconds'],timeout_seconds=45,max_retries=q['request_retries'],user_agent=f'fcp-rgfca-v2-geometry-recovery-{a.shard_index}/1.0')
 outs=[]; audits=[]
 for sr in shard.itertuples(index=False):
  parsed=[]; so=set(); sp=set(); err=''; exhausted=False; pages=0
  for page in range(1,q['maximum_pages_per_species']+1):
   p=h9_query_for_species(int(sr.inat_taxon_id),per_page=q['per_page'],maximum_positional_accuracy_m=q['maximum_positional_accuracy_m'],flowering_term_id=q['flowering_term_id'],flowering_term_value_id=q['flowering_term_value_id'],allowed_photo_licenses=tuple(q['allowed_photo_licenses']))
   p['order_by']=q['stable_order_by']; p['order']=q['stable_order']; p['page']=page; pages+=1
   try: payload=client.observations(p)
   except Exception as e: err=f'{type(e).__name__}:{str(e)[:180]}'; break
   raw=payload.get('results') or []
   for ob in raw:
    if not isinstance(ob,Mapping): continue
    x=parse_h9_observation(ob,expected_taxon_id=int(sr.inat_taxon_id),maximum_positional_accuracy_m=q['maximum_positional_accuracy_m'],allowed_photo_licenses=allowed)
    if x is None: continue
    oid,pid=int(x['observation_id']),int(x['photo_id'])
    if oid in exobs or pid in exphoto or oid in so or pid in sp: continue
    so.add(oid); sp.add(pid); parsed.append(x)
   fr=pd.DataFrame(parsed); cap=observer_cap(fr,q['observer_cap']) if len(fr) else fr
   if len(cap)>=q['stop_when_after_observer_cap_reaches']: break
   total=payload.get('total_results')
   if len(raw)<q['per_page'] or (total is not None and page*q['per_page']>=int(total)): exhausted=True; break
  fr=pd.DataFrame(parsed); cap=observer_cap(fr,q['observer_cap']) if len(fr) else fr
  if len(cap):
   cap=cap.copy();
   if 'row_hash' not in cap: cap['row_hash']=[row_hash(o,p) for o,p in zip(cap.observation_id,cap.photo_id)]
   sel=cap.sort_values('row_hash',kind='mergesort').head(c['photos_per_species']).copy().reset_index(drop=True); sel['sample_role']=sr.sample_role; sel['species_order']=int(sr.species_order); sel['photo_order']=range(1,len(sel)+1)
   fields=old['allowed_saved_fields']; sel=sel[fields]; outs.append(sel)
  else: sel=cap
  audits.append({'species_order':int(sr.species_order),'sample_role':sr.sample_role,'species':sr.species,'inat_taxon_id':int(sr.inat_taxon_id),'pages_requested':pages,'after_observer_cap':int(len(cap)),'retained':int(len(sel)),'full_fixed_300':bool(len(sel)==c['photos_per_species']),'api_exhausted':exhausted,'terminal_request_error':err})
 out=pd.concat(outs,ignore_index=True) if outs else pd.DataFrame(columns=old['allowed_saved_fields']); aud=pd.DataFrame(audits).sort_values('species_order')
 a.output_dir.mkdir(parents=True,exist_ok=True); op=a.output_dir/f'rgfca_v2_geometry_recovery_photos_{a.shard_index:02d}.csv'; oa=a.output_dir/f'rgfca_v2_geometry_recovery_audit_{a.shard_index:02d}.csv'; om=a.output_dir/f'rgfca_v2_geometry_recovery_{a.shard_index:02d}.json'
 out.to_csv(op,index=False); aud.to_csv(oa,index=False); om.write_text(json.dumps({'protocol':c['protocol'],'shard_index':a.shard_index,'shard_count':a.shard_count,'full_species':int(aud.full_fixed_300.sum()),'retained_photos':int(len(out)),'terminal_errors':int(aud.terminal_request_error.fillna('').astype(str).str.len().gt(0).sum()),'image_pixels_opened':False,'flower_colour_opened':False,'lineage':{'contract_sha256':sha(C),'photos_sha256':sha(op),'audit_sha256':sha(oa)}},indent=2)+'\n')
if __name__=='__main__': main()
