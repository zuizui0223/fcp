#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]; C=ROOT/'docs/supporting/rgfca_sharedness_v2_geometry_transport_recovery_contract_v1.json'
def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--input-root',type=Path,required=True); ap.add_argument('--output-photos',type=Path,required=True); ap.add_argument('--output-audit',type=Path,required=True); ap.add_argument('--output-json',type=Path,required=True); a=ap.parse_args(); c=json.loads(C.read_text())
 ps=[]; au=[]; shards=[]
 for i in range(6):
  p=list(a.input_root.rglob(f'rgfca_v2_geometry_recovery_photos_{i:02d}.csv')); u=list(a.input_root.rglob(f'rgfca_v2_geometry_recovery_audit_{i:02d}.csv')); m=list(a.input_root.rglob(f'rgfca_v2_geometry_recovery_{i:02d}.json'))
  if len(p)!=1 or len(u)!=1 or len(m)!=1: raise RuntimeError(f'missing shard {i}')
  mm=json.loads(m[0].read_text());
  if mm['protocol']!=c['protocol'] or mm['image_pixels_opened'] or mm['flower_colour_opened']: raise RuntimeError('bad shard')
  ps.append(pd.read_csv(p[0])); au.append(pd.read_csv(u[0])); shards.append({'shard':i,'photos':sha(p[0]),'audit':sha(u[0])})
 photos=pd.concat(ps,ignore_index=True).sort_values(['species_order','photo_order'],kind='mergesort').reset_index(drop=True); audit=pd.concat(au,ignore_index=True).sort_values('species_order',kind='mergesort').reset_index(drop=True)
 full=(len(audit)==150 and audit.inat_taxon_id.nunique()==150 and audit.full_fixed_300.astype(bool).all() and len(photos)==45000 and photos.photo_id.nunique()==45000 and photos.observation_id.nunique()==45000 and photos.groupby('inat_taxon_id').size().eq(300).all() and audit.terminal_request_error.fillna('').astype(str).str.len().eq(0).all())
 a.output_photos.parent.mkdir(parents=True,exist_ok=True); a.output_json.parent.mkdir(parents=True,exist_ok=True); photos.to_csv(a.output_photos,index=False); audit.to_csv(a.output_audit,index=False)
 r={'protocol':c['protocol'],'status':'complete_transport_recovered_metadata_geometry_150x300' if full else 'transport_recovery_failed_incomplete_frame','complete_fixed_frame':bool(full),'species':int(len(audit)),'training_species':int(audit.loc[audit.sample_role=='training','inat_taxon_id'].nunique()),'evaluation_species':int(audit.loc[audit.sample_role=='evaluation','inat_taxon_id'].nunique()),'retained_photos':int(len(photos)),'photos_per_species':300,'species_below_target':int((~audit.full_fixed_300.astype(bool)).sum()),'terminal_request_error_species':int(audit.terminal_request_error.fillna('').astype(str).str.len().gt(0).sum()),'image_pixels_opened':False,'flower_colour_opened':False,'synthetic_sharedness_qualification_permitted':bool(full),'image_acquisition_permitted':False,'ecological_claim_changed':False,'lineage':{'contract_sha256':sha(C),'combined_photos_sha256':sha(a.output_photos),'combined_audit_sha256':sha(a.output_audit),'shards':shards},'claim_ceiling':c['claim_ceiling']}
 a.output_json.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n'); print(json.dumps(r,indent=2))
if __name__=='__main__': main()
