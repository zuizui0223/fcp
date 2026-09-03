#!/usr/bin/env python3
"""Freeze the 12,000-photo geometry-balanced preimage reservoir without opening pixels."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path

SALT="jbi-random-photo-atlas-v1"; SEED=20260903; NLON=24; NLAT=12; TARGET=12000

def h(*x: object) -> str:
    return hashlib.sha256((SALT+'\x1f'+'\x1f'.join(map(str,x))).encode()).hexdigest()

def cell(lat: float, lon: float) -> tuple[int,int,str]:
    j=min(NLON-1,max(0,int(math.floor(((lon+180.0)%360.0)/360.0*NLON))))
    u=(math.sin(math.radians(lat))+1.0)/2.0
    i=min(NLAT-1,max(0,int(math.floor(u*NLAT))))
    return i,j,f"E{i:02d}_{j:02d}"

def read(path: Path) -> list[dict[str,str]]:
    with path.open(encoding='utf-8',newline='') as f: return list(csv.DictReader(f))

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('--metadata',type=Path,required=True); ap.add_argument('--metadata-qc',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--qc-out',type=Path,required=True); a=ap.parse_args()
    meta=read(a.metadata); mq=json.loads(a.metadata_qc.read_text())
    a.out.parent.mkdir(parents=True,exist_ok=True)
    status='pass_preimage_reservoir'; reason=''
    selected=[]
    if mq.get('status')!='complete' or not mq.get('metadata_complete'):
        status='not_evaluable_metadata_incomplete'; reason='closed-window metadata universe did not complete'
    elif len(meta)<TARGET:
        status='not_evaluable_insufficient_closed_window_metadata'; reason=f'eligible={len(meta)} < frozen_target={TARGET}'
    else:
        by={}
        for r in meta:
            lat=float(r['latitude']); lon=float(r['longitude']); ii,jj,cid=cell(lat,lon)
            q=dict(r); q['grid_lat_index']=ii; q['grid_lon_index']=jj; q['grid_cell_id']=cid
            q['preimage_priority']=h(SEED,'row',r['observation_id'],r['photo_id'])
            by.setdefault(cid,[]).append(q)
        for cid in by: by[cid].sort(key=lambda r:(r['preimage_priority'],int(r['observation_id'])))
        cells=sorted(by,key=lambda c:h(SEED,'cell',c)); curs={c:0 for c in cells}
        while len(selected)<TARGET:
            progress=False
            for cid in cells:
                k=curs[cid]
                if k<len(by[cid]):
                    selected.append(by[cid][k]); curs[cid]=k+1; progress=True
                    if len(selected)>=TARGET: break
            if not progress: break
        if len(selected)!=TARGET:
            status='not_evaluable_reservoir_selection'; reason=f'selected={len(selected)}'
    fields=[]
    if selected:
        for idx,r in enumerate(selected):
            r['reservoir_rank']=idx+1
            r['measurement_id']=f"random-atlas-v1-{int(r['photo_id']):012d}"
            r['species_blind_id']='sp-'+hashlib.sha256((SALT+'\x1f'+r['species_key']).encode()).hexdigest()[:16]
            r['image_filename']=f"{r['measurement_id']}.jpg"
        fields=list(selected[0])
    else:
        fields=['measurement_id','species_blind_id','image_filename','observation_id','photo_id','species_key','latitude','longitude','photo_url']
    with a.out.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n'); w.writeheader(); w.writerows(selected)
    counts={}
    for r in selected: counts[r['grid_cell_id']]=counts.get(r['grid_cell_id'],0)+1
    qc={
      'status':status,'reason':reason,'protocol':'jbi-random-photo-atlas-preimage-reservoir-v1','seed':SEED,'salt':SALT,
      'grid':{'longitude_bins':NLON,'sin_latitude_bins':NLAT,'equal_area':True},'frozen_target':TARGET,
      'eligible_metadata':len(meta),'selected':len(selected),'occupied_cells_metadata':len({cell(float(r['latitude']),float(r['longitude']))[2] for r in meta}) if meta else 0,
      'occupied_cells_selected':len(counts),'min_selected_per_occupied_cell':min(counts.values()) if counts else 0,'max_selected_per_cell':max(counts.values()) if counts else 0,
      'unique_species_keys_selected':len({r['species_key'] for r in selected}),'pixel_accessed':False,
      'selection_rule':'deterministic hash-ordered round robin across frozen 24x12 lon-by-sin(lat) cells; no colour or species eligibility used'
    }
    a.qc_out.write_text(json.dumps(qc,indent=2,sort_keys=True)+'\n'); print(json.dumps(qc,indent=2,sort_keys=True)); return 0
if __name__=='__main__': raise SystemExit(main())
