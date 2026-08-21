#!/usr/bin/env python3
"""Download deterministic coordinate samples for focal and matched-control species."""
from __future__ import annotations
import argparse,csv,json,math,time,urllib.parse,urllib.request
from pathlib import Path

API='https://api.gbif.org/v1/occurrence/search'
FIELDS=['canonical_name','family','role','focal_species','match_level','gbif_key','decimalLatitude','decimalLongitude','year','basisOfRecord']

def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as h:return list(csv.DictReader(h))

def request(url,timeout,retries):
    err=None
    for i in range(retries):
        try:
            req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'fcp-niche/1.0'})
            with urllib.request.urlopen(req,timeout=timeout) as r:return json.load(r)
        except Exception as e:
            err=e; time.sleep(2**i)
    raise RuntimeError(str(err))

def species_table(focals,controls):
    out={}
    for r in rows(focals):
        n=r.get('canonical_name','').strip()
        if n: out[(n,'focal',n)]={'canonical_name':n,'family':r.get('family',''),'role':'focal','focal_species':n,'match_level':'self'}
    if controls:
        for r in rows(controls):
            n=r.get('control_candidate','').strip(); f=r.get('focal_species','').strip()
            if n: out[(n,'control',f)]={'canonical_name':n,'family':r.get('control_family',''),'role':'control','focal_species':f,'match_level':r.get('match_level','')}
    return list(out.values())

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--focals',required=True); ap.add_argument('--controls')
    ap.add_argument('--out',required=True); ap.add_argument('--qc-out',required=True); ap.add_argument('--limit',type=int,default=300)
    ap.add_argument('--timeout',type=int,default=45); ap.add_argument('--retries',type=int,default=4); ap.add_argument('--delay',type=float,default=.2)
    a=ap.parse_args(); spp=species_table(a.focals,a.controls); output=[]; errors=[]
    for i,s in enumerate(spp):
        params={'scientificName':s['canonical_name'],'hasCoordinate':'true','occurrenceStatus':'present','limit':min(300,max(1,a.limit))}
        try: data=request(API+'?'+urllib.parse.urlencode(params),a.timeout,a.retries)
        except RuntimeError as e: errors.append({'canonical_name':s['canonical_name'],'error':str(e)[:300]}); continue
        seen=set()
        for x in data.get('results',[]):
            try: lat=float(x.get('decimalLatitude')); lon=float(x.get('decimalLongitude'))
            except (TypeError,ValueError): continue
            if not(math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180): continue
            if abs(lat)<1e-9 and abs(lon)<1e-9: continue
            cell=(round(lat,3),round(lon,3))
            if cell in seen: continue
            seen.add(cell)
            output.append({**s,'gbif_key':x.get('key',''),'decimalLatitude':lat,'decimalLongitude':lon,'year':x.get('year',''),'basisOfRecord':x.get('basisOfRecord','')})
        if i+1<len(spp):time.sleep(a.delay)
    p=Path(a.out);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=FIELDS);w.writeheader();w.writerows(output)
    counts={}
    for r in output:counts[r['canonical_name']]=counts.get(r['canonical_name'],0)+1
    qc={'species_requested':len(spp),'species_with_coordinates':len(counts),'coordinate_rows':len(output),'species_ge20':sum(v>=20 for v in counts.values()),'failed_requests':len(errors),'errors':errors,'method':'GBIF occurrence search, deterministic first-page sample, coordinates deduplicated at 0.001 degrees'}
    Path(a.qc_out).write_text(json.dumps(qc,indent=2)+'\n');print(json.dumps(qc))
if __name__=='__main__':main()
