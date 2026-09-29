#!/usr/bin/env python3
from __future__ import annotations

import argparse, csv, json, math, urllib.parse, urllib.request
from collections import defaultdict
from pathlib import Path


def hav(lat1,lon1,lat2,lon2):
    r=6371.0088
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=math.radians(lat2-lat1); dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1,math.sqrt(a)))


def max_span(coords):
    pts=list(dict.fromkeys(coords))
    if len(pts)<2: return 0.0
    best=0.0
    for i,(a,b) in enumerate(pts):
        for c,d in pts[i+1:]:
            best=max(best,hav(a,b,c,d))
    return best


def ttt(path:Path,out:Path):
    by=defaultdict(lambda:defaultdict(list))
    nrows=0
    with path.open(newline='',encoding='utf-8-sig') as f:
        r=csv.DictReader(f)
        req={'AccSpeciesName','Latitude','Longitude','SiteName','Trait','Value'}
        missing=req-set(r.fieldnames or [])
        if missing: raise RuntimeError(f'TTT missing columns {sorted(missing)}; got {r.fieldnames}')
        for row in r:
            nrows+=1
            sp=(row.get('AccSpeciesName') or '').strip()
            tr=(row.get('Trait') or '').strip()
            try:
                lat=float(row['Latitude']); lon=float(row['Longitude']); float(row['Value'])
            except Exception:
                continue
            if not sp or not tr or not math.isfinite(lat) or not math.isfinite(lon):
                continue
            by[tr][sp].append((lat,lon,(row.get('SiteName') or '').strip()))
    rows=[]
    for tr,species in sorted(by.items()):
        total_species=len(species)
        spatial=0
        for sp,obs in species.items():
            coords=[(a,b) for a,b,_ in obs]
            sites=set((s if s else f'{a:.4f},{b:.4f}') for a,b,s in obs)
            if len(obs)>=10 and len(sites)>=3 and max_span(coords)>=50:
                spatial+=1
        rows.append({'trait':tr,'species_total':total_species,'spatial_eligible_species':spatial,
                     'paired_source_gate': total_species>=20 and spatial>=5})
    passed=[x for x in rows if x['paired_source_gate']]
    result={
      'source':'Tundra Trait Team archived cleaned dataset',
      'status':'TTT_PAIRED_SOURCE_PREFLIGHT_PASS' if passed else 'HOLD_TTT_NO_TRAIT_PASSES_PAIRED_SOURCE_GATE',
      'rows_read':nrows,
      'traits_seen':len(rows),
      'traits_passing_paired_source_gate':len(passed),
      'gate':{'species_total_min':20,'spatial_eligible_species_min':5,'obs_per_species_min':10,'sites_per_species_min':3,'span_km_min':50},
      'trait_coverage':rows,
      'outcomes_computed':False,
      'note':'Coverage only. No phylogenetic or spatial turnover coefficient was calculated.'
    }
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def fetch_json(url):
    req=urllib.request.Request(url,headers={'User-Agent':'fcp-spatiotemporal-preflight/0.1','Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=60) as h:
        return json.loads(h.read().decode('utf-8'))


def mangal(out:Path,max_pages:int=100):
    base='https://mangal.io/api/v2/network'
    nets=[]
    for page in range(1,max_pages+1):
        url=base+'?'+urllib.parse.urlencode({'count':1000,'page':page})
        data=fetch_json(url)
        if isinstance(data,dict):
            if 'data' in data and isinstance(data['data'],list): data=data['data']
            else: data=[data]
        if not data: break
        nets.extend(data)
        if len(data)<1000: break
    if not nets: raise RuntimeError('Mangal network API returned zero rows')
    def present(v): return v is not None and str(v).strip() not in ('','null','None')
    geo=0; dated=0; both=0; dataset_ids=set(); fields=set()
    for n in nets:
        fields.update(n.keys())
        lat=n.get('latitude'); lon=n.get('longitude'); date=n.get('date')
        g=present(lat) and present(lon); d=present(date)
        geo+=int(g); dated+=int(d); both+=int(g and d)
        ds=n.get('dataset_id')
        if ds is None and isinstance(n.get('dataset'),dict): ds=n['dataset'].get('id')
        if present(ds): dataset_ids.add(str(ds))
    result={
      'source':'Mangal API v2 network endpoint',
      'status':'MANGAL_NETWORK_METADATA_PREFLIGHT_PASS' if geo>=100 and both>=50 else 'HOLD_MANGAL_NETWORK_METADATA_COVERAGE',
      'network_rows_retrieved':len(nets),
      'networks_with_coordinates':geo,
      'networks_with_date':dated,
      'networks_with_coordinates_and_date':both,
      'distinct_dataset_ids_observed':len(dataset_ids),
      'fields_observed':sorted(fields),
      'gate':{'networks_with_coordinates_min':100,'networks_with_coordinates_and_date_min':50},
      'interaction_values_opened':False,
      'rewiring_computed':False,
      'note':'Network metadata only. Interactions/nodes were not queried in this preflight.'
    }
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('ttt'); a.add_argument('--csv',type=Path,required=True); a.add_argument('--out',type=Path,required=True)
    b=sub.add_parser('mangal'); b.add_argument('--out',type=Path,required=True); b.add_argument('--max-pages',type=int,default=100)
    x=ap.parse_args()
    if x.cmd=='ttt': ttt(x.csv,x.out)
    else: mangal(x.out,x.max_pages)

if __name__=='__main__': main()
