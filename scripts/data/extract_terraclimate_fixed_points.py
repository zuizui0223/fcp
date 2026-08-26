#!/usr/bin/env python3
from __future__ import annotations
import argparse, io, json, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
import requests

DEFAULT_VARS=['tmin','tmax','ppt','def','vpd']
DEFAULT_BASE='http://thredds.northwestknowledge.net:8080/thredds/ncss/agg_terraclimate_{var}_1950_CurrentYear_GLOBE.nc'

def request_one(var,lat,lon,start_year,end_year,base,timeout,retries):
    params={'var':var,'latitude':f'{lat:.8f}','longitude':f'{lon:.8f}',
            'time_start':f'{start_year}-01-01T00:00:00Z','time_end':f'{end_year}-12-31T23:59:59Z','accept':'csv'}
    last=None
    for attempt in range(retries):
        try:
            r=requests.get(base.format(var=var),params=params,timeout=timeout)
            r.raise_for_status()
            d=pd.read_csv(io.StringIO(r.text))
            if d.empty: raise RuntimeError('empty NCSS response')
            time_col=next((c for c in d.columns if c.lower() in {'date','time'} or c.lower().startswith('date')),None)
            value_col=next((c for c in d.columns if c==var or c.lower().startswith(var.lower()+'[')),None)
            if time_col is None or value_col is None: raise RuntimeError(f'unexpected columns: {list(d.columns)}')
            z=pd.DataFrame({'time':pd.to_datetime(d[time_col],utc=True,errors='coerce'),var:pd.to_numeric(d[value_col],errors='coerce')}).dropna(subset=['time'])
            z=z.loc[z.time.dt.year.between(start_year,end_year)].drop_duplicates('time').sort_values('time')
            return z
        except Exception as e:
            last=e
            if attempt+1<retries: time.sleep(min(2**attempt,16))
    raise RuntimeError(f'{var}@({lat},{lon}) failed after {retries} attempts: {last}')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--points',required=True); p.add_argument('--out',required=True); p.add_argument('--qc-out',required=True); p.add_argument('--cache-dir',required=True)
    p.add_argument('--start-year',type=int,default=1958); p.add_argument('--end-year',type=int,default=2025); p.add_argument('--variables',nargs='+',default=DEFAULT_VARS)
    p.add_argument('--base-template',default=DEFAULT_BASE); p.add_argument('--workers',type=int,default=12); p.add_argument('--timeout',type=int,default=120); p.add_argument('--retries',type=int,default=4)
    a=p.parse_args(); pts=pd.read_csv(a.points)
    required={'canonical_name','family','selection_rank','terraclimate_lat','terraclimate_lon','tc_lat_idx','tc_lon_idx'}
    if not required<=set(pts): raise SystemExit(f'missing point columns: {sorted(required-set(pts))}')
    cache=Path(a.cache_dir); cache.mkdir(parents=True,exist_ok=True)
    unique=pts[['tc_lat_idx','tc_lon_idx','terraclimate_lat','terraclimate_lon']].drop_duplicates().sort_values(['tc_lat_idx','tc_lon_idx'])
    jobs=[]
    for r in unique.itertuples(index=False):
        for var in a.variables:
            cp=cache/f'{var}_{int(r.tc_lat_idx)}_{int(r.tc_lon_idx)}_{a.start_year}_{a.end_year}.csv'
            jobs.append((var,float(r.terraclimate_lat),float(r.terraclimate_lon),int(r.tc_lat_idx),int(r.tc_lon_idx),cp))
    failures=[]
    def run(job):
        var,lat,lon,ii,jj,cp=job
        if cp.exists() and cp.stat().st_size>50:
            return job,'cached'
        z=request_one(var,lat,lon,a.start_year,a.end_year,a.base_template,a.timeout,a.retries)
        z.to_csv(cp,index=False); return job,'downloaded'
    counts={'cached':0,'downloaded':0}
    with ThreadPoolExecutor(max_workers=max(1,a.workers)) as ex:
        futs={ex.submit(run,j):j for j in jobs}
        for i,f in enumerate(as_completed(futs),1):
            try:
                _,status=f.result(); counts[status]+=1
            except Exception as e:
                j=futs[f]; failures.append({'var':j[0],'tc_lat_idx':j[3],'tc_lon_idx':j[4],'error':str(e)[:500]})
            if i%250==0 or i==len(futs): print({'jobs_done':i,'jobs_total':len(futs),'failures':len(failures)},flush=True)
    if failures:
        Path(a.qc_out).write_text(json.dumps({'status':'failed','failures':failures,'counts':counts},indent=2)+'\n')
        raise SystemExit(f'{len(failures)} TerraClimate point-variable requests failed')
    pieces=[]
    for r in unique.itertuples(index=False):
        merged=None
        for var in a.variables:
            cp=cache/f'{var}_{int(r.tc_lat_idx)}_{int(r.tc_lon_idx)}_{a.start_year}_{a.end_year}.csv'
            z=pd.read_csv(cp,parse_dates=['time'])
            merged=z if merged is None else merged.merge(z,on='time',how='outer',validate='one_to_one')
        merged['tc_lat_idx']=int(r.tc_lat_idx); merged['tc_lon_idx']=int(r.tc_lon_idx)
        merged['terraclimate_lat']=float(r.terraclimate_lat); merged['terraclimate_lon']=float(r.terraclimate_lon)
        pieces.append(merged)
    climate=pd.concat(pieces,ignore_index=True)
    out=pts.merge(climate,on=['tc_lat_idx','tc_lon_idx','terraclimate_lat','terraclimate_lon'],how='left',validate='many_to_many')
    out['time']=pd.to_datetime(out.time,utc=True); out['year']=out.time.dt.year; out['month']=out.time.dt.month
    expected_months=(a.end_year-a.start_year+1)*12
    per=out.groupby(['canonical_name','selection_rank']).time.nunique()
    qc={'status':'complete','protocol':'terraclimate-fixed-points-v1','terraclimate_version_intent':'v1.1/current aggregated service; do not mix with legacy v1.0','start_year':a.start_year,'end_year':a.end_year,'variables':a.variables,'fixed_points':len(pts),'unique_grid_cells':len(unique),'expected_months_per_point':expected_months,'min_months_per_point':int(per.min()),'max_months_per_point':int(per.max()),'rows':len(out),'requests':len(jobs),'cache_counts':counts,'source':a.base_template,'method_guard':'fixed spatial locations are followed through climate time; GBIF occurrence year is not used as the temporal sampling process'}
    if per.min()<expected_months: qc['warning']='some point-month histories are incomplete'
    op=Path(a.out); op.parent.mkdir(parents=True,exist_ok=True); out.to_csv(op,index=False); Path(a.qc_out).write_text(json.dumps(qc,indent=2)+'\n'); print(json.dumps(qc,indent=2))
if __name__=='__main__': main()
