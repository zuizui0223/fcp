#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, itertools
from pathlib import Path
import numpy as np, pandas as pd

RUNNER=Path('scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py')
INPUT=Path('data/derived/rgfca_reserve_replication_measured_photos_v1.csv')
REF=Path('results/rgfca_global_signal_recurrence_20260909/reference_loadings.csv')
UP=Path('results/rgfca_global_signal_recurrence_20260909/recurrence_realizations.csv')
OUT=Path('results/rgfca_flower_side_survivor_season_match_20260909')
RX,RY,MODE=-23,11,1
OBS_SEED=20260909; NREP=200; REPORT=500_000.0; TOL=1e-12; MIN_MATCHED_SPECIES=10
SCENARIOS=[('primary_season30_year5',30.0,5.0),('sensitivity_season30_no_year',30.0,None),('sensitivity_season30_year2',30.0,2.0)]

def load(path):
    spec=importlib.util.spec_from_file_location('m',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def circular_mean_day(doys):
    a=2*np.pi*(np.asarray(doys,float)-1.0)/365.25;s=np.mean(np.sin(a));c=np.mean(np.cos(a));ang=np.arctan2(s,c)
    if ang<0:ang+=2*np.pi
    return float(1.0+ang*365.25/(2*np.pi))

def circular_day_distance(a,b):
    d=abs(float(a)-float(b));return min(d,365.25-d)

def main():
    OUT.mkdir(parents=True,exist_ok=True);m=load(RUNNER);df,sp,obs,dv,x,y=m.load_measurements(INPUT)
    raw=pd.read_csv(INPUT,usecols=['measurement_id','observed_on']);dmap=raw.drop_duplicates('measurement_id').set_index('measurement_id')['observed_on'];dates=pd.to_datetime(df['measurement_id'].map(dmap),errors='coerce')
    if dates.isna().any():raise RuntimeError(f'missing dates: {int(dates.isna().sum())}')
    years=dates.dt.year.to_numpy(float);doys=dates.dt.dayofyear.to_numpy(float)
    fv=np.column_stack([df[f'palette_count_{c}'].to_numpy(float)/df['flower_total12'].to_numpy(float) for c in m.PALETTE]);bv=np.column_stack([df[f'background_palette_count_{c}'].to_numpy(float)/df['background_total12'].to_numpy(float) for c in m.PALETTE]);assert np.max(np.abs(dv-(fv-bv)))<TOL
    groups,_,_=m.make_engine(sp,obs,dv,x,y);ref=pd.read_csv(REF,index_col=0).to_numpy(float)[:,:5];perms=list(itertools.permutations(range(5)));rng=np.random.default_rng(OBS_SEED)
    def vectors_and_dates(inv,elig):
        order=np.argsort(inv,kind='stable');s=inv[order];starts=np.r_[0,np.flatnonzero(np.diff(s))+1];ends=np.r_[starts[1:],len(order)];ids=[];ds=[];fs=[];bs=[];state_doy=[];state_year=[]
        for a,b in zip(starts,ends):
            sid=int(s[a]);
            if not elig[sid]:continue
            idx=order[a:b];state_doy.append(circular_mean_day(doys[idx]));state_year.append(float(np.mean(years[idx])));oo=obs[idx];uo=np.unique(oo);draw=rng.choice(uo,size=len(uo),replace=True);u=rng.random((len(uo),len(idx)));u[draw[:,None]!=oo[None,:]]=np.inf;ch=idx[np.argmin(u,axis=1)];ids.append(sid);ds.append(dv[ch].mean(0));fs.append(fv[ch].mean(0));bs.append(bv[ch].mean(0))
        return np.asarray(ids,np.int32),np.vstack(ds),np.vstack(fs),np.vstack(bs),np.asarray(state_doy),np.asarray(state_year)
    def decomp(ids,dm,fm,bm,sdoy,syear,xc,yc,ssp):
        sid=ssp[ids];order=np.argsort(sid,kind='stable');sid=sid[order];ids=ids[order];dm,fm,bm=dm[order],fm[order],bm[order];sdoy=sdoy[order];syear=syear[order];starts=np.r_[0,np.flatnonzero(np.diff(sid))+1];ends=np.r_[starts[1:],len(sid)];cov=np.zeros((12,12));adm=[];ns=nst=0
        for a,b in zip(starts,ends):
            if b-a<2:continue
            ii=ids[a:b];xx,yy=xc[ii],yc[ii]
            if np.sqrt((xx[:,None]-xx[None,:])**2+(yy[:,None]-yy[None,:])**2).max()<m.SEPARATION_M:continue
            zd=dm[a:b]-dm[a:b].mean(0);zf=fm[a:b]-fm[a:b].mean(0);zb=bm[a:b]-bm[a:b].mean(0);assert np.max(np.abs(zd-(zf-zb)))<TOL;k=b-a;cov+=(zd.T@zd)/k;adm.append((int(sid[a]),ii.copy(),zf,sdoy[a:b].copy(),syear[a:b].copy()));ns+=1;nst+=k
        cov/=ns;ev,ec=np.linalg.eigh(cov);ix=np.argsort(ev)[::-1];ev,ec=ev[ix],ec[:,ix];return ns,nst,ec,ev/ev.sum(),adm
    rows=[];checks=[];species_rows=[]
    for rep in range(1,NREP+1):
        sx=float(rng.uniform(0,m.CELL_M));sy=float(rng.uniform(0,m.CELL_M));inv,elig,xc,yc,ssp=groups(sx,sy);ids,dm,fm,bm,sdoy,syear=vectors_and_dates(inv,elig);ns,nst,ec,expl,adm=decomp(ids,dm,fm,bm,sdoy,syear,xc,yc,ssp);first=ec[:,:5];sim=np.abs(ref.T@first);perm=max(perms,key=lambda p:sum(sim[i,p[i]] for i in range(5)));matched=np.column_stack([first[:,perm[i]] for i in range(5)]);cos=[]
        for i in range(5):
            dot=float(ref[:,i]@matched[:,i]);
            if dot<0:matched[:,i]*=-1;dot=-dot
            cos.append(dot)
        checks.append({'rep':rep,'n_species':ns,'n_states':nst,**{f'M{i+1}_cosine':cos[i] for i in range(5)},**{f'M{i+1}_matched_run_rank':perm[i]+1 for i in range(5)},**{f'M{i+1}_explained':expl[perm[i]] for i in range(5)}});bysp={}
        for spid,ii,zf,dd,yy in adm:
            sf=zf@matched;bysp[int(spid)]=[{'stateid':int(stateid),'x':float(xc[stateid]),'y':float(yc[stateid]),'rx':int(np.floor(xc[stateid]/REPORT)),'ry':int(np.floor(yc[stateid]/REPORT)),'score':float(sf[j,MODE]),'doy':float(dd[j]),'year':float(yy[j])} for j,stateid in enumerate(ii)]
        for scenario,day_max,year_max in SCENARIOS:
            sp_contrasts=[]
            for spid,recs in bysp.items():
                cands=[r for r in recs if r['rx']==RX and r['ry']==RY]
                if not cands:continue
                cand_contrasts=[]
                for c in cands:
                    comps=[]
                    for q in recs:
                        if q['rx']==RX and q['ry']==RY:continue
                        if np.hypot(q['x']-c['x'],q['y']-c['y'])<m.SEPARATION_M:continue
                        if circular_day_distance(q['doy'],c['doy'])>day_max:continue
                        if year_max is not None and abs(q['year']-c['year'])>year_max:continue
                        comps.append(q['score'])
                    if comps:cand_contrasts.append(c['score']-float(np.mean(comps)))
                if cand_contrasts:
                    val=float(np.mean(cand_contrasts));sp_contrasts.append(val);species_rows.append({'rep':rep,'scenario':scenario,'spid':spid,'n_candidate_states':len(cands),'n_matched_candidate_states':len(cand_contrasts),'species_contrast':val})
            nsp=len(sp_contrasts);supported=nsp>=MIN_MATCHED_SPECIES;rows.append({'rep':rep,'scenario':scenario,'n_matched_species':nsp,'supported':int(supported),'matched_contrast':float(np.mean(sp_contrasts)) if supported else np.nan})
    checks=pd.DataFrame(checks);up=pd.read_csv(UP);cols=['rep','n_species','n_states']+[z for i in range(1,6) for z in (f'M{i}_cosine',f'M{i}_matched_run_rank',f'M{i}_explained')];mm=checks[cols].merge(up[cols],on='rep',suffixes=('_n','_u'));md=0.0
    for c in cols[1:]:md=max(md,float(np.max(np.abs(mm[f'{c}_n'].to_numpy(float)-mm[f'{c}_u'].to_numpy(float)))))
    if md>TOL:raise RuntimeError(f'upstream mismatch {md}')
    res=pd.DataFrame(rows);spr=pd.DataFrame(species_rows);summaries=[]
    for scenario,_,_ in SCENARIOS:
        q=res[(res.scenario==scenario)&(res.supported==1)].copy();vals=q.matched_contrast.dropna();summaries.append({'scenario':scenario,'supported_reps':len(q),'matched_species_median':q.n_matched_species.median() if len(q) else np.nan,'matched_species_min':q.n_matched_species.min() if len(q) else np.nan,'matched_species_max':q.n_matched_species.max() if len(q) else np.nan,'contrast_median':vals.median() if len(vals) else np.nan,'contrast_q025':vals.quantile(.025) if len(vals) else np.nan,'contrast_q975':vals.quantile(.975) if len(vals) else np.nan,'positive_fraction':(vals>0).mean() if len(vals) else np.nan,'robust_gate':bool(scenario=='primary_season30_year5' and len(q)>=100 and (vals>0).mean()>=.70),'upstream_maxdiff':md})
    pd.DataFrame(summaries).to_csv(OUT/'season_match_summary.csv',index=False);res.to_csv(OUT/'season_match_realizations.csv',index=False);spr.to_csv(OUT/'season_match_species_contrasts.csv',index=False);checks.to_csv(OUT/'upstream_mode_reproduction_check.csv',index=False)
if __name__=='__main__':main()
