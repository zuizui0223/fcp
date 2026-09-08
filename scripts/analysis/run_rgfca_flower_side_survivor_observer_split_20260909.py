#!/usr/bin/env python3
from __future__ import annotations
import hashlib, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

RUNNER=Path('scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py')
INPUT=Path('data/derived/rgfca_reserve_replication_measured_photos_v1.csv')
UP=Path('results/rgfca_global_signal_recurrence_20260909/recurrence_realizations.csv')
OUT=Path('results/rgfca_flower_side_survivor_observer_split_20260909')
RX,RY=-23,11; REPORT=500_000.; MINSP=5; TOL=1e-12
SEEDS={'A':20260912,'B':20260913}

def load(path):
    spec=importlib.util.spec_from_file_location('m',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    m=load(RUNNER);df,sp,obs,dv,x,y=m.load_measurements(INPUT)
    fv=np.column_stack([df[f'palette_count_{c}'].to_numpy(float)/df['flower_total12'].to_numpy(float) for c in m.PALETTE])
    bv=np.column_stack([df[f'background_palette_count_{c}'].to_numpy(float)/df['background_total12'].to_numpy(float) for c in m.PALETTE]);assert np.max(np.abs(dv-(fv-bv)))<TOL
    oids=df['observer_id'].astype(str).to_numpy();halves=np.array(['A' if int(hashlib.sha256(s.encode()).hexdigest(),16)%2==0 else 'B' for s in oids])
    up=pd.read_csv(UP);m2cols=[f'M2_loading_{c}' for c in m.PALETTE];rows=[]
    for half in ('A','B'):
        idx=np.flatnonzero(halves==half);spp=sp[idx];xx=x[idx];yy=y[idx];ff=fv[idx];oo_raw=oids[idx]
        oo=pd.factorize(pd.Series(oo_raw),sort=True)[0].astype(np.int32);rng=np.random.default_rng(SEEDS[half]);base=128;cx_offset=64;cy_offset=64;observer_base=int(oo.max())+2
        def build(shift_x,shift_y):
            cx=np.floor((xx+shift_x)/m.CELL_M).astype(np.int16);cy=np.floor((yy+shift_y)/m.CELL_M).astype(np.int16);key=((spp.astype(np.int64)*base+(cx+cx_offset).astype(np.int64))*base+(cy+cy_offset).astype(np.int64));unique_key,first,inv,count=np.unique(key,return_index=True,return_inverse=True,return_counts=True);so=key*observer_base+oo.astype(np.int64);uso=np.unique(so);sk=uso//observer_base;si=np.searchsorted(unique_key,sk);nobs=np.bincount(si,minlength=len(unique_key)).astype(np.int16);elig=(count>=m.MIN_PHOTOS)&(nobs>=m.MIN_OBSERVERS);xc=np.bincount(inv,weights=xx,minlength=len(unique_key))/count;yc=np.bincount(inv,weights=yy,minlength=len(unique_key))/count;ssp=spp[first];order=np.argsort(inv,kind='stable');sinv=inv[order];starts=np.r_[0,np.flatnonzero(np.diff(sinv))+1];ends=np.r_[starts[1:],len(order)];ids=[];fm=[]
            for a,b in zip(starts,ends):
                sid=int(sinv[a]);
                if not elig[sid]:continue
                ii=order[a:b];state_obs=oo[ii];uo=np.unique(state_obs);draw=rng.choice(uo,size=len(uo),replace=True);u=rng.random((len(uo),len(ii)));u[draw[:,None]!=state_obs[None,:]]=np.inf;ch=ii[np.argmin(u,axis=1)];ids.append(sid);fm.append(ff[ch].mean(0))
            return np.asarray(ids,np.int32),np.vstack(fm),xc,yc,ssp
        for _,ur in up.iterrows():
            rep=int(ur.rep);ids,fm,xc,yc,ssp=build(float(ur.shift_x_m),float(ur.shift_y_m));sid=ssp[ids];order=np.argsort(sid,kind='stable');sid=sid[order];ids=ids[order];fm=fm[order];starts=np.r_[0,np.flatnonzero(np.diff(sid))+1];ends=np.r_[starts[1:],len(sid)];scores=[];m2=ur[m2cols].to_numpy(float)
            for a,b in zip(starts,ends):
                if b-a<2:continue
                ii=ids[a:b];gx=xc[ii];gy=yc[ii];d=np.sqrt((gx[:,None]-gx[None,:])**2+(gy[:,None]-gy[None,:])**2)
                if d.max()<m.SEPARATION_M:continue
                zf=fm[a:b]-fm[a:b].mean(0);sc=zf@m2
                for j,stateid in enumerate(ii):
                    rx=int(np.floor(xc[stateid]/REPORT));ry=int(np.floor(yc[stateid]/REPORT))
                    if rx==RX and ry==RY:scores.append({'spid':int(sid[a]),'score':float(sc[j])})
            if scores:
                z=pd.DataFrame(scores).groupby('spid',as_index=False).score.mean();nsp=z.spid.nunique()
                if nsp>=MINSP:
                    val=float(z.score.mean());rows.append({'rep':rep,'half':half,'n_species':nsp,'flower_score':val,'positive':int(val>0)})
    r=pd.DataFrame(rows);summ=[]
    for half in ('A','B'):
        g=r[r.half==half].copy();summ.append({'half':half,'support_reps':len(g),'n_species_median':g.n_species.median(),'n_species_min':g.n_species.min(),'n_species_max':g.n_species.max(),'median_flower_score':g.flower_score.median(),'q025':g.flower_score.quantile(.025),'q975':g.flower_score.quantile(.975),'positive_fraction':g.positive.mean(),'passes_half_gate':bool(len(g)>=100 and g.positive.mean()>=.70)})
    s=pd.DataFrame(summ);s['observer_disjoint_robust']=bool(s.passes_half_gate.all());r.to_csv(OUT/'observer_split_realizations.csv',index=False);s.to_csv(OUT/'observer_split_summary.csv',index=False)
if __name__=='__main__':main()
