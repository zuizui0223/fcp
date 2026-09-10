#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
DISC=ROOT/'data/derived/global_monte_carlo_measured_photos_v1.csv'
RES=ROOT/'data/derived/rgfca_reserve_replication_measured_photos_v1.csv'
MET=ROOT/'results/rgfca_worldclim_old34_transfer_step7d_20260910/species_environment_metrics.csv'
OUT=ROOT/'results/rgfca_state_construction_null_step7f_20260910'; OUT.mkdir(parents=True,exist_ok=True)
MORPHS=['white','yellow_orange','red_pink','blue_purple']; N=999; SEED=20260910
OBS={'discovery':9.37182419530179,'reserve':9.661780091846277}


def bseries(s):
    if pd.api.types.is_bool_dtype(s): return s.fillna(False).astype(bool)
    return s.fillna('').astype(str).str.lower().isin(['true','1','yes'])


def hav_pairs(lat,lon):
    i,j=np.triu_indices(len(lat),1); la1=np.radians(lat[i]); la2=np.radians(lat[j]); dl=np.radians(lon[j]-lon[i]); dp=la2-la1
    a=np.sin(dp/2)**2+np.cos(la1)*np.cos(la2)*np.sin(dl/2)**2
    d=6371.0088*2*np.arcsin(np.minimum(1,np.sqrt(a))); return i,j,d


def load(path):
    use=['species','morph','global_classifiable','latitude','longitude','observer_id']
    d=pd.read_csv(path,usecols=use); keep=bseries(d.global_classifiable)&d.morph.isin(MORPHS)
    return d.loc[keep].dropna(subset=['latitude','longitude']).copy()


def classify_bank(g,rng):
    n=len(g); counts=g.morph.value_counts(); p=np.array([counts.get(m,0)/n for m in MORPHS],float)
    if n<40 or np.sort(p)[-2]<.10: return None
    lat=g.latitude.to_numpy(float); lon=g.longitude.to_numpy(float); morph=g.morph.astype(str).to_numpy(); obs=g.observer_id.astype(str).to_numpy(); ii,jj,dist=hav_pairs(lat,lon)
    cod=(dist<=100.0)&(obs[ii]!=obs[jj])
    cal=np.empty(N); ev_delta=np.empty(N); ev_C=np.zeros(N,bool)
    for b in range(N):
        mp=morph[rng.permutation(n)]; diff=mp[ii]!=mp[jj]
        cal[b]=np.median(dist[diff])-np.median(dist[~diff])
    for b in range(N):
        mp=morph[rng.permutation(n)]; diff=mp[ii]!=mp[jj]; ev_delta[b]=np.median(dist[diff])-np.median(dist[~diff])
        cross=cod&diff; endpoints=np.unique(np.concatenate([ii[cross],jj[cross]])) if cross.any() else np.array([],int); observers=np.unique(np.concatenate([obs[ii[cross]],obs[jj[cross]]])) if cross.any() else np.array([],str)
        ev_C[b]=bool(cod.sum()>=30 and cross.sum()>=5 and len(endpoints)>=4 and len(observers)>=3)
    pvals=np.array([(1+np.sum(cal>=x-1e-12))/(N+1) for x in ev_delta]); ev_S=(ev_delta>0)&(pvals<=.05)
    return ev_C,ev_S


def run_tranche(name,path,seed):
    d=load(path); met=pd.read_csv(MET); met=met[(met.tranche==name)&(met.metric_status=='complete')][['species','hv3d_rarefied20']].dropna(); hv=dict(zip(met.species,met.hv3d_rarefied20))
    rng=np.random.default_rng(seed); banks={}
    for sp,g in d.groupby('species',sort=True):
        if sp not in hv: continue
        r=classify_bank(g.reset_index(drop=True),rng)
        if r is not None: banks[sp]=r
    deltas=[]; nC=[]; nS=[]
    species=sorted(banks)
    for b in range(N):
        cvals=[]; svals=[]
        for sp in species:
            C,S=banks[sp][0][b],banks[sp][1][b]
            if C and not S: cvals.append(hv[sp])
            elif S and not C: svals.append(hv[sp])
        if cvals and svals:
            deltas.append(float(np.median(svals)-np.median(cvals))); nC.append(len(cvals)); nS.append(len(svals))
    arr=np.asarray(deltas,float); obs=OBS[name]
    p=float((1+np.sum(arr>=obs-1e-12))/(len(arr)+1)) if len(arr) else 1.0
    out={'tranche':name,'eligible_species':len(species),'valid_null_worlds':int(len(arr)),'observed_H0_delta':obs,'null_median':float(np.median(arr)) if len(arr) else None,'null_q005':float(np.quantile(arr,.005)) if len(arr) else None,'null_q025':float(np.quantile(arr,.025)) if len(arr) else None,'null_q975':float(np.quantile(arr,.975)) if len(arr) else None,'null_q995':float(np.quantile(arr,.995)) if len(arr) else None,'one_sided_p':p,'null_pure_C_median':float(np.median(nC)) if nC else None,'null_pure_S_median':float(np.median(nS)) if nS else None}
    pd.DataFrame({'null_H0_delta':arr,'pure_C_count':nC,'pure_S_count':nS}).to_csv(OUT/f'{name}_null_worlds.csv',index=False)
    return out


def main():
    ds=run_tranche('discovery',DISC,SEED+1); rs=run_tranche('reserve',RES,SEED+2)
    above=lambda z: z['valid_null_worlds']>0 and z['observed_H0_delta']>z['null_q975'] and z['one_sided_p']<=.05
    within=lambda z: z['valid_null_worlds']>0 and z['observed_H0_delta']<=z['null_q975']
    cls='construction_insufficient' if above(ds) and above(rs) else ('construction_sufficient' if within(ds) and within(rs) else 'mixed_or_unresolved')
    result={'analysis':'rgfca_state_construction_null_step7f','discovery':ds,'reserve':rs,'diagnostic_class':cls,'boundary':'Post-result classifier-null diagnostic; no causal claim.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# RGFCA Step 7F — state-construction null diagnostic','']
    for label,z in [('Discovery',ds),('Reserve',rs)]:
        lines += [f'## {label}',f"- eligible species: **{z['eligible_species']}**",f"- valid null worlds: **{z['valid_null_worlds']} / {N}**",f"- observed H0 delta: **{z['observed_H0_delta']:.4f}**",f"- classifier-null median: **{z['null_median']:.4f}**",f"- classifier-null 95% interval: **[{z['null_q025']:.4f}, {z['null_q975']:.4f}]**",f"- classifier-null 99% interval: **[{z['null_q005']:.4f}, {z['null_q995']:.4f}]**",f"- one-sided p(null >= observed): **{z['one_sided_p']:.6g}**",f"- median null pure-state counts: **{z['null_pure_C_median']:.1f} C* / {z['null_pure_S_median']:.1f} S***",'']
    lines += [f'## Diagnostic class\n\n**{cls}**','', 'Morph labels are shuffled within species while coordinates, observers, morph counts, and real hypervolume values are fixed.']
    (OUT/'RESULT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__': main()
