#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, itertools, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

NREP=200
OBS_SEED=20260909
NULL_SEED=20260911
TOL=1e-12
PROTOCOL_COMMIT='106ffd9986a5a2316a0909df14fb26f0dee28296'
EXPECTED_INPUT_SHA256='0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6'
EXPECTED_LAND_SHA256='3991a92bfb1c1e7d6f89d11750c455370d1633338f8838b4dbe76296bc8212f4'
EXPECTED_REFERENCE_SHA256='d9dfddc4cc18052921c2b2a4262ab09da903277757f9ab251222627282c172ad'
EXPECTED_UPSTREAM_SHA256='131dd3f36d4f6fcdfb53b312c129973cf7949c792813a7072411d31b391fa4e4'


def load_module(path: Path):
    spec=importlib.util.spec_from_file_location('rgfca_recurrence_runner', path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--recurrence-runner', type=Path, default=Path('scripts/analysis/run_rgfca_global_signal_recurrence_20260909.py'))
    p.add_argument('--input', type=Path, default=Path('data/derived/rgfca_reserve_replication_measured_photos_v1.csv'))
    p.add_argument('--land-metadata', type=Path, default=Path('data/derived/rgfca_reserve_island_mainland_metadata.csv'))
    p.add_argument('--reference-loadings', type=Path, default=Path('results/rgfca_global_signal_recurrence_20260909/reference_loadings.csv'))
    p.add_argument('--upstream-realizations', type=Path, default=Path('results/rgfca_global_signal_recurrence_20260909/recurrence_realizations.csv'))
    p.add_argument('--outdir', type=Path, default=Path('results/rgfca_recurrent_signal_ecological_overlay_20260909'))
    p.add_argument('--reps', type=int, default=NREP)
    p.add_argument('--seed', type=int, default=OBS_SEED)
    p.add_argument('--null-seed', type=int, default=NULL_SEED)
    return p.parse_args()


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_inputs(args):
    expected={args.input:EXPECTED_INPUT_SHA256,args.land_metadata:EXPECTED_LAND_SHA256,args.reference_loadings:EXPECTED_REFERENCE_SHA256,args.upstream_realizations:EXPECTED_UPSTREAM_SHA256}
    for path,want in expected.items():
        got=sha256(path)
        if got!=want:
            raise RuntimeError(f'input SHA-256 mismatch for {path}: {got} != {want}')


def slope(x, y, w):
    den=float(np.sum(w*x*x))
    return float(np.sum(w*x*y)/den) if den>0 else np.nan


def main():
    args=parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    verify_inputs(args)
    mod=load_module(args.recurrence_runner)
    df, sp, obs, d_values, x, y = mod.load_measurements(args.input)
    f_values=np.column_stack([df[f'palette_count_{c}'].to_numpy(float)/df['flower_total12'].to_numpy(float) for c in mod.PALETTE])
    b_values=np.column_stack([df[f'background_palette_count_{c}'].to_numpy(float)/df['background_total12'].to_numpy(float) for c in mod.PALETTE])
    if np.max(np.abs(d_values-(f_values-b_values)))>TOL:
        raise RuntimeError('photo-level D=F-B identity failed')

    land=pd.read_csv(args.land_metadata, usecols=['measurement_id','island_primary'])
    if land['measurement_id'].duplicated().any():
        raise RuntimeError('land metadata has duplicate measurement_id')
    lm=df[['measurement_id']].merge(land,on='measurement_id',how='left',validate='one_to_one')
    if lm['island_primary'].isna().any():
        raise RuntimeError('missing island classification for eligible photo')
    island=lm['island_primary'].astype(bool).to_numpy(float)
    abs_lat10=np.abs(df['latitude'].to_numpy(float))/10.0

    groups_for_shift, _, _ = mod.make_engine(sp, obs, d_values, x, y)
    reference=pd.read_csv(args.reference_loadings,index_col=0).to_numpy(float)[:,:5]
    permutations=list(itertools.permutations(range(5)))
    obs_rng=np.random.default_rng(args.seed)
    null_rng=np.random.default_rng(args.null_seed)

    def state_vectors_triplet(inv, eligible):
        order=np.argsort(inv,kind='stable'); sorted_inv=inv[order]
        starts=np.r_[0,np.flatnonzero(np.diff(sorted_inv))+1]; ends=np.r_[starts[1:],len(order)]
        state_ids=[]; ds=[]; fs=[]; bs=[]
        for start,end in zip(starts,ends):
            sid=int(sorted_inv[start])
            if not eligible[sid]: continue
            idx=order[start:end]
            state_obs=obs[idx]; unique_obs=np.unique(state_obs)
            drawn_obs=obs_rng.choice(unique_obs,size=len(unique_obs),replace=True)
            u=obs_rng.random((len(unique_obs),len(idx)))
            u[drawn_obs[:,None]!=state_obs[None,:]]=np.inf
            chosen=idx[np.argmin(u,axis=1)]
            state_ids.append(sid)
            ds.append(d_values[chosen].mean(axis=0))
            fs.append(f_values[chosen].mean(axis=0))
            bs.append(b_values[chosen].mean(axis=0))
        return np.asarray(state_ids,np.int32),np.vstack(ds),np.vstack(fs),np.vstack(bs)

    def state_predictors(inv):
        counts=np.bincount(inv).astype(float)
        lat=np.bincount(inv,weights=abs_lat10,minlength=len(counts))/counts
        isl=np.bincount(inv,weights=island,minlength=len(counts))/counts
        return lat,isl

    def decompose(state_ids,dm,fm,bm,xc,yc,state_species):
        species_id=state_species[state_ids]
        order=np.argsort(species_id,kind='stable')
        species_id=species_id[order]; state_ids=state_ids[order]
        dm,fm,bm=dm[order],fm[order],bm[order]
        starts=np.r_[0,np.flatnonzero(np.diff(species_id))+1]; ends=np.r_[starts[1:],len(species_id)]
        covariance=np.zeros((12,12)); admitted=[]; n_species=0; n_states=0
        for start,end in zip(starts,ends):
            if end-start<2: continue
            ids=state_ids[start:end]
            xx,yy=xc[ids],yc[ids]
            dd=np.sqrt((xx[:,None]-xx[None,:])**2+(yy[:,None]-yy[None,:])**2)
            if dd.max()<mod.SEPARATION_M: continue
            zd=dm[start:end]-dm[start:end].mean(axis=0,keepdims=True)
            zf=fm[start:end]-fm[start:end].mean(axis=0,keepdims=True)
            zb=bm[start:end]-bm[start:end].mean(axis=0,keepdims=True)
            if np.max(np.abs(zd-(zf-zb)))>TOL: raise RuntimeError('state D=F-B failed')
            k=end-start; covariance+=(zd.T@zd)/k
            admitted.append((int(species_id[start]),ids.copy(),zd,zf,zb))
            n_species+=1; n_states+=k
        covariance/=n_species
        evals,evec=np.linalg.eigh(covariance); idx=np.argsort(evals)[::-1]
        evals=evals[idx]; evec=evec[:,idx]; explained=evals/evals.sum()
        return n_species,n_states,evec,explained,admitted

    rows=[]; upstream_checks=[]
    for rep in range(1,args.reps+1):
        sx=float(obs_rng.uniform(0,mod.CELL_M)); sy=float(obs_rng.uniform(0,mod.CELL_M))
        inv,eligible,xc,yc,state_species=groups_for_shift(sx,sy)
        state_lat,state_island=state_predictors(inv)
        state_ids,dm,fm,bm=state_vectors_triplet(inv,eligible)
        n_species,n_states,evec,explained,admitted=decompose(state_ids,dm,fm,bm,xc,yc,state_species)
        first5=evec[:,:5]; sim=np.abs(reference.T@first5)
        perm=max(permutations,key=lambda p:sum(sim[i,p[i]] for i in range(5)))
        matched=np.column_stack([first5[:,perm[i]] for i in range(5)])
        cos=[]
        for i in range(5):
            dot=float(reference[:,i]@matched[:,i])
            if dot<0: matched[:,i]*=-1; dot=-dot
            cos.append(dot)
        upstream_checks.append({'rep':rep,'n_species':n_species,'n_states':n_states,
            **{f'M{i+1}_cosine':cos[i] for i in range(5)},
            **{f'M{i+1}_matched_run_rank':perm[i]+1 for i in range(5)},
            **{f'M{i+1}_explained':explained[perm[i]] for i in range(5)}})

        data={('abs_latitude','x'):[],('abs_latitude','xn'):[],('island_share','x'):[],('island_share','xn'):[],('w','w'):[],('species','id'):[]}
        comp_store={(comp,m):[] for comp in ('F','D','B') for m in range(3)}
        informative={'abs_latitude':0,'island_share':0}
        for spid,ids,zd,zf,zb in admitted:
            k=len(ids)
            lat=state_lat[ids]; isl=state_island[ids]
            latc=lat-lat.mean(); islc=isl-isl.mean()
            if np.sum(latc*latc)>0: informative['abs_latitude']+=1
            if np.sum(islc*islc)>0: informative['island_share']+=1
            latn=lat[null_rng.permutation(k)]; latn=latn-latn.mean()
            isln=isl[null_rng.permutation(k)]; isln=isln-isln.mean()
            data[('abs_latitude','x')].append(latc); data[('abs_latitude','xn')].append(latn)
            data[('island_share','x')].append(islc); data[('island_share','xn')].append(isln)
            data[('w','w')].append(np.full(k,1.0/k)); data[('species','id')].append(np.full(k,spid))
            sd=zd@matched; sf=zf@matched; sb=zb@matched
            for m in range(3):
                comp_store[('F',m)].append(sf[:,m]); comp_store[('D',m)].append(sd[:,m]); comp_store[('B',m)].append(sb[:,m])
        w=np.concatenate(data[('w','w')])
        for pred in ('abs_latitude','island_share'):
            xx=np.concatenate(data[(pred,'x')]); xn=np.concatenate(data[(pred,'xn')])
            for comp in ('F','D','B'):
                for m in range(3):
                    yy=np.concatenate(comp_store[(comp,m)])
                    rows.append({'rep':rep,'mode':f'M{m+1}','component':comp,'predictor':pred,
                        'observed_slope':slope(xx,yy,w),'null_slope':slope(xn,yy,w),
                        'n_informative_species':informative[pred], 'n_species':n_species,'n_states':n_states})

    realization=pd.DataFrame(rows)
    checks=pd.DataFrame(upstream_checks)
    upstream=pd.read_csv(args.upstream_realizations)
    check_cols=['rep','n_species','n_states']+[x for i in range(1,6) for x in (f'M{i}_cosine',f'M{i}_matched_run_rank',f'M{i}_explained')]
    merged=checks[check_cols].merge(upstream[check_cols],on='rep',suffixes=('_new','_up'))
    maxdiff=0.0
    for col in check_cols[1:]:
        maxdiff=max(maxdiff,float(np.max(np.abs(merged[f'{col}_new'].to_numpy(float)-merged[f'{col}_up'].to_numpy(float)))))
    if maxdiff>TOL:
        raise RuntimeError(f'upstream recurrence reproduction mismatch: {maxdiff}')

    summary=[]
    for keys,g in realization.groupby(['mode','component','predictor'],sort=True):
        mode,comp,pred=keys; obs_s=g['observed_slope']; nul=g['null_slope']; obsmed=float(obs_s.median())
        p=(1+int((np.abs(nul.to_numpy())>=abs(obsmed)).sum()))/(len(nul)+1)
        summary.append({'mode':mode,'component':comp,'predictor':pred,
            'observed_median':obsmed,'observed_q025':float(obs_s.quantile(.025)),'observed_q975':float(obs_s.quantile(.975)),
            'fraction_observed_gt0':float((obs_s>0).mean()),
            'null_median':float(nul.median()),'null_q025':float(nul.quantile(.025)),'null_q975':float(nul.quantile(.975)),
            'observed_minus_null_median':float((obs_s-nul).median()),
            'p_like_two_sided':p,
            'informative_species_median':float(g['n_informative_species'].median()),
            'informative_species_min':int(g['n_informative_species'].min()),
            'informative_species_max':int(g['n_informative_species'].max())})
    summary=pd.DataFrame(summary)
    global_summary=pd.DataFrame([{'n_reps':args.reps,'observed_seed':args.seed,'null_seed':args.null_seed,
        'n_eligible_photos':len(df),'n_eligible_species':df['species'].nunique(),'upstream_max_abs_diff':maxdiff,
        'species_median':checks['n_species'].median(),'states_median':checks['n_states'].median()}])
    realization.to_csv(args.outdir/'ecological_overlay_realizations.csv',index=False)
    summary.to_csv(args.outdir/'ecological_overlay_summary.csv',index=False)
    checks.to_csv(args.outdir/'upstream_mode_reproduction_check.csv',index=False)
    global_summary.to_csv(args.outdir/'ecological_overlay_global_summary.csv',index=False)
    print(global_summary.to_string(index=False))
    print('\nPRIMARY F COMPONENT')
    print(summary[summary.component.eq('F')].to_string(index=False))
    print('\nALL strongest by p-like')
    print(summary.sort_values(['p_like_two_sided','mode','component','predictor']).head(18).to_string(index=False))

if __name__=='__main__': main()
