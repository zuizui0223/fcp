#!/usr/bin/env python3
"""Source-locked Sisyrinchium M-species line-conditioned clustering and congener exposure.

A retrospective fixed-inventory observational benchmark; NOT species independence,
spatial exchangeability evidence, or a causal selection test.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
from collections import OrderedDict
import numpy as np
import pandas as pd

SOURCE_MD5='008ec196aa5acb769dbfd36e835e32ba'
FREQUENCY_SHA256='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
R=9999
SEEDS={'exposure':2026101309,'distance':2026101310,'clumping':2026101311,'latest_exposure':2026101320,'latest_clumping':2026101321,'latest_distance':2026101322}


def load_source(path:Path)->pd.DataFrame:
    bs=Path(path).read_bytes()
    if hashlib.sha256(bs).hexdigest()!=FREQUENCY_SHA256:
        raise ValueError('ORIGINAL_SISYRINCHIUM_FREQUENCY_SHA256_MISMATCH')
    d=pd.read_csv(path,sep='\t')
    cols={'year','population','line','quad','distance','humidity','v','w','l','s'}
    if set(d.columns)!=cols or len(d)!=512:raise ValueError('SOURCE_SCHEMA_OR_ROWS_CHANGED')
    for col in ['v','w','l','s']:
        if d[col].isna().any() or (d[col]<0).any() or (d[col]%1!=0).any():
            raise ValueError('INVALID_SPECIES_COUNT_'+col)
        d[col]=d[col].astype(int)
    if d[['v','w','l','s']].sum().to_dict()!={'v':5081,'w':8723,'l':4463,'s':29}:
        raise ValueError('SOURCE_SPECIES_INVENTORIES_CHANGED')
    if d.duplicated(['year','population','line','quad']).any():
        raise ValueError('DUPLICATE_SOURCE_QUADRAT')
    d['n']=d.v+d.w
    denom=d.n+d.l+d.s
    d['a']=np.divide(d.l+d.s,denom,out=np.full(len(d),np.nan),where=denom>0)
    return d


def informative_exposure_line(n,v,a,min_line_n=10)->bool:
    n=np.asarray(n,dtype=int);v=np.asarray(v,dtype=int);a=np.asarray(a,dtype=float)
    occ=n>0
    return bool(np.count_nonzero(occ)>=2 and sum(n)>=min_line_n and sum(v[occ])>=1
                and sum(n[occ]-v[occ])>=1 and np.all(np.isfinite(a[occ]))
                and np.ptp(a[occ])>1e-12)


def line_score(n,v,a):
    n=np.asarray(n,dtype=float);v=np.asarray(v,dtype=float);a=np.asarray(a,dtype=float)
    occ=n>0;n=n[occ];v=v[occ];a=a[occ]
    N=n.sum();V=v.sum();p=V/N
    if V<=0 or V>=N:raise ValueError('NO_WITHIN_LINE_COLOUR_VARIATION')
    return float(((v-n*p)*a).sum()/(N*p*(1-p)))


def distance_residual(n,d,a):
    n=np.asarray(n,dtype=float);d=np.asarray(d,dtype=float);a=np.asarray(a,dtype=float)
    if (n<0).any() or not np.all(np.isfinite(a)) or not np.all(np.isfinite(d)):
        raise ValueError('INVALID_DISTANCE_OR_EXPOSURE')
    if sum(n)<=0:raise ValueError('EMPTY_LINE')
    avgd=float(np.average(d,weights=n));avga=float(np.average(a,weights=n))
    cov=float(np.dot(n*(d-avgd),a-avga));var=float(np.dot(n*(d-avgd),d-avgd))
    beta=cov/var if var>0 else 0.0
    return a-avga-beta*(d-avgd)


def B_and_line_null(n,v,line):
    n=np.asarray(n,dtype=int);v=np.asarray(v,dtype=int);line=np.asarray(line).astype(str)
    if n.shape!=v.shape or n.shape!=line.shape or (n<0).any() or (v<0).any() or (v>n).any():
        raise ValueError('BAD_SOURCE_COUNTS')
    occ=n>0;n=n[occ];v=v[occ];line=line[occ]
    N=int(n.sum());V=int(v.sum());W=N-V
    if N<=1 or V==0 or W==0:raise ValueError('MONOMORPHIC_OR_INSUFFICIENT')
    p=V/N
    den=N*p*(1-p)
    B=float(np.sum((v-n*p)**2/n)/den)
    between=0.0;within_expected=0.0
    for name in sorted(set(line)):
        is_line=line==name
        nn=n[is_line];vv=v[is_line]
        Nl=int(nn.sum());Vl=int(vv.sum());pl=Vl/Nl;ml=len(nn)
        between+=Nl*(pl-p)**2/den
        if Nl>1 and ml>1:
            within_expected+=(Nl*pl*(1-pl)*(ml-1)/(Nl-1))/den
    return B,float(between),float(between+within_expected)


def _stratify(n,v,a,line,distance=None,min_line_n=10):
    n=np.asarray(n,dtype=int);v=np.asarray(v,dtype=int);a=np.asarray(a,dtype=float)
    line=np.asarray(line).astype(str)
    if distance is not None:distance=np.asarray(distance,dtype=float)
    groups=[]
    for name in sorted(set(line)):
        ix=np.flatnonzero((line==name)&(n>0))
        nn=n[ix];vv=v[ix];aa=a[ix]
        if len(ix)==0:continue
        info=informative_exposure_line(nn,vv,aa,min_line_n=min_line_n)
        dist_ok=False;ar=None
        if distance is not None and info and len(ix)>=3 and np.ptp(distance[ix])>0:
            ar=distance_residual(nn,distance[ix],aa)
            dist_ok=np.ptp(ar)>1e-12
        Nl=int(nn.sum());Vl=int(vv.sum());p=Vl/Nl
        groups.append(dict(label=name,ix=ix,n=nn,v=vv,a=aa,N=Nl,V=Vl,W=Nl-Vl,
                           p=p,info=info,distance_ok=dist_ok,a_residual=ar,
                           D=Nl*p*(1-p) if 0<p<1 else 0.0))
    return groups


def primary_support(groups):
    usable=[g for g in groups if g['info']]
    N=sum(g['N'] for g in usable)
    V=sum(g['V'] for g in usable)
    return dict(n_lines=len(usable),N=N,V=V,W=N-V,
                status='ELIGIBLE' if len(usable)>=2 and N>=30 and V>=5 and N-V>=5 else 'HOLD_LINE_OR_EXPOSURE_SUPPORT')


def distance_support(groups):
    usable=[g for g in groups if g['distance_ok']]
    N=sum(g['N'] for g in usable);V=sum(g['V'] for g in usable)
    return dict(n_lines=len(usable),N=N,V=V,W=N-V,
                status='ELIGIBLE' if len(usable)>=2 and N>=30 and V>=5 and N-V>=5 else 'HOLD_DISTANCE_RESIDUAL_SUPPORT')


def clumping_support(groups):
    usable=[g for g in groups if len(g['n'])>=2 and g['V']>0 and g['W']>0]
    return dict(n_lines=len(usable),status='ELIGIBLE' if len(usable)>=2 else 'HOLD_WITHINLINE_COLOUR_SUPPORT')


def plusone(observed, null,upper=True):
    null=np.asarray(null,dtype=float)
    return float((1+np.count_nonzero(null>=observed-1e-12 if upper else null<=observed+1e-12))/(len(null)+1))


def holm_two(p1,p2):
    p1=float(p1);p2=float(p2)
    if p1<=p2:return (min(1,2*p1),max(p1,p2))
    return (max(p1,p2),min(1,2*p2))


def site_conditioned(n,v,a,line,R=R,seed=None,min_line_n=10,distance=None):
    n=np.asarray(n,dtype=int);v=np.asarray(v,dtype=int);a=np.asarray(a,dtype=float)
    B,between,E=B_and_line_null(n,v,line)
    groups=_stratify(n,v,a,line,distance=distance,min_line_n=min_line_n)
    rng=np.random.default_rng(seed)
    N=int(n.sum());V=int(v.sum());p=V/N;den=N*p*(1-p)
    Bsim=np.zeros(R,dtype=float); Tsim=np.zeros(R,dtype=float);Dsum=0.0;Tnum=0.0
    Tdsim=np.zeros(R,dtype=float);Dd=0.0;Tdnum=0.0
    for g in groups:
        nn=g['n'];vv=g['v'];Nl=g['N'];Vl=g['V'];pl=g['p']
        if 0<Vl<Nl and len(nn)>=2:
            x=rng.multivariate_hypergeometric(nn,Vl,size=R,method='marginals')
            if (x<0).any() or (x>nn).any() or not (x.sum(axis=1)==Vl).all():
                raise AssertionError('COLOUR_AND_PLOT_INVENTORY_NOT_PRESERVED')
            Bsim+=np.sum((x-nn*p)**2/nn,axis=1)/den
        else:
            x=None
            Bsim+=np.sum((vv-nn*p)**2/nn)/den
        if g['info']:
            Dsum+=g['D']
            Tnum+=float(np.dot(vv-nn*pl,g['a']))
            if x is not None:
                Tsim+=np.dot(x-nn*pl,g['a'])
        if g['distance_ok']:
            Dd+=g['D']
            Tdnum+=float(np.dot(vv-nn*pl,g['a_residual']))
            if x is not None:
                Tdsim+=np.dot(x-nn*pl,g['a_residual'])
    out={'B_observed':B,'line_between':between,'line_null_EB':E,'delta_line':B-E,
         'B_null':Bsim,'n_replicates':R,
         'T_observed':Tnum/Dsum if Dsum else None,'T_null':Tsim/Dsum if Dsum else None,
         'T_distance_observed':Tdnum/Dd if Dd else None,'T_distance_null':Tdsim/Dd if Dd else None,
         'primary_support':primary_support(groups),'distance_support':distance_support(groups),
         'clumping_support':clumping_support(groups),'line_groups':groups}
    return out


def _source_summary(d):
    return dict(n_field_quadrat_records=len(d),M_v=int(d.v.sum()),M_w=int(d.w.sum()),
                other_L=int(d.l.sum()),other_S=int(d.s.sum()),n_zero_target_quadrats=int((d.n==0).sum()),
                n_source_siteyears=int(d.groupby(['year','population']).ngroups),
                original_tab_sha256=FREQUENCY_SHA256)


def _evaluate_records(d,mode='all',R=R):
    rows=[]
    for (year,site),g in d.groupby(['year','population'],sort=True):
        N=int(g.n.sum());V=int(g.v.sum());W=N-V;m=int((g.n>0).sum())
        status= 'ELIGIBLE' if N>=20 and V>=5 and W>=5 and m>=2 else ('HOLD_TOO_FEW_M' if N<20 else 'HOLD_SITE_MONO_OR_COLOR_SUPPORT' if min(V,W)<5 else 'HOLD_ONE_PLOT')
        rows.append(dict(year=year,site=site,siteyear=f'{year}|{site}',N=N,V=V,W=W,m=m,source_quadrats=len(g),
                         other_species_n=int(g.l.sum()+g.s.sum()),status_site=status,df=g))
    if mode=='latest':
        latest={}
        for r in rows:
            if r['site'] not in latest or r['year']>latest[r['site']]['year']:
                latest[r['site']]=r
        rows=sorted(latest.values(),key=lambda x:x['siteyear'])
    return rows


def analysis(d, R=R,latest=False,save_null=False):
    context=_evaluate_records(d,'latest' if latest else 'all',R)
    rngs={k:np.random.default_rng(SEEDS[k if not latest else 'latest_'+k]) for k in ['exposure','clumping','distance']}
    records=[];draws={'exposure':[],'clumping':[],'distance':[]}
    for s in context:
        row={k:v for k,v in s.items() if k!='df'}
        if s['status_site']!='ELIGIBLE':
            row.update(primary_status='HOLD_SITE',clumping_status='HOLD_SITE',distance_status='HOLD_SITE')
            records.append(row);continue
        g=s['df']
        n=g.n.to_numpy(dtype=int);v=g.v.to_numpy(dtype=int);a=g.a.to_numpy(dtype=float)
        l=g.line.to_numpy();distance=g.distance.to_numpy(dtype=float)
        groups=_stratify(n,v,a,l,distance=distance,min_line_n=10)
        psp=primary_support(groups);csp=clumping_support(groups);dsp=distance_support(groups)
        B,Bb,E=B_and_line_null(n,v,l)
        row.update(B_original=B,B_between_recorded_lines=Bb,B_conditioned_E=E,
                   residual_delta=B-E,primary_support_nlines=psp['n_lines'],
                   primary_support_N=psp['N'],primary_status=psp['status'],
                   distance_support_nlines=dsp['n_lines'],distance_support_N=dsp['N'],
                   distance_status=dsp['status'],clumping_support_nlines=csp['n_lines'],
                   clumping_status=csp['status'])
        # Source-matched line shuffles, but separate fixed seeds for all preplanned comparisons.
        for kind in draws:
            if kind=='exposure' and psp['status']!='ELIGIBLE':continue
            if kind=='clumping' and csp['status']!='ELIGIBLE':continue
            if kind=='distance' and dsp['status']!='ELIGIBLE':continue
            childseed=int(rngs[kind].integers(1,2**32-1))
            out=site_conditioned(n,v,a,l,R=R,seed=childseed,distance=distance,min_line_n=10)
            if kind=='exposure':
                row['T_congener_purple_minus_white']=out['T_observed']
                draws[kind].append((row['siteyear'],out['T_null']))
            elif kind=='clumping':
                draws[kind].append((row['siteyear'],out['B_null']-out['line_null_EB']))
            else:
                row['T_congener_after_linear_distance']=out['T_distance_observed']
                draws[kind].append((row['siteyear'],out['T_distance_null']))
        records.append(row)
    sums={}
    for kind in draws:
        sitecol={'exposure':'T_congener_purple_minus_white','clumping':'residual_delta',
                 'distance':'T_congener_after_linear_distance'}[kind]
        statcol={'exposure':'primary_status','clumping':'clumping_status','distance':'distance_status'}[kind]
        valid=[r for r in records if r.get(statcol)=='ELIGIBLE']
        names=[r['siteyear'] for r in valid]
        if len(valid)<3:
            sums[kind]={'status':'HOLD_LT3_SITEYEARS','n':len(valid),'names':names}
            continue
        array_by_name=dict(draws[kind]);
        assert set(array_by_name)==set(names)
        x=np.array([r[sitecol] for r in valid])
        sim=np.stack([array_by_name[n] for n in names])
        null=sim.mean(axis=0)
        obs=float(x.mean());p=plusone(obs,null)
        summary={'status':'ELIGIBLE','n':len(valid),'n_named_sites':len(set(r['site'] for r in valid)),
                 'names':names,'estimate':obs,'null_mean':float(np.mean(null)),
                 'p_one_sided':p,'n_positive_siteyears':int(np.count_nonzero(x>0)),
                 'siteyear_estimates':{r['siteyear']:float(r[sitecol]) for r in valid},
                 'null_CI95':[float(v) for v in np.quantile(null,[.025,.975])]}
        if kind=='clumping':
            summary['mean_B_original']=float(np.mean([r['B_original'] for r in valid]))
            summary['mean_between_recorded_lines']=float(np.mean([r['B_between_recorded_lines'] for r in valid]))
            summary['mean_expected_within_lines']=float(np.mean([r['B_conditioned_E']-r['B_between_recorded_lines'] for r in valid]))
            summary['n_effective_siteyears']=len(valid)
        # Posthoc (not prespecified inferential) one-siteyear-out fragility ledger.
        loo=[]
        if len(valid)>=4:
            for i, r in enumerate(valid):
                rest=np.delete(x,i);rnull=np.delete(sim,i,axis=0).mean(axis=0)
                loo.append({'dropped_siteyear':r['siteyear'],'estimate':float(rest.mean()),'p_one_sided':plusone(float(rest.mean()),rnull)})
        summary['posthoc_leave_one_siteyear_out']=loo
        sums[kind]=summary
        if save_null:summary['_null_replicates']=null.tolist()
    if sums['exposure']['status']=='ELIGIBLE' and sums['clumping']['status']=='ELIGIBLE':
        pa,pc=holm_two(sums['exposure']['p_one_sided'],sums['clumping']['p_one_sided'])
        sums['exposure']['holm_two_primary_family']=pa
        sums['clumping']['holm_two_primary_family']=pc
    return {'results':sums,'siteyear_ledger':records,'source':_source_summary(d),
            'analysis_scope':'retrospective_source_exposed_fixed_original_sites_line_conditioned_no_causal_selection'}


def _roundobj(x):
    if isinstance(x,float):return round(x,14)
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:_roundobj(v) for k,v in x.items()}
    if isinstance(x,list):return [_roundobj(v) for v in x]
    return x


def main():
    a=argparse.ArgumentParser();a.add_argument('--source',default='source/frequency.txt');a.add_argument('--outdir',default='results');a.add_argument('--reps',type=int,default=R)
    args=a.parse_args(); d=load_source(Path(args.source));outdir=Path(args.outdir);outdir.mkdir(exist_ok=True,parents=True)
    main=analysis(d,R=args.reps);latest=analysis(d,R=args.reps,latest=True)
    result={'schema':'fcp_sisyrinchium_line_conditional_cogeneric_exposure_20261013_v1',
            'protocol':'docs/FCP_SISYRINCHIUM_WITHIN_LINE_COGENERIC_EXPOSURE_PROTOCOL_20261013.md',
            'R':args.reps,'seeds':SEEDS,'full':main,'latest_year_per_site_sensitivity':latest,
            'no_new_independent_plants':True,'observational_not_causal':True,
            'source_paper_already_reported_cross_site_congener_association':True}
    (outdir/'result.json').write_text(json.dumps(_roundobj(result),indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8')
    rows=pd.DataFrame([{k:v for k,v in x.items()} for x in main['siteyear_ledger']])
    rows.to_csv(outdir/'siteyear_support_and_effects.csv',index=False,float_format='%.16g')
    for kind,summary in main['results'].items():
        print(json.dumps({'endpoint':kind,'status':summary['status'],'eligible':summary['n'],'estimate':summary.get('estimate'),
                          'p':summary.get('p_one_sided'),'holm':summary.get('holm_two_primary_family'),
                          'positive':summary.get('n_positive_siteyears'),'sites':summary.get('names')},ensure_ascii=False))
    print(json.dumps({'latest':{k:{'n':v['n'],'estimate':v.get('estimate'),'p':v.get('p_one_sided')} for k,v in latest['results'].items()}},ensure_ascii=False))

if __name__=='__main__':main()
