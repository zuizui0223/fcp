#!/usr/bin/env python3
"""Retrospective Shannon chain-rule partition using published group summaries.

No population records are reconstructed. Bounds are conservative deterministic
rounding enclosures, not sampling confidence intervals or tests of selection.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Sequence

COLORS=('white','pink','purple')
EPS=1e-11


def entropy(p: Sequence[float]) -> float:
    x=[float(v) for v in p]
    if not x or any(not math.isfinite(v) or v < -EPS or v > 1+EPS for v in x):
        raise ValueError('invalid probability vector')
    if abs(sum(x)-1)>EPS:
        raise ValueError('probabilities do not sum to one')
    return -math.fsum(v*math.log(v) for v in x if v>0)


def vertices(lo: Sequence[float], hi: Sequence[float]) -> list[list[float]]:
    """Vertices of box intersect sum(p)=1; dimension three for this source."""
    lo=list(map(float,lo));hi=list(map(float,hi));d=len(lo)
    if d!=3 or len(hi)!=d or any(not math.isfinite(v) for v in lo+hi):
        raise ValueError('three finite dimensions required')
    if any(l<0 or u>1 or l>u+EPS for l,u in zip(lo,hi)) or sum(lo)>1+EPS or sum(hi)<1-EPS:
        raise ValueError('infeasible probability bounds')
    out=[]
    for free in range(d):
        fixed=[i for i in range(d) if i!=free]
        for options in itertools.product((0,1),repeat=d-1):
            p=[0.]*d
            for j,o in zip(fixed,options):p[j]=(lo[j],hi[j])[o]
            p[free]=1-math.fsum(p)
            if all(lo[i]-EPS<=p[i]<=hi[i]+EPS for i in range(d)):
                p=[max(0.,min(1.,v)) for v in p]
                if not any(max(abs(a-b) for a,b in zip(p,v))<1e-12 for v in out):out.append(p)
    if not out:raise ValueError('empty simplex intersection')
    return sorted(out)


def entropy_bounds(lo: Sequence[float], hi: Sequence[float]) -> list[float]:
    """Concave entropy min at a vertex; max by bounded equalization (KKT)."""
    v=vertices(lo,hi)
    left,right=0.,1.
    for _ in range(90):
        t=(left+right)/2
        if sum(max(l,min(u,t)) for l,u in zip(lo,hi))>1:right=t
        else:left=t
    mx=[max(l,min(u,(left+right)/2)) for l,u in zip(lo,hi)]
    return [min(entropy(p) for p in v),entropy(mx)]


def _integer(v: float, name: str) -> int:
    x=float(v)
    if not math.isfinite(x) or x<0 or not x.is_integer():raise ValueError(name+' must be a nonnegative integer')
    return int(x)


def decompose(n: int, npoly: int, mono: Sequence[int], q: Sequence[float], alpha: float) -> dict:
    n=_integer(n,'n');npoly=_integer(npoly,'npoly');mono=[_integer(v,'mono') for v in mono]
    if not n or not npoly or len(mono)!=3 or npoly+sum(mono)!=n:raise ValueError('inconsistent stand counts')
    w=npoly/n;hp=entropy(q);alpha=float(alpha)
    if not math.isfinite(alpha) or alpha < -EPS or alpha>w*hp+EPS:
        raise ValueError('mean within entropy incompatible with group summaries')
    mu=[w*float(p)+c/n for p,c in zip(q,mono)]
    gamma=entropy(mu)
    bm=sum(mono)/n*entropy([v/sum(mono) for v in mono]) if sum(mono) else 0.
    bp=w*hp-alpha
    bs=gamma-w*hp-bm
    beta=gamma-alpha
    return dict(gamma=gamma,alpha=alpha,beta=beta,poly=bp,mono=bm,status=bs,
                poly_share=bp/beta if beta>1e-14 else None,mu=mu,poly_mean_entropy=alpha/w,
                chain_error=beta-(bp+bm+bs))


def evaluate(row: dict, freq_halfwidth: float = .0005) -> dict:
    if freq_halfwidth not in (.0005,.001):raise ValueError('undeclared frequency tolerance')
    n=_integer(row['n_all'],'n_all');np=_integer(row['n_poly'],'n_poly')
    mono=[_integer(row['mono_'+c],'mono') for c in COLORS]
    if np<=0 or np+sum(mono)!=n:raise ValueError('inconsistent group sizes')
    w=np/n;c=[v/n for v in mono]
    q0=[float(row['poly_'+k]) for k in COLORS]
    mu0=[float(row['all_'+k]) for k in COLORS]
    lo=[max(0.,q-freq_halfwidth,(mu-freq_halfwidth-cc)/w) for q,mu,cc in zip(q0,mu0,c)]
    hi=[min(1.,q+freq_halfwidth,(mu+freq_halfwidth-cc)/w) for q,mu,cc in zip(q0,mu0,c)]
    vv=vertices(lo,hi)
    rep=[sum(v[k] for v in vv)/len(vv) for k in range(3)]
    a=float(row['mean_H_all']);al=max(0.,a-.005);ah=a+.005
    hp=entropy_bounds(lo,hi)
    hg=entropy_bounds([w*p+cc for p,cc in zip(lo,c)],[w*p+cc for p,cc in zip(hi,c)])
    if ah>w*hp[0]+EPS:raise ValueError('additional within-entropy feasibility constraints needed')
    point=decompose(n,np,mono,rep,a)
    bm=point['mono'];b=[hg[0]-ah,hg[1]-al];p=[w*hp[0]-ah,w*hp[1]-al]
    if b[0]<=0:raise ValueError('beta denominator not bounded away from zero')
    # Independent extrema make these OUTER enclosures; no sharpness is asserted.
    # h appears in numerator and denominator, so plain interval division is conservative.
    enc={'gamma':hg,'alpha':[al,ah],'beta':b,'poly':p,'mono':[bm,bm],
         'status':[max(0.,hg[0]-w*hp[1]-bm),hg[1]-w*hp[0]-bm],
         'poly_share':[max(0.,p[0]/b[1]),min(1.,p[1]/b[0])]}
    return {'year':int(row['year']),'n_all':n,'n_poly':np,'n_mono':sum(mono),
            'literal_poly_sum':sum(q0),'literal_mixture_minus_all':
            [w*q+cc-mu for q,cc,mu in zip(q0,c,mu0)],
            'poly_lower':lo,'poly_upper':hi,'feasible_vertices':vv,
            'representative_poly_mean':rep,'representative':point,
            'enclosures':enc,'status':'OK','freq_halfwidth':freq_halfwidth,
            'interval_type':('CONSERVATIVE_ROUNDING_OUTER_ENCLOSURE_NOT_CI' if freq_halfwidth==.0005 else 'PUBLICATION_DISCREPANCY_SENSITIVITY_NOT_CI'),
            'majority_among_poly_throughout_outer_enclosure':enc['poly_share'][0]>.5}


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args()
    source=a.root/'data/conditional_summary.csv'
    with source.open() as f:rows=list(csv.DictReader(f))
    out={'schema':'fcp_hesperis_conditional_entropy_hierarchy_v1','date':'2026-10-12',
         'source_doi':'10.1093/aob/mcag035','parent_commit':'2128286fc39a1da7207f9b8f294af698fc3de0f8',
         'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'cases':[],'discrepancy_sensitivity':[evaluate(r,freq_halfwidth=.001) for r in rows],
         'source_status':'PUBLISHED_ROUNDED_SUMMARIES_ONLY_RAW_COUNTS_NOT_RECOVERED',
         'representative_not_recovered_true_mean':True,'raw_stand_records_reconstructed':False,
         'hypothesis_test_performed':False,'finite_sample_bias_corrected':False,
         'paired_year_trend_estimated':False,'parent_results_changed':False,
         'interpretation':'Descriptive chain-rule allocation: not selection, drift, genetics, or independent replication.'}
    for r in rows:
        try:out['cases'].append(evaluate(r))
        except ValueError as e:
            out['cases'].append({'year':int(r['year']),'status':'HOLD_INCONSISTENT_PUBLISHED_SUMMARIES','reason':str(e)})
    d=a.root/'results';d.mkdir(exist_ok=True)
    (d/'result.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    with (d/'hierarchy.csv').open('w',newline='') as f:
        keys=['analysis','year','n_all','n_poly','gamma','alpha','beta','poly','mono','status','poly_share','poly_share_lo','poly_share_hi']
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');w.writeheader()
        for label in ['cases','discrepancy_sensitivity']:
          for x in out[label]:
            if x['status']!='OK':continue
            w.writerow({'analysis':label,**{k:x[k] for k in ['year','n_all','n_poly']},
                        **{k:x['representative'][k] for k in ['gamma','alpha','beta','poly','mono','status','poly_share']},
                        'poly_share_lo':x['enclosures']['poly_share'][0],'poly_share_hi':x['enclosures']['poly_share'][1]})
    for x in out['cases']+out['discrepancy_sensitivity']:
        if x['status']!='OK':
            print(x);continue
        print(x['year'],'poly share:',x['representative']['poly_share'],'enclosure:',x['enclosures']['poly_share'],
              'components:',{k:x['representative'][k] for k in ['beta','poly','mono','status']})

if __name__=='__main__': main()
