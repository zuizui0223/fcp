#!/usr/bin/env python3
"""Independent Fraction/Decimal check; never imports hierarchy.py.

Checks arithmetic on the same source transcription, not independent field data.
"""
import csv, itertools, json, math
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=Decimal
COLORS=('white','pink','purple')
checks=0;maxerr=0.

def close(a,b,tol=2e-12):
    global checks,maxerr
    e=abs(float(a)-float(b));maxerr=max(maxerr,e);checks+=1
    assert e<=tol,(a,b,e)

def truth(x):
    global checks
    checks+=1; assert x

def d(v):
    v=F(v);return D(v.numerator)/D(v.denominator)

def H(p):
    return -sum((d(x)*d(x).ln() for x in p if x),D(0))

def corner_points(l,u):
    pts=set()
    for fixed in itertools.combinations(range(3),2):
        f=next(i for i in range(3) if i not in fixed)
        for v in itertools.product(*[(l[i],u[i]) for i in fixed]):
            q=[F(0)]*3
            for i,x in zip(fixed,v):q[i]=x
            q[f]=1-sum(q)
            if all(l[i]<=q[i]<=u[i] for i in range(3)):pts.add(tuple(q))
    return sorted(pts)

def extremes(l,u):
    vertices=corner_points(l,u)
    # Independent finite KKT enumeration: each coordinate is lower/free/upper.
    candidates=list(vertices)
    for active in itertools.product((-1,0,1),repeat=3):
        free=[i for i,a in enumerate(active) if a==0]
        if not free:continue
        q=[l[i] if a==-1 else u[i] if a==1 else F(0) for i,a in enumerate(active)]
        level=(1-sum(q))/len(free)
        for i in free:q[i]=level
        if all(l[i]<=q[i]<=u[i] for i in range(3)):candidates.append(tuple(q))
    return min(H(q) for q in vertices),max(H(q) for q in candidates)

def quantities(q,n,np,mono,a):
    w=F(np,n);mu=[w*q[i]+F(mono[i],n) for i in range(3)]
    gamma=H(mu);hp=H(q);alpha=d(a)
    bm=d(F(sum(mono),n))*H([F(v,sum(mono)) for v in mono]) if sum(mono) else D(0)
    beta=gamma-alpha;bp=d(w)*hp-alpha;bs=gamma-d(w)*hp-bm
    return dict(gamma=gamma,alpha=alpha,beta=beta,poly=bp,mono=bm,status=bs,poly_share=bp/beta)

def check_case(row,saved,eps):
    n=int(row['n_all']);np=int(row['n_poly']);m=[int(row['mono_'+k]) for k in COLORS];w=F(np,n)
    q0=[F(row['poly_'+k]) for k in COLORS];mu=[F(row['all_'+k]) for k in COLORS]
    l=[max(F(0),q0[i]-eps,(mu[i]-eps-F(m[i],n))/w) for i in range(3)]
    u=[min(F(1),q0[i]+eps,(mu[i]+eps-F(m[i],n))/w) for i in range(3)]
    v=corner_points(l,u)
    if not v:
        truth(saved['status']=='HOLD_INCONSISTENT_PUBLISHED_SUMMARIES');return
    truth(saved['status']=='OK');truth(len(v)==len(saved['feasible_vertices']))
    q=[sum(x[i] for x in v)/len(v) for i in range(3)]
    for a,b in zip(q,saved['representative_poly_mean']):close(d(a),b)
    a=F(row['mean_H_all']);al=a-F(5,1000);ah=a+F(5,1000)
    point=quantities(q,n,np,m,a)
    for k,x in point.items():close(x,saved['representative'][k])
    close(point['beta'],point['poly']+point['mono']+point['status'])
    hp=extremes(l,u);hg=extremes([w*x+F(c,n) for x,c in zip(l,m)],[w*x+F(c,n) for x,c in zip(u,m)])
    b=(hg[0]-d(ah),hg[1]-d(al));p=(d(w)*hp[0]-d(ah),d(w)*hp[1]-d(al));bm=point['mono']
    expected=dict(gamma=hg,alpha=(d(al),d(ah)),beta=b,poly=p,mono=(bm,bm),
                  status=(max(D(0),hg[0]-d(w)*hp[1]-bm),hg[1]-d(w)*hp[0]-bm),
                  poly_share=(p[0]/b[1],p[1]/b[0]))
    for k,pair in expected.items():
        for x,y in zip(pair,saved['enclosures'][k]):close(x,y)
    # All vertices + deterministic segments between every pair of vertices.
    for x,y in itertools.combinations_with_replacement(v,2):
        for j in range(11):
            qq=[(x[k]*(10-j)+y[k]*j)/10 for k in range(3)]
            for aa in (al,a,ah):
                stats=quantities(qq,n,np,m,aa)
                for k,val in stats.items():
                    low,high=saved['enclosures'][k]
                    truth(float(low)-2e-12<=float(val)<=float(high)+2e-12)
    truth(saved['enclosures']['poly_share'][0]>.5)

with localcontext() as context:
    context.prec=60
    with (ROOT/'data/conditional_summary.csv').open() as f: rows=list(csv.DictReader(f))
    result=json.loads((ROOT/'results/result.json').read_text())
    for row,saved in zip(rows,result['cases']):check_case(row,saved,F(5,10000))
    for row,saved in zip(rows,result['discrepancy_sensitivity']):check_case(row,saved,F(1,1000))
    # exact 2022 inconsistency and minimal common half-width necessary for white
    n=112;np=102;gap=F('0.1735')-F(np,n)*F('0.1905')
    truth(gap==F(1,112000));min_half=(F('.174')-F(np,n)*F('.190'))/(1+F(np,n))
    out={'status':'PASS','checks':checks,'maximum_absolute_difference':maxerr,'decimal_digits':60,
         'same_transcription_not_independent_source_reextraction':True,
         'strict_2022_gap_all_mean':str(d(gap)),'strict_2022_gap_fraction':str(gap),
         'minimum_common_frequency_halfwidth_required_by_white_2022':str(d(min_half)),
         'checks_include_independent_KKT_entropy_extrema':True,'outer_enclosures_not_sampling_CI':True}
    (ROOT/'results/independent_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))
