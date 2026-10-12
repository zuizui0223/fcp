#!/usr/bin/env python3
"""Observation-subsampling-only null; not drift, adaptation or a corrected effect size.

Agrawal (2020) MGF + sampling-without-replacement convex comparison + Jensen.
No independence between stands is needed. Incomplete stands must be simple
random samples of at least n0; complete finite censuses have zero sampling error.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from parent_entropy import entropy, summary_result

INPUT_SHA='5c1a2fa4b0a7b6cf7f79f1fda10e8e253036af1b0af725fdeced85d0783cd5be'
SCENARIOS=(1,25,50,75,100)
MARGIN=1e-12


def positive_integer(x: int) -> bool:
    return isinstance(x,int) and not isinstance(x,bool) and x>0


def tail_upper(gap_nats: float, min_subsample_n: int, k: int=3) -> float:
    """Conservative analytic p-value upper bound, NOT an exact null probability."""
    if not math.isfinite(gap_nats) or gap_nats<0:
        raise ValueError('gap must be nonnegative and finite')
    if not positive_integer(min_subsample_n) or not positive_integer(k) or k<2:
        raise ValueError('positive integer n and integer k>=2 required')
    x=min_subsample_n*gap_nats; a=k-1
    if x<=a: return 1.0
    return min(1.0, math.exp(-x+a*(1+math.log(x/a)))*(1+1e-12))


def minimum_n(gap_nats: float, k: int=3, alpha: float=.05/3) -> int:
    if not math.isfinite(gap_nats) or gap_nats<=0 or not 0<alpha<1:
        raise ValueError('positive finite gap and alpha in (0,1) required')
    hi=1
    while tail_upper(gap_nats,hi,k)>alpha:
        hi*=2
        if hi>10**8: raise ValueError('threshold beyond supported range')
    lo=0
    while hi-lo>1:
        mid=(lo+hi)//2
        if tail_upper(gap_nats,mid,k)<=alpha:hi=mid
        else:lo=mid
    return hi


def validate_design(sample_n: int, census_N: int, n0: int=100) -> str:
    if not all(positive_integer(x) for x in (sample_n,census_N,n0)) or sample_n>census_N:
        raise ValueError('invalid sample/census size')
    if sample_n==census_N:return 'census'
    if sample_n<n0:raise ValueError('incomplete stand below assumed sampling minimum')
    return 'subsample'


def kl(q, p) -> float:
    entropy(q);entropy(p)
    if len(q)!=len(p):raise ValueError('different alphabets')
    if any(a>0 and b==0 for a,b in zip(q,p)):return math.inf
    return math.fsum(a*math.log(a/b) for a,b in zip(q,p) if a>0)


def gap(rows, weights) -> float:
    if not rows or len(rows)!=len(weights):raise ValueError('incompatible rows and weights')
    entropy(weights)
    if len({len(r) for r in rows})!=1:raise ValueError('different alphabets')
    means=[math.fsum(w*r[j] for w,r in zip(weights,rows)) for j in range(len(rows[0]))]
    val=entropy(means)-math.fsum(w*entropy(r) for w,r in zip(weights,rows))
    if val<-1e-10:raise ValueError('negative Jensen gap')
    return max(0.,val)


def compositions(n: int, k: int):
    if k==1:
        yield (n,);return
    for first in range(n+1):
        for rest in compositions(n-first,k-1):yield (first,)+rest


def multinomial(n: int, p):
    out=[]
    for c in compositions(n,len(p)):
        if any(v and prob==0 for v,prob in zip(c,p)):continue
        prob=math.factorial(n)/math.prod(math.factorial(v) for v in c)
        prob*=math.prod(a**v for a,v in zip(p,c))
        out.append((tuple(v/n for v in c),prob))
    return out


def hypergeometric(n: int, counts):
    den=math.comb(sum(counts),n);out=[]
    for c in compositions(n,len(counts)):
        if any(v>N for v,N in zip(c,counts)):continue
        prob=math.prod(math.comb(N,v) for N,v in zip(counts,c))/den
        out.append((tuple(v/n for v in c),prob))
    return out


def exact_checks() -> dict:
    """Enumerated finite controls; they supplement rather than prove the theorem."""
    checks=0; max_mass_error=0.
    def verify(cond):
        nonlocal checks
        checks+=1
        if not cond:raise AssertionError('exact finite-control violation')
    for p in ((.2,.3,.5),(.01,.09,.9),(1/3,1/3,1/3)):
        for n in range(1,13):
            dist=multinomial(n,p)
            err=abs(sum(pr for _,pr in dist)-1);max_mass_error=max(max_mass_error,err)
            verify(err<1e-12)
            for f in (0.,.2,.5,.8):
                mgf=math.fsum(pr*math.exp(f*n*kl(q,p)) for q,pr in dist)
                verify(mgf<=(1-f)**-2+1e-11)
    for counts in ((4,3,5),(0,2,4),(6,6,0)):
        N=sum(counts);p=tuple(v/N for v in counts)
        for n in range(1,N+1):
            a=hypergeometric(n,counts);b=multinomial(n,p)
            verify(abs(sum(pr for _,pr in a)-1)<1e-12)
            for f in (.2,.5,.8):
                wo=math.fsum(pr*math.exp(f*n*kl(q,p)) for q,pr in a)
                wr=math.fsum(pr*math.exp(f*n*kl(q,p)) for q,pr in b)
                verify(wo<=wr+1e-11)
                verify(wo<=(1-f)**-2+1e-11)
            if n==N:verify(len(a)==1 and abs(kl(a[0][0],p))<1e-12)
    # Perfectly dependent, opposite deviations across two binary stands.
    for n in range(1,21):
        worlds=[(gap([q,tuple(reversed(q))],[.5,.5]),pr) for q,pr in multinomial(n,(.5,.5))]
        for cut in (.01,.05,.1,.2,.4,.65):
            tail=math.fsum(pr for val,pr in worlds if val>=cut)
            verify(tail<=tail_upper(cut,n,2)+1e-11)
        # One stand is an exactly observed census with true 50:50 composition.
        for q,pr in multinomial(n,(.5,.5)):
            verify(gap([q,(.5,.5)],[.5,.5])<=.5*kl(q,(.5,.5))+1e-12)
    return {'status':'PASS','checks':checks,'max_multinomial_mass_error':max_mass_error,
            'coverage':['multinomial_MGF','without_replacement_convex_order','perfect_dependence','census_mixture'],
            'empirical_enumerations_not_a_substitute_for_the_proof':True}


def run(root: Path) -> dict:
    src=root/'data/hesperis_published_summary.csv'
    if hashlib.sha256(src.read_bytes()).hexdigest()!=INPUT_SHA:raise ValueError('parent source hash mismatch')
    with src.open() as f:rows=list(csv.DictReader(f))
    out=[]
    for row in rows:
        parent=summary_result(row)
        b=max(0.,parent['rounding_bounds']['H_between'][0]-MARGIN)
        cases=[]
        for n0 in SCENARIOS:
            u=tail_upper(b,n0,3)
            cases.append({'min_incomplete_stand_sample_n':n0,'tail_upper':u,
                          'bonferroni_three_years_upper':min(1.,3*u),
                          'expected_null_gap_upper_nats':math.log1p(2/n0)})
        out.append({'year':parent['year'],'n_stands':parent['n_stands'],
                    'observed_gap_nats':parent['H_between'],
                    'rounding_lower_gap_nats':b,
                    'observed_between_fraction_UNCORRECTED':parent['between_fraction'],
                    'minimum_n0_for_three_year_bonferroni':minimum_n(b),
                    'scenarios':cases})
    result={'schema':'fcp_hesperis_observation_subsampling_bound_v1',
            'status':'CONDITIONAL_SUMMARY_BASED_BOUND_COMPLETE_RAW_COUNTS_UNRECOVERED',
            'parent_commit':'2128286fc39a1da7207f9b8f294af698fc3de0f8',
            'input_sha256':INPUT_SHA,'k':3,'logarithm':'natural','numerical_gap_margin':MARGIN,
            'primary_min_subsample_n':100,'requires_reported_design_to_hold':True,
            'individual_sampling_compliance_audited':False,
            'between_stand_independence_required':False,
            'between_year_independence_required':False,
            'null':'identical finite census colour compositions; observation subsampling only',
            'not_null':'a common superpopulation with demographic drift or finite population formation',
            'effect_size_bias_correction_performed':False,
            'no_reconstructed_raw_counts':True,'prior_results_changed':False,
            'years':out,'exact_controls':exact_checks()}
    dest=root/'results';dest.mkdir(exist_ok=True)
    (dest/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    r=run(ap.parse_args().root)
    print(json.dumps({'years':r['years'],'exact_controls':r['exact_controls']},indent=2))
