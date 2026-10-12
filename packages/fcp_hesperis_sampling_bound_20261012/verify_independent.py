#!/usr/bin/env python3
"""60-digit source-to-bound recomputation; imports no production analysis."""
import csv
import hashlib
import itertools
import json
from decimal import Decimal, localcontext
from pathlib import Path

D=Decimal

def verify(root):
    saved=json.loads((root/'results/result.json').read_text())
    raw=root/'data/hesperis_published_summary.csv'
    rows=list(csv.DictReader(raw.open()))
    checks=0;max_abs=0.;max_relative=0.
    def check(cond):
        nonlocal checks
        checks+=1
        if not cond:raise AssertionError('independent verification failed')
    def close(actual,expected):
        nonlocal max_abs,max_relative
        a=D(str(actual));err=abs(a-expected)
        max_abs=max(max_abs,float(err))
        if expected:max_relative=max(max_relative,float(err/abs(expected)))
        check(err<=max(D('1e-20'),abs(expected)*D('1e-10')))
    def H(p):return -sum(x*x.ln() for x in p if x)
    def U(b,n):
        x=D(n)*b
        return D(1) if x<=2 else (-x+D(2)*(D(1)+(x/2).ln())).exp()
    check(hashlib.sha256(raw.read_bytes()).hexdigest()==saved['input_sha256'])
    check(len(rows)==len(saved['years'])==3)
    with localcontext() as ctx:
        ctx.prec=60
        for row,s in zip(rows,saved['years']):
            check(int(row['year'])==s['year']);check(int(row['n_stands'])==s['n_stands'])
            check(row['n_stands']==row['alpha_n_stands'])
            p=[D(row[k]) for k in ('white','pink','purple')];alpha=D(row['mean_shannon'])
            close(s['observed_gap_nats'],H(p)-alpha)
            close(s['observed_between_fraction_UNCORRECTED'],1-alpha/H(p))
            lower=[x-D('.0005') for x in p];upper=[x+D('.0005') for x in p]
            hs=[]
            for i in range(3):
                other=[j for j in range(3) if j!=i]
                for ends in itertools.product((0,1),repeat=2):
                    q=[D(0)]*3
                    for j,e in zip(other,ends):q[j]=(lower,upper)[e][j]
                    q[i]=1-sum(q)
                    if lower[i]<=q[i]<=upper[i]:hs.append(H(q))
            b=min(hs)-alpha-D('.005')-D('1e-12')
            close(s['rounding_lower_gap_nats'],b)
            for c in s['scenarios']:
                n=c['min_incomplete_stand_sample_n'];u=U(b,n)
                close(c['tail_upper'],min(D(1),u*(1+D('1e-12'))))
                close(c['bonferroni_three_years_upper'],min(D(1),3*u*(1+D('1e-12'))))
                close(c['expected_null_gap_upper_nats'],(1+D(2)/n).ln())
                # Check upper bound is not rounded anti-conservatively.
                check(D(str(c['tail_upper']))>=u)
            ns=next(n for n in range(1,10001) if U(b,n)*(1+D('1e-12'))<=D('.05')/3)
            check(ns==s['minimum_n0_for_three_year_bonferroni'])
        check(saved['between_stand_independence_required'] is False)
        check(saved['effect_size_bias_correction_performed'] is False)
        check(saved['individual_sampling_compliance_audited'] is False)
    result={'status':'PASS','checks':checks,'decimal_precision':60,
            'max_abs_difference':max_abs,'max_relative_difference':max_relative,
            'production_code_imported':False,'raw_census_recovered':False,
            'scope':'input hash, three source/rounding partitions, all 15 bounds and thresholds; not independent field sampling validation'}
    (root/'results/independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':verify(Path(__file__).resolve().parent)
