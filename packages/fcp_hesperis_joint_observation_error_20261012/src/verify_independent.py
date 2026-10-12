#!/usr/bin/env python3
"""Independent 85-digit Decimal replay; does NOT import joint_bound.py.
Verifies source identity, tail formula, Fano variation, inverse thresholds and
frozen JSON/CSV. This is independent arithmetic, not a new proof/field sample.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from decimal import Decimal,localcontext
from pathlib import Path

BASE=Path(__file__).resolve().parents[1]
PR173='980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1'
D=lambda x: Decimal(str(x))


def min1(value):return min(D(1),value)


def info_loss(e):
    if e==0:return D(0)
    if e>=D(2)/D(3):return D(3).ln()
    return -e*e.ln()-(1-e)*(1-e).ln()+e*D(2).ln()


def bound_raw(b,n,e,eta=D(0)):
    rem=b-info_loss(e)
    if rem<=0 or n*rem<=2: u=D(1)
    else: u=min1((-n*rem + 2 + 2*(n*rem/D(2)).ln()).exp())
    return min1(u+eta),rem


def upper_three(b,n,e,eta=D(0)):
    p,_=bound_raw(b,D(n),e,eta)
    return min1(D(3)*p)


def critical(b,n,alpha=D('.05')):
    if upper_three(b,n,D(0))>alpha:return None
    low,high=D(0),D(2)/D(3)
    for _ in range(145):
        mid=(low+high)/2
        if upper_three(b,n,mid)<=alpha:low=mid
        else:high=mid
    return low


def close(a,b,atol=5e-13,rtol=6e-12):
    return abs(float(a)-float(b)) <= atol+rtol*abs(float(a))


def run(root:Path)->dict:
    checks=0;largest=0;largest_relative=0
    def eq(label,got,want,atol=5e-13,rtol=6e-12):
        nonlocal checks,largest,largest_relative
        if not close(got,want,atol,rtol):
            raise ValueError(f'{label}: Decimal {got} vs production {want}')
        diff=abs(float(got)-float(want));largest=max(largest,diff)
        if float(got)>1e-13:largest_relative=max(largest_relative,diff/float(got))
        checks+=1
    root=Path(root)
    src=(root/'data'/'parent_entropy_result.json').read_bytes()
    if hashlib.sha256(src).hexdigest()!=PR173:raise ValueError('PR173 byte hash mismatch')
    parent=json.loads(src)
    obs=json.loads((root/'results'/'result.json').read_text())
    if obs['parent_PR173_sha256']!=PR173:raise ValueError('incorrect result source')
    if obs['measured_misclassification_rate'] is not None or obs['raw_stand_counts_recovered']:
        raise ValueError('biological observation overclaim')
    if obs['n_years']!=3 or len(parent['hesperis'])!=3:raise ValueError('different annual opportunity')
    if obs['classification_error_semantics'].find('pathwise')<0:raise ValueError('classification assumption lost')
    with localcontext() as ctx:
        ctx.prec=85
        for r in parent['hesperis']:
            yr=r['year']; b=D(r['rounding_bounds']['H_between'][0]); t=obs['annual'][str(yr)]
            eq('B_lower',b,t['rounding_conservative_between_lower_nats'],atol=0)
            for n in [25,50,75,100]:
                cr=critical(b,n)
                actual=t['critical_error_n'+str(n)]
                if cr is None:
                    if actual is not None: raise ValueError('critical 0% should HOLD')
                    checks+=1
                else:
                    eq('critical_eps',cr,actual,atol=2e-14)
                for ep,s in t['scenario'][str(n)].items():
                    e=D(ep);up,rem=bound_raw(b,D(n),e)
                    u3=upper_three(b,n,e)
                    eq('Fano',info_loss(e),s['fano_bound_nats'])
                    eq('reduced_gap',rem,s['test_gap_after_classification_budget_nats'])
                    eq('adjusted_upper',u3,s['adjusted_bound'])
                    eq('raw_upper',up,s['unadjusted_joint_upper'])
                    eq('epsilon',e,s['epsilon_upper'])
                    eq('correct_sample_gap',max(D(0),rem),s['minimum_correctly_classified_sample_gap_nats'])
                    if (u3 <= D('.05')) != s['reject_joint_null_under_pathwise_budget_at_0_05']:
                        raise ValueError('decision mismatch')
                    checks+=1
        # Regression to separately frozen zero-error bound, rounded upward in PR174.
        frozen={2021:D('5.7051216422e-8'),2022:D('4.7620579400e-7'),2023:D('9.1683820904e-7')}
        for r in parent['hesperis']:
            y=r['year'];b=D(r['rounding_bounds']['H_between'][0]);p=upper_three(b,100,D(0))
            if abs(p-frozen[y]) > frozen[y]*D('.00005'):
                raise ValueError('cannot reproduce prior PR174 upper bound')
            checks+=1
        # Extra independent probability-budget unions, no independence multiplication.
        for eta in [D('.001'),D('.005'),D('.01')]:
            for r in parent['hesperis']:
                b=D(r['rounding_bounds']['H_between'][0]);e=D('.01')
                corrected=upper_three(b,100,e,eta)
                if corrected < upper_three(b,100,e,D(0)):
                    raise ValueError('probabilistic error-budget failure must not help')
                checks+=1
    csvpath=root/'results'/'joint_sensitivity.csv'
    with csvpath.open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=3*4*10:raise ValueError('CSV missing sensitivity rows')
    for row in rows:
        yr=int(row['year']);n=int(row['min_sample_n0']);e=row['epsilon_upper']
        val=obs['annual'][str(yr)]['scenario'][str(n)][e]
        eq('CSV joint upper',val['adjusted_bound'],row['p_upper_three_year_Bonferroni'])
    return {'status':'PASS','checks':checks,'decimal_digits':85,
            'max_absolute_difference':largest,'max_relative_difference':largest_relative,
            'source_sha256':PR173,
            'not_independent_publisher_data_or_measured_error':True,
            'proof_kind':'conditional_event_inclusion_fano_plus_mgf_and_union',
            'n_annual':3,'n_scenario_rows':len(rows)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=BASE)
    args=ap.parse_args()
    record=run(args.root)
    (args.root/'results'/'independent_verification.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps(record))

if __name__=='__main__':main()
