#!/usr/bin/env python3
"""Independent 70-digit Decimal recreation; never imports information_bound."""
from __future__ import annotations
import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path

SOURCE_SHA = '980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1'

def fano(e,k=3):
    one=Decimal(1)
    ceiling=Decimal(k).ln()
    if e>=one-one/Decimal(k):
        return ceiling
    if e==0:
        return Decimal(0)
    return -e*e.ln()-(one-e)*(one-e).ln()+e*Decimal(k-1).ln()

def root(b):
    lo,hi=Decimal(0),Decimal(2)/Decimal(3)
    for _ in range(240):
        mid=(lo+hi)/2
        if fano(mid)<b:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2

def verify(parent,generated):
    with localcontext() as ctx:
        ctx.prec=75
        b=Path(parent).read_bytes()
        assert hashlib.sha256(b).hexdigest()==SOURCE_SHA, 'parent SHA mismatch'
        p=json.loads(b);g=json.loads(Path(generated).read_bytes())
        assert p['schema']=='fcp_hesperis_summary_shannon_partition_v1'
        assert g['schema']=='fcp_hesperis_error_information_sensitivity_v1'
        assert g['parent']['parent_sha256']==SOURCE_SHA
        checks=[]
        for src,got in zip(p['hesperis'],g['annual']):
            assert int(src['year'])==int(got['year'])
            lo=Decimal(str(src['rounding_bounds']['H_between'][0]))
            hi=Decimal(str(src['rounding_bounds']['H_between'][1]))
            center=Decimal(str(src['H_between']))
            Hgmax=Decimal(str(src['rounding_bounds']['H_total'][1]))
            assert lo <= center <= hi
            assert abs(root(lo)-Decimal(str(got['threshold_error_rate'])))<Decimal('1e-14')
            checks.append(('root',str(src['year'])))
            for case in got['scenarios']:
                e=Decimal(str(case['error_rate_upper']))
                phi=fano(e)
                bottom=max(Decimal(0),lo-phi)
                upper=min(Decimal(3).ln(),Hgmax+phi)
                fraction=bottom/upper
                for key,target in [('fano_information_ceiling_nats',phi),('guaranteed_true_sample_B_min',bottom),('guaranteed_true_sample_fraction_min',fraction)]:
                    assert abs(Decimal(str(case[key]))-target)<Decimal('1e-13'),(src['year'],e,key)
                    checks.append((str(src['year']),str(e),key))
        result={'status':'PASS','n_checks':len(checks),'precision':ctx.prec,
                'independent_implementation':True,'parent_sha256':SOURCE_SHA,
                'raw_stand_counts_recovered':False,'measured_colour_misclassification_rate':None}
        return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--parent',type=Path,required=True);a.add_argument('--result',type=Path,required=True);a.add_argument('--outdir',type=Path,required=True);v=a.parse_args()
    out=verify(v.parent,v.result)
    v.outdir.mkdir(parents=True,exist_ok=True)
    (v.outdir/'independent_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))
