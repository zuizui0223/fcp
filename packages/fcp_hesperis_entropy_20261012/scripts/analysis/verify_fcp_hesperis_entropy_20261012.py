#!/usr/bin/env python3
"""Independent 50-digit Decimal check. Does not import the analysis module."""
import argparse,csv,hashlib,itertools,json,math
from decimal import Decimal,localcontext
from pathlib import Path


def verify(root):
    root=Path(root)
    result=json.loads((root/'results/result.json').read_text())
    checks=0; maximum=0.
    def close(got,want):
        nonlocal checks,maximum
        delta=abs(float(got)-float(want)); checks+=1; maximum=max(maximum,delta)
        if not math.isfinite(delta) or delta>2e-12: raise ValueError((got,str(want),delta))
    def H(p): return -sum((x*x.ln() for x in p if x),Decimal(0))
    with localcontext() as ctx:
        ctx.prec=50
        with (root/'data/hesperis_published_summary.csv').open() as f: source=list(csv.DictReader(f))
        for raw,found in zip(source,result['hesperis'],strict=True):
            p=[Decimal(raw[k]) for k in ('white','pink','purple')]
            alpha=Decimal(raw['mean_shannon']); gamma=H(p); beta=gamma-alpha
            for key,want in [('H_total',gamma),('H_within',alpha),('H_between',beta),
                             ('between_fraction',beta/gamma),('effective_beta',beta.exp())]:
                close(found[key],want)
            # For these widely separated probability boxes, H is increasing
            # under transfer from largest to smallest. Both extrema are vertices.
            eps=Decimal('0.0005'); vertices=[]
            for free in range(3):
                fixed=[i for i in range(3) if i!=free]
                for signs in itertools.product((-1,1),repeat=2):
                    q=[Decimal(0)]*3
                    for i,s in zip(fixed,signs):q[i]=p[i]+s*eps
                    q[free]=1-sum(q)
                    if p[free]-eps<=q[free]<=p[free]+eps:vertices.append(H(q))
            glo,ghi=min(vertices),max(vertices)
            alo,ahi=alpha-Decimal('.005'),alpha+Decimal('.005')
            for key,values in [('H_total',[glo,ghi]),('H_within',[alo,ahi]),
                ('H_between',[glo-ahi,ghi-alo]),('between_fraction',[1-ahi/glo,1-alo/ghi]),
                ('effective_beta',[(glo-ahi).exp(),(ghi-alo).exp()])]:
                for got,want in zip(found['rounding_bounds'][key],values):close(got,want)
            # Independent feasible grid checks; reject non-simplex pseudo-data.
            for a in range(11):
                for b in range(11):
                    q=[p[0]-eps+Decimal(a)*eps/5,p[1]-eps+Decimal(b)*eps/5,Decimal(0)]
                    q[2]=1-q[0]-q[1]
                    if p[2]-eps<=q[2]<=p[2]+eps:
                        value=H(q); checks+=1
                        if not glo-Decimal('1e-45')<=value<=ghi+Decimal('1e-45'):
                            raise ValueError('independent grid outside entropy envelope')
        with (root/'data/linanthus_counts.csv').open() as f: counts=list(csv.DictReader(f))
        for out in result['linanthus_same_metric']:
            rows=[r for r in counts if int(r['year'])==out['year']]
            n=[Decimal(r['blue'])+Decimal(r['white']) for r in rows]
            checks+=1
            if int(sum(n))!=out['N']:raise ValueError('count discrepancy')
            if min(n)<20:
                if out['status']!='HOLD_NO_PLANTS' and out['status']!='HOLD_SPARSE_SECTIONS':
                    raise ValueError('sparse year should be HOLD')
                continue
            p=[[Decimal(r['blue'])/nn,Decimal(r['white'])/nn] for r,nn in zip(rows,n)]
            w=[Decimal(1)/4]*4 if out['weighting']=='equal_section' else [nn/sum(n) for nn in n]
            mean=[sum(ww*pp[k] for ww,pp in zip(w,p)) for k in range(2)]
            alpha=sum(ww*H(pp) for ww,pp in zip(w,p)); gamma=H(mean); beta=gamma-alpha
            # KL expression verifies that alpha subtraction is a true partition.
            kl=sum(ww*sum(pp[k]*(pp[k]/mean[k]).ln() for k in range(2) if pp[k]) for ww,pp in zip(w,p))
            for key,want in [('H_total',gamma),('H_within',alpha),('H_between',beta),
                ('between_fraction',beta/gamma),('effective_beta',beta.exp())]:close(out[key],want)
            close(out['H_between'],kl)
    receipt={'status':'PASS','checks':checks,'maximum_absolute_difference':maximum,
        'precision_decimal_digits':50,'same_source_transcription_not_independent_source_recovery':True,
        'rounding_envelopes_not_sampling_confidence_intervals':True,
        'result_sha256':hashlib.sha256((root/'results/result.json').read_bytes()).hexdigest()}
    (root/'results/independent_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt));return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    verify(p.parse_args().root)
