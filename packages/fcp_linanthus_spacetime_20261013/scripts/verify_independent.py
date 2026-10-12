#!/usr/bin/env python3
"""Independent exact-rational audit; imports no production analysis code."""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from fractions import Fraction
from collections import defaultdict
from pathlib import Path

def check(path,results):
    if hashlib.sha256(path.read_bytes()).hexdigest()!='fa408c7ad6e521a8410252d903b024a057b1ee4f7e27406927e77fcc22bf4ac7':
        raise ValueError('Source mismatch')
    with path.open(newline='') as f:
        rows=list(csv.DictReader(f))
    vv={};nn={};counts={};sides={};section_names=set()
    for r in rows:
        y=int(r['year']);s=r['transect']+'_'+r['side'];section_names.add(s)
        n=int(r['white'])+int(r['blue'])
        counts[y]=counts.get(y,0)+n
        vv[(y,s)]=Fraction(int(r['white']),n) if n else None
        nn[(y,s)]=n
    assert len(vv)==28 and len(section_names)==4
    years=sorted([y for y in counts if min(nn[(y,s)] for s in section_names)>=20])
    years100=sorted([y for y in counts if min(nn[(y,s)] for s in section_names)>=100])
    if years!=[1992,1993,1995,1996,1998] or years100!=[1992,1993,1995,1998]: raise AssertionError('Eligibility drift')
    audit={'status':'PASS','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'independent_method':'exact fractions for raw frequencies, variance decomposition and left-out-year predictions',
      'checks':0,'maximum_absolute_difference':0.,'max_relative_difference':0.,
      'no_independent_biological_data':True,'all_count_free_yields_missing_not_zero':True}
    def eq(name,got,expected):
        if isinstance(got,Fraction):got=float(got)
        if isinstance(expected,Fraction):expected=float(expected)
        if not math.isclose(got,expected,abs_tol=1e-12,rel_tol=1e-12):
            raise AssertionError((name,got,expected))
        audit['checks']+=1;audit['maximum_absolute_difference']=max(audit['maximum_absolute_difference'],abs(got-expected))
        if expected: audit['max_relative_difference']=max(audit['max_relative_difference'],abs(got-expected)/abs(expected))
    assert results['count_by_year']=={str(y):n for y,n in counts.items()}
    eq('allplants',sum(counts.values()),results['n_counted_plant_records_all_years'])
    for threshold,k in [(20,'primary'),(100,'sensitivity_min100')]:
        ys=years if threshold==20 else years100
        sub=[vv[(y,s)] for y in ys for s in sorted(section_names)]
        mu=sum(sub,Fraction(0))/len(sub)
        sm={s:sum((vv[(y,s)] for y in ys),Fraction(0))/len(ys) for s in section_names}
        ym={y:sum((vv[(y,s)] for s in section_names),Fraction(0))/4 for y in ys}
        vt=sum(((p-mu)**2 for p in sub),Fraction(0))/len(sub)
        vs=sum(((p-mu)**2 for p in sm.values()),Fraction(0))/4
        vy=sum(((p-mu)**2 for p in ym.values()),Fraction(0))/len(ys)
        vr=sum(((vv[(y,s)]-sm[s]-ym[y]+mu)**2 for y in ys for s in section_names),Fraction(0))/(4*len(ys))
        assert vt==vs+vy+vr
        part=results[k]['variance_partition']
        for name,val in [('total_variance',vt),('spatial_variance',vs),('temporal_variance',vy),('interaction_variance',vr),('spatial_fraction',vs/vt),('temporal_fraction',vy/vt),('interaction_fraction',vr/vt)]:eq(f'{k}.{name}',val,part[name])
        for s,p in sm.items():eq(f'{k}.site.{s}',p,part['section_white_fraction_means'][s])
        for y,p in ym.items():eq(f'{k}.year.{y}',p,part['year_white_fraction_means'][str(y)])
        annual={x['year']:x for x in results[k]['annual_side_contrasts']}
        preds={(x['year'],x['section'],x['model']):x for x in results[k]['cv_predictions']}
        sumloss={m:Fraction(0) for m in ['GLOBAL','SIDE','SECTION']}
        for y in ys:
            tr=[vv[(t,s)] for t in ys if t!=y for s in sorted(section_names)]
            full=sum(tr,Fraction(0))/len(tr)
            sd={side:sum((vv[(t,s)] for t in ys if t!=y for s in section_names if s.endswith('_'+side)),Fraction(0))/(2*(len(ys)-1)) for side in ('east','west')}
            st={s:sum((vv[(t,s)] for t in ys if t!=y),Fraction(0))/(len(ys)-1) for s in section_names}
            for s in sorted(section_names):
                values={'GLOBAL':full,'SIDE':sd[s.split('_')[1]],'SECTION':st[s]}
                for mod,pr in values.items():
                    saved=preds[(y,s,mod)]
                    eq(f'{k}.predict.{y}.{s}.{mod}',pr,saved['prediction'])
                    sq=(vv[(y,s)]-pr)**2
                    eq(f'{k}.squared_error.{y}.{s}.{mod}',sq,saved['error_sq'])
                    eq(f'{k}.obs.{y}.{s}.{mod}',vv[(y,s)],saved['observed_white_fraction'])
                    sumloss[mod]+=sq
            for side in ('east','west'):
                twos=[s for s in section_names if s.endswith('_'+side)]
                p=sum((vv[(y,s)] for s in twos),Fraction(0))/2
                eq(f'{k}.annual.{y}.{side}',p,annual[y][side]['equal_transect_white_fraction'])
        denom=4*len(ys)
        for mod,sse in sumloss.items():eq(f'{k}.meanloss.{mod}',sse/denom,results[k]['mean_frequency_MSE'][mod])
        gain=(sumloss['GLOBAL']-sumloss['SIDE'])/sumloss['GLOBAL']
        eq(f'{k}.relative_gain',gain,results[k]['relative_error_reduction_from_side'])
    eq('censusratio',max(counts[y] for y in years)/min(counts[y] for y in years),results['annual_abundance_max_over_min_main'])
    return audit

def main():
    ap=argparse.ArgumentParser()
    base=Path(__file__).resolve().parents[1]
    ap.add_argument('--source',type=Path,default=base/'source/linanthus_counts.csv')
    ap.add_argument('--result',type=Path,default=base/'results/result.json')
    ap.add_argument('--out',type=Path,default=base/'results/independent_verification.json')
    args=ap.parse_args()
    z=check(args.source,json.loads(args.result.read_text()))
    args.out.write_text(json.dumps(z,indent=2,sort_keys=True)+'\n')
    print(json.dumps(z))
if __name__=='__main__':main()
