#!/usr/bin/env python3
"""Independent csv/Fraction reconstruction; never imports production analysis."""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import numpy as np

SITES=['G_1','H_1','H_3','GPS-143','GPS-160']
KEY=['COUNTRY','Year','GPS-code','Plant_ID']
MD5='a1200e0197b282ee126700ae234eb95e'
SHA='723c7a03b6ffa010e6093996d8203b6770ab1ee1bc6ddc21ec9c3d5b47d3b8c9'


def close(got, expected, label, tol=1e-10):
    if not math.isclose(float(got),float(expected),abs_tol=tol,rel_tol=tol):
        raise AssertionError(f'{label}: {got} != {expected}')


def verify(source,outdir):
    b=source.read_bytes()
    assert len(b)==345652 and hashlib.md5(b).hexdigest()==MD5 and hashlib.sha256(b).hexdigest()==SHA
    raw=list(csv.DictReader(b.decode('cp1252').splitlines()))
    assert len(raw)==1323
    grouped=defaultdict(list)
    for x in raw:
        assert all(x[k].strip() for k in KEY)
        grouped[tuple(x[k] for k in KEY)].append(x)
    led=list(csv.DictReader((outdir/'plant_patch_ledger.csv').open()))
    assert len(led)==len(grouped)
    result=json.loads((outdir/'result.json').read_text())
    checks=0
    fractions={};lengths={};days={};morphs={}
    for row in led:
        key=tuple(row[k] for k in KEY)
        g=grouped[key]
        y=sum(x['Robbed_flower'].strip().lower()=='yes' for x in g)
        n=sum(x['Robbed_flower'].strip().lower()=='no' for x in g)
        colors={x['Flower_morph'].strip().lower() for x in g if x['Flower_morph'].strip()}
        c='CONFLICT' if len(colors)>1 else next(iter(colors)) if len(colors)==1 else 'MISSING'
        dates={x['Date'].strip() for x in g if x['Date'].strip()}
        assert len(g)==int(row['n_flower_rows']);checks+=1
        assert y==int(row['n_yes']);checks+=1
        assert n==int(row['n_no']);checks+=1
        assert c==row['morph'];checks+=1
        assert len(dates)==int(row['n_dates']);checks+=1
        if y+n:
            f=Fraction(y,y+n)
            close(f,row['robbed_fraction'],f'{key} robbed fraction');checks+=1
            fractions[key]=f
        else:
            assert row['robbed_fraction']=='' or row['robbed_fraction'].lower()=='nan';checks+=1
        vv=[Decimal(x['Length_flower']) for x in g if x['Length_flower'].strip() not in ('','NA')]
        if vv:
            mm=sum(vv)/len(vv)
            close(mm,row['length_plant_mean'],f'{key} length mean',tol=1e-9);checks+=1
            lengths[key]=float(mm)
        days[key]=dates;morphs[key]=c
    def strata(endpoint='robbed'):
        arrays=[]
        vals=fractions if endpoint=='robbed' else lengths
        for site in SITES:
            v=[(key,vals[key]) for key in sorted(vals) if key[1]=='2024' and key[2]==site and morphs[key]=='violet']
            p=[(key,vals[key]) for key in sorted(vals) if key[1]=='2024' and key[2]==site and morphs[key]=='pink']
            if len(v)<5 or len(p)<5:raise AssertionError(f'insufficient {site}')
            arrays.append((site,v,p))
        return arrays
    raw_s=strata('robbed')
    for site,v,p in raw_s:
        saved=next(x for x in result['main_site_contrasts'] if x['site']==site)
        assert len(v)==saved['n_violet_groups'];checks+=1
        assert len(p)==saved['n_pink_groups'];checks+=1
        ev=sum(x for _,x in v)/len(v);ep=sum(x for _,x in p)/len(p)
        close(ev,saved['violet_mean'],site+' violet');checks+=1
        close(ep,saved['pink_mean'],site+' pink');checks+=1
        close(ev-ep,saved['difference_violet_minus_pink'],site+' contrast');checks+=1
    point=sum((sum(x for _,x in v)/len(v)-sum(x for _,x in p)/len(p) for _,v,p in raw_s),Fraction(0))/len(raw_s)
    close(point,result['main_paired_colour_reference']['site_equal_estimate'],'primary exact Fraction T');checks+=1
    # Independent permutation using source-ledger enumerated site-grouped arrays.
    rng=np.random.default_rng(2026101301)
    vals=[]
    for site,v,p in raw_s:
        vals.append((np.array([float(x) for _,x in v]),np.array([float(x) for _,x in p])))
    exceeds=0
    for _ in range(9999):
        tmp=[]
        for v,p in vals:
            ix=rng.permutation(len(v)+len(p))
            x=np.r_[v,p]
            tmp.append(x[ix[:len(v)]].mean()-x[ix[len(v):]].mean())
        if abs(np.mean(tmp))>=abs(float(point))-1e-12:exceeds+=1
    gotp=(exceeds+1)/10000
    close(gotp,result['main_paired_colour_reference']['p_two_sided'],'primary source-only RNG permutation',tol=1e-12);checks+=1
    # Rebuild date sensitivity using independently obtained plant dates.
    valid_by_day=defaultdict(lambda:{'violet':[],'pink':[]})
    for key,fr in fractions.items():
        if key[1]!='2024' or key[2] not in SITES or morphs[key] not in ('violet','pink') or len(days[key])!=1:continue
        d=next(iter(days[key]))
        valid_by_day[key[2],d][morphs[key]].append(float(fr))
    parts=[]
    for (site,date),mp in sorted(valid_by_day.items()):
        if len(mp['violet'])>=2 and len(mp['pink'])>=2:
            parts.append((site,date,np.array(mp['violet']),np.array(mp['pink'])))
    ndate=sum(len(v)+len(p) for _,_,v,p in parts)
    assert ndate==result['date_blocked_sensitivity']['n_supported_groups'];checks+=1
    assert len(set(x[0] for x in parts))==result['date_blocked_sensitivity']['n_supported_sites'];checks+=1
    # Validate independent date sensitivity statistic, allowing floating rounding.
    dx=defaultdict(list)
    for site,date,v,p in parts:dx[site].append(v.mean()-p.mean())
    datet=float(np.mean([np.mean(dx[s]) for s in SITES if s in dx]))
    close(datet,result['date_blocked_sensitivity']['site_equal_estimate'],'date-blocked T');checks+=1
    len_s=strata('length')
    length_t=float(np.mean([np.mean([x for _,x in v])-np.mean([x for _,x in p]) for site,v,p in len_s]))
    close(length_t,result['secondary_flower_length_mm']['within_site_reference']['site_equal_estimate'],'length T');checks+=1
    qc=result['quality']
    assert qc['n_source_flower_rows']==1323 and qc['n_2024_primary_site_flower_rows']==977;checks+=1
    assert qc['n_grouped_2024_five_sites']==226;checks+=1
    assert qc['source']['SHA256']==SHA;checks+=1
    resultobj={'schema':'fcp_digitalis_robbed_source_independent_verification_v1','status':'PASS','checks':checks,
               'source_md5':MD5,'source_sha256':SHA,'source_n_flower_rows':1323,
               'n_grouped_units':len(grouped),'main_fraction_exact_numerator':point.numerator,
               'main_fraction_exact_denominator':point.denominator,
               'main_fraction_float':float(point),
               'permutation_two_sided_p':gotp,'verified_date_overlap_units':ndate,
               'morphology_main_mm_diff':length_t,
               'verification_boundary':'Separate csv.DictReader + rational arithmetic and independently reconstructed permutation; same original study and source, not independent biological replication'}
    (outdir/'independent_verification.json').write_text(json.dumps(resultobj,indent=2)+'\n')
    return resultobj

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    a=ap.parse_args();print(json.dumps(verify(a.source,a.outdir),indent=2))
