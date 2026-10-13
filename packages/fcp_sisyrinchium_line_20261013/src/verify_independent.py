#!/usr/bin/env python3
"""Independent exact-rational verification from original source TXT; NO import of within_line.

Reconstructs 16 siteyears, all eligible line-conditioned original numerators,
analytic null expectations, plus synthetic finite-label exact enumeration.
"""
from __future__ import annotations
import argparse,csv,json,hashlib,math,itertools
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

EXPECTED='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
SCOREKEY={'exposure':'T_congener_purple_minus_white','distance':'T_congener_after_linear_distance','clumping':'residual_delta'}
R=9999


def original_rows(path):
    b=Path(path).read_bytes()
    assert hashlib.sha256(b).hexdigest()==EXPECTED
    allrows=list(csv.DictReader(b.decode('utf-8-sig').splitlines(),delimiter='\t'))
    assert len(allrows)==512
    for r in allrows:
        for k in ['v','w','l','s','line','distance']:
            r[k]=int(r[k])
        r['n']=r['v']+r['w'];r['a']=F(r['l']+r['s'],r['n']+r['l']+r['s']) if r['n']+r['l']+r['s']>0 else None
    return allrows


def groups_from_rows(rows):
    siteyears=defaultdict(list)
    for r in rows:siteyears[(r['year'],r['population'])].append(r)
    result=[]
    for (y,s),r in sorted(siteyears.items()):
        lines=defaultdict(list)
        for a in r:
            if a['n']>0:lines[a['line']].append(a)
        groups=[]
        N=sum(x['n'] for x in r);V=sum(x['v'] for x in r);W=N-V
        for lid,plots in sorted(lines.items()):
            Nl=sum(x['n'] for x in plots);Vl=sum(x['v'] for x in plots);Wl=Nl-Vl
            pp=F(Vl,Nl)
            aa=[x['a'] for x in plots]
            spread=(max(aa)-min(aa)) if aa else F(0)
            info=(len(plots)>=2 and Nl>=10 and Vl>=1 and Wl>=1 and spread>0)
            Ds=F(Vl*Wl,Nl)
            da=None;dok=False
            if info and len(plots)>=3 and len(set(x['distance'] for x in plots))>=2:
                mean_a=sum(F(x['n'])*x['a'] for x in plots)/Nl
                mean_d=sum(F(x['n'])*x['distance'] for x in plots)/Nl
                cov=sum(F(x['n'])*(x['distance']-mean_d)*(x['a']-mean_a) for x in plots)
                va=sum(F(x['n'])*(x['distance']-mean_d)**2 for x in plots)
                beta=cov/va if va else F(0)
                da=[x['a']-mean_a-beta*(x['distance']-mean_d) for x in plots]
                dok=(max(da)-min(da)>F(1,10**12))
            groups.append(dict(line=lid,plots=plots,N=Nl,V=Vl,W=Wl,p=pp,D=Ds,info=info,distance_ok=dok,a_residual=da))
        record=dict(siteyear=y+'|'+s,site=s,year=y,N=N,V=V,W=W,n_quadrats=len(r),n_occupied=sum(x['n']>0 for x in r),
                    lines=groups,base_eligible=N>=20 and V>=5 and W>=5 and sum(x['n']>0 for x in r)>=2)
        if record['base_eligible']:
            p=F(V,N);den=N*p*(1-p)
            b_num=sum((F(x['v'])-F(x['n'])*p)**2/x['n'] for x in r if x['n']>0)
            bb_num=sum(g['N']*(g['p']-p)**2 for g in groups)
            bexp_num=bb_num+sum(g['D']*F(len(g['plots'])-1,g['N']-1) for g in groups if len(g['plots'])>1 and g['N']>1)
            record.update(B=b_num/den,B_between=bb_num/den,E=bexp_num/den,delta=(b_num-bexp_num)/den)
            expr=[g for g in groups if g['info']]
            Nd=sum(g['N'] for g in expr); Vd=sum(g['V'] for g in expr)
            record['primary_eligible']=len(expr)>=2 and Nd>=30 and Vd>=5 and Nd-Vd>=5
            record['primary_nlines']=len(expr);record['primary_N']=Nd
            if record['primary_eligible']:
                numer=sum(sum((F(x['v'])-F(x['n'])*g['p'])*x['a'] for x in g['plots']) for g in expr)
                D=sum(g['D'] for g in expr)
                record['primary_T']=numer/D
            dr=[g for g in groups if g['distance_ok']]
            Nr=sum(g['N'] for g in dr);Vr=sum(g['V'] for g in dr)
            record['distance_eligible']=len(dr)>=2 and Nr>=30 and Vr>=5 and Nr-Vr>=5
            record['distance_nlines']=len(dr);record['distance_N']=Nr
            if record['distance_eligible']:
                numer=sum(sum((F(x['v'])-F(x['n'])*g['p'])*a for x,a in zip(g['plots'],g['a_residual'])) for g in dr)
                D=sum(g['D'] for g in dr);record['distance_T']=numer/D
            record['clumping_eligible']=sum(len(g['plots'])>=2 and g['V']>0 and g['W']>0 for g in groups)>=2
        else:
            record.update(primary_eligible=False,distance_eligible=False,clumping_eligible=False)
        result.append(record)
    return result


def microcase_exact():
    # N=4, V=2, two equal plots -> E[B]=(m-1)/(N-1)=1/3.
    N=4;V=2;n=(2,2);p=F(V,N)
    space=list(itertools.combinations(range(N),V));bs=[];ts=[]
    for inds in space:
        v1=sum(x<2 for x in inds);v2=V-v1
        B=sum((F(v)-k*p)**2/k for v,k in zip([v1,v2],n))/(N*p*(1-p))
        exposure=[F(1,2),F(0)];score=sum((F(v)-k*p)*a for v,k,a in zip([v1,v2],n,exposure))/(N*p*(1-p))
        bs.append(B);ts.append(score)
    assert sum(bs)/len(bs)==F(1,3)
    assert sum(ts)==0
    assert bs.count(F(1))==2 and bs.count(F(0))==4
    return 3


def verify(source, result_path, out):
    rows=original_rows(source);snap=json.loads(Path(result_path).read_text())
    assert snap['schema']=='fcp_sisyrinchium_line_conditional_cogeneric_exposure_20261013_v1'
    assert snap['R']==R
    checks=microcase_exact()
    sites=groups_from_rows(rows)
    assert len(sites)==16
    d={(r['siteyear']):r for r in snap['full']['siteyear_ledger']}
    assert len(d)==16
    exact_source_totals={'n_field_quadrat_records':512,'M_v':5081,'M_w':8723,'other_L':4463,'other_S':29}
    for k,v in exact_source_totals.items():assert snap['full']['source'][k]==v;checks+=1
    eligible=defaultdict(list)
    differences=[]
    for rec in sites:
        q=d[rec['siteyear']]
        for col in ['N','V','W','m','source_quadrats']:
            orig={'m':rec['n_occupied'],'source_quadrats':rec['n_quadrats']}.get(col,rec.get(col))
            assert q[col]==orig,(rec['siteyear'],col,q.get(col),orig)
            checks+=1
        assert (q['status_site']=='ELIGIBLE')==rec['base_eligible'];checks+=1
        if not rec['base_eligible']:
            continue
        vals=[('B_original','B'),('B_between_recorded_lines','B_between'),('B_conditioned_E','E'),('residual_delta','delta')]
        for src,dest in vals:
            got=float(rec[dest]);actual=q[src]
            diff=abs(got-actual)
            assert diff<1e-12,(rec['siteyear'],src,got,actual)
            differences.append(diff);checks+=1
        for k,k2,k3 in [('primary','primary_eligible','primary_T'),('distance','distance_eligible','distance_T')]:
            enabled=rec[k2];f=k+'_status'
            assert (q[f]=='ELIGIBLE')==enabled,(rec['siteyear'],f,enabled,q[f]);checks+=1
            if enabled:
                name='T_congener_purple_minus_white' if k=='primary' else 'T_congener_after_linear_distance'
                diff=abs(float(rec[k3])-q[name]); assert diff<1e-12;differences.append(diff);checks+=1
                eligible['exposure' if k=='primary' else 'distance'].append(float(rec[k3]))
        assert (q['clumping_status']=='ELIGIBLE')==rec['clumping_eligible'];checks+=1
        if rec['clumping_eligible']:eligible['clumping'].append(float(rec['delta']))
    for name,a in eligible.items():
        summary=snap['full']['results'][name]
        assert summary['status']=='ELIGIBLE'
        assert summary['n']==len(a)
        assert abs(summary['estimate']-sum(a)/len(a))<1e-12
        checks+=3
        differences.append(abs(summary['estimate']-sum(a)/len(a)))
    assert len(eligible['exposure'])==8 and len(eligible['clumping'])==11 and len(eligible['distance'])==8
    for name in ['exposure','clumping','distance']:
        assert len(snap['full']['results'][name]['siteyear_estimates'])==len(eligible[name]);checks+=1
    assert snap['full']['results']['distance']['posthoc_leave_one_siteyear_out']
    assert snap['full']['results']['exposure']['p_one_sided']<=.001
    assert snap['full']['results']['clumping']['p_one_sided']<=.001
    assert snap['source_paper_already_reported_cross_site_congener_association'] is True
    out.mkdir(exist_ok=True,parents=True)
    receipt={'status':'PASS','original_independent_source_sha256':EXPECTED,
             'independent_method':'csv.DictReader and fractions.Fraction including exact observed site×line and analytical null, microcase exact enumeration; does not import main script',
             'checks':checks,'max_absolute_difference':max(differences) if differences else 0,
             'eligible_siteyears':{k:len(v) for k,v in eligible.items()},
             'simulation_crosscheck':'source-conditioned permutation tail must be rerun by audited primary code; no bogus independent field data',
             'congener_species_are_not_white_M_morph':True,'not_causal_test':True}
    (out/'independent_verification.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt))
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default='source/frequency.txt');p.add_argument('--result',default='results/result.json');p.add_argument('--outdir',default='results');a=p.parse_args();verify(Path(a.source),Path(a.result),Path(a.outdir))
