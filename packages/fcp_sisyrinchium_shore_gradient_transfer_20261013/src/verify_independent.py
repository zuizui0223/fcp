#!/usr/bin/env python3
"""Independent direct source reader and positive-Hessian Newton verification.

Does NOT import the production transport.py or posthoc_standardized.py.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq

REQ_MD5='008ec196aa5acb769dbfd36e835e32ba'
REQ_SHA='4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
REQ_TXT='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
LAMBDA=10

def original_independent(path):
    data=Path(path).read_bytes()
    assert len(data)==127502 and hashlib.md5(data).hexdigest()==REQ_MD5 and hashlib.sha256(data).hexdigest()==REQ_SHA
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert z.testzip() is None
        member=z.read('data/frequency.txt')
    assert hashlib.sha256(member).hexdigest()==REQ_TXT
    lines=member.decode('utf-8-sig').replace('\r','\n').splitlines()
    header=lines[0].split('\t');assert header==['year','population','line','quad','distance','humidity','v','w','l','s']
    out={}
    for line in lines[1:]:
        if not line:continue
        c=line.split('\t');assert len(c)==len(header)
        d=dict(zip(header,c))
        key=(d['year'],d['population'],int(d['line']),d['quad'])
        assert key not in out
        counts={k:int(d[k]) for k in ('v','w','l','s')}
        assert all(v>=0 for v in counts.values())
        out[key]={'year':key[0],'population':key[1],'siteyear':key[0]+'|'+key[1], 'line':key[2], 'quad':key[3],
                  'distance':int(d['distance']),'n':counts['v']+counts['w'],**counts}
    assert len(out)==512
    assert [sum(z[k] for z in out.values()) for k in ('v','w','l','s')]==[5081,8723,4463,29]
    return out

def independent_site_ridge(rows):
    n=np.array([z['n'] for z in rows],float)
    v=np.array([z['v'] for z in rows],float)
    dist=np.array([z['distance'] for z in rows],float)
    assert len(rows)>=3 and n.sum()>25 and v.sum()>=5 and (n-v).sum()>=5
    center=np.dot(dist,n)/n.sum()
    sd=math.sqrt(np.dot((dist-center)**2,n)/n.sum())
    assert sd>1e-11
    x=(dist-center)/sd
    X=np.column_stack([np.ones(len(x)),x])
    intercept=math.log((v.sum()+.5)/(n.sum()-v.sum()+.5))
    theta=np.array([intercept,0.],float)
    def fun(theta):
        xx=X@theta
        return float(np.dot(n,np.logaddexp(0,xx))-np.dot(v,xx)
             +np.logaddexp(0,theta[0])-.5*theta[0]+5*theta[1]**2)
    for _ in range(60):
        z=X@theta
        p=expit(z)
        delta=n*p-v
        g=X.T@delta + np.array([expit(theta[0])-.5, LAMBDA*theta[1]])
        if max(abs(g))<1e-8:break
        weights=n*p*(1-p)
        h=X.T@(X*weights[:,None])
        h[0,0]+=expit(theta[0])*(1-expit(theta[0]))
        h[1,1]+=LAMBDA
        step=np.linalg.solve(h,g)
        f0=fun(theta)
        a=1.0
        while a>1e-8 and fun(theta-a*step)>f0+1e-10:a*=.5
        theta-=a*step
    p=expit(X@theta);res=X.T@(n*p-v)+np.array([expit(theta[0])-.5,LAMBDA*theta[1]])
    assert max(abs(res))<1e-6
    return {'per_distance':float(theta[1]/sd),'standardized':float(theta[1]),'sd':float(sd), 'n_M':int(n.sum())}

def target_base_intercept(rows,coefficient):
    d=np.array([r['distance'] for r in rows],float);n=np.array([r['n'] for r in rows],float);v=np.array([r['v'] for r in rows],float)
    mu=np.dot(n,d)/n.sum()
    contrast=coefficient*(d-mu)
    def f(a):return np.dot(n,expit(a+contrast))-v.sum()+expit(a)-.5
    intercept=brentq(f,-70,70,xtol=1e-13)
    assert abs(f(intercept))<1e-7
    return float(intercept),float(mu)

def compare(a,b,checks,tol=3e-9):
    err=abs(float(a)-float(b));assert err<=tol,(a,b,err)
    checks.append(err)

def verify(original_zip,primary_result,posthoc_result,outdir):
    orig=original_independent(original_zip)
    p=json.loads(Path(primary_result).read_text());post=json.loads(Path(posthoc_result).read_text())
    assert p['schema']=='fcp_sisyrinchium_shoreline_gradient_cross_population_transport_v1'
    assert post['status']=='POSTHOC_DIAGNOSTIC_NOT_PRIMARY'
    sites=sorted(p['source_fit']);assert len(sites)==6
    rows_by_site={s:[r for r in orig.values() if r['siteyear']==s and r['n']>0 and r['distance']>=0] for s in sites}
    source_fit={s:independent_site_ridge(rows_by_site[s]) for s in sites}
    checks=[]
    for s in sites:
        for k,refk in [('per_distance','beta_per_distance'),('standardized','beta_standardized'),('sd','weighted_distance_sd')]:
            compare(source_fit[s][k],p['source_fit'][s][refk],checks)
        assert source_fit[s]['n_M']==p['source_fit'][s]['n_M']
        checks.append(0.)
    bykey={(r['year'],r['population'],r['line'],r['quad']):r for r in p['predictions']}
    zh={(r['siteyear'],int(r['line']),r['quad']):r for r in csv.DictReader((Path(posthoc_result).parent/'posthoc_distance_scale_predictions.csv').open(newline=''))}
    assert len(bykey)==len(zh)==168
    full_numeric={s:defaultdict(float) for s in sites}
    fold_count=0
    for fold in p['folds']:
        s=fold['siteyear'];line=fold['heldout_line'];others=[x for x in sites if x!=s]
        assert s not in others and others==fold['source_sites']
        median=float(np.median([source_fit[x]['per_distance'] for x in others]))
        mean=float(np.mean([source_fit[x]['per_distance'] for x in others]))
        zmed=float(np.median([source_fit[x]['standardized'] for x in others]))
        zmean=float(np.mean([source_fit[x]['standardized'] for x in others]))
        for lhs,rhs in [(median,fold['borrow_median_beta_per_distance']),(mean,fold['borrow_mean_beta_per_distance'])]:compare(lhs,rhs,checks)
        train=[x for x in rows_by_site[s] if x['line']!=line]
        test=[x for x in rows_by_site[s] if x['line']==line and (x['year'],x['population'],x['line'],x['quad']) in bykey]
        assert sum(x['n'] for x in train)==fold['train_n_M']
        assert sum(x['n'] for x in test)==fold['test_n_M']
        local_slope=independent_site_ridge(train)['per_distance']
        compare(local_slope,fold['target_training_local_beta'],checks)
        a,mu=target_base_intercept(train,median)
        a2,mu2=target_base_intercept(train,mean)
        for observed,original in [(a,fold['target_intercept_calibration']),(a2,fold['target_intercept_calibration_mean']),(mu,fold['target_train_distance_center'])]:compare(observed,original,checks)
        assert abs(mu-mu2)<1e-10
        n=np.array([x['n'] for x in train],float);d=np.array([x['distance'] for x in train],float)
        target_sd=math.sqrt(np.dot(n,(d-mu)**2)/n.sum())
        baz,cz=target_base_intercept(train,zmed/target_sd)
        baz2,cz2=target_base_intercept(train,zmean/target_sd)
        for item in test:
            key=(item['year'],item['population'],item['line'],item['quad'])
            res=bykey[key]
            qmed=float(np.clip(expit(a+median*(item['distance']-mu)),.01,.99))
            qmean=float(np.clip(expit(a2+mean*(item['distance']-mu)),.01,.99))
            qz=float(np.clip(expit(baz+zmed/target_sd*(item['distance']-cz)),.01,.99))
            qz2=float(np.clip(expit(baz2+zmean/target_sd*(item['distance']-cz2)),.01,.99))
            q2=zh[(s,line,item['quad'])]
            for calc,saved in [(qmed,res['p_BORROW_MEDIAN']),(qmean,res['p_BORROW_MEAN']),
                               (qz,float(q2['p_BORROW_Z_MEDIAN'])),(qz2,float(q2['p_BORROW_Z_MEAN']))]:compare(calc,saved,checks)
            full_numeric[s]['n']+=item['n']
            for name,q in [('BASE',res['p_BASE']),('LOCAL_DIST',res['p_LOCAL_DIST']),('BORROW_MEDIAN',qmed),('BORROW_MEAN',qmean),('BORROW_Z_MEDIAN',qz),('BORROW_Z_MEAN',qz2)]:
                qq=max(.01,min(.99,q));full_numeric[s][name]+=(-item['v']*math.log(qq)-item['w']*math.log1p(-qq))
            assert item['v']==res['v'] and item['w']==res['w'] and item['n']==res['n']
            checks.append(0.)
        fold_count+=1
    assert fold_count==32 and sum(x['n'] for x in full_numeric.values())==6692
    for s in sites:
        sr=next(x for x in p['sites'] if x['siteyear']==s)
        postr=next(x for x in post['sites'] if x['siteyear']==s)
        for name in ['BASE','LOCAL_DIST','BORROW_MEDIAN','BORROW_MEAN']:
            compare(full_numeric[s][name]/full_numeric[s]['n'],sr['logloss_'+name],checks)
        for name in ['BASE','LOCAL_DIST','BORROW_Z_MEDIAN','BORROW_Z_MEAN']:
            compare(full_numeric[s][name]/full_numeric[s]['n'],postr['score'][name],checks)
    # Scientific score identity, independent from production aggregation code.
    for name in ['BASE','LOCAL_DIST','BORROW_MEDIAN','BORROW_MEAN']:
        expected=sum(full_numeric[s][name]/full_numeric[s]['n'] for s in sites)/6
        compare(expected,p['main']['mean']['logloss_'+name],checks)
    for name in ['BASE','LOCAL_DIST','BORROW_Z_MEDIAN','BORROW_Z_MEAN']:
        expected=sum(full_numeric[s][name]/full_numeric[s]['n'] for s in sites)/6
        compare(expected,post['main']['mean_logloss'][name],checks)
    assert sum(source_fit[s]['per_distance']<0 for s in sites)==5
    assert source_fit['Y13|banpaku']['per_distance']>0
    assert sum(1 for s in sites if next(x for x in p['sites'] if x['siteyear']==s)['gain_local_vs_base']>0)==6
    receipt={'schema':'fcp_sisyrinchium_gradient_transfer_independent_newton_verification_v1',
             'status':'PASS','method':'Standalone original CR-delimited field TXT read + site-ridge Newton with exact Hessian + independent target intercept root + all holdout predictions',
             'n_independent_arithmetic_checks':len(checks),'max_absolute_numeric_difference':max(checks),
             'source_verified_by_hash':True,'no_production_module_import':True,
             'n_source_sites':6,'n_outer_sites':6,'n_inner_heldout_lines':fold_count,'n_heldout_quadrats':len(bykey),
             'n_heldout_M_plants':int(sum(x['n'] for x in full_numeric.values())),
             'source_sign_count_negative':5,'source_sign_count_positive':1,
             'no_geographically_new_field_data':True,'not_selection_or_causality':True}
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/'independent_verification.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,ensure_ascii=False))
    return receipt

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',required=True);ap.add_argument('--result',required=True);ap.add_argument('--posthoc',required=True);ap.add_argument('--outdir',required=True)
    z=ap.parse_args()
    verify(z.source,z.result,z.posthoc,z.outdir)
