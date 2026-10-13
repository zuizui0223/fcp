#!/usr/bin/env python3
"""Independent raw-source and Newton-Hessian regression verifier; never import predictive.py."""
import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import expit

ORIGINAL_ARCHIVE_MD5='008ec196aa5acb769dbfd36e835e32ba'
ORIGINAL_ARCHIVE_SHA256='4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
ORIGINAL_FREQUENCY_SHA256='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
FEATURES={'BASE':[], 'DIST':['distance'], 'CONGEN':['other_share'], 'BOTH':['distance','other_share']}


def read_original(path):
    raw=Path(path).read_bytes()
    assert len(raw)==127502
    assert hashlib.md5(raw).hexdigest()==ORIGINAL_ARCHIVE_MD5
    assert hashlib.sha256(raw).hexdigest()==ORIGINAL_ARCHIVE_SHA256
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        assert z.testzip() is None
        b=z.read('data/frequency.txt')
    assert hashlib.sha256(b).hexdigest()==ORIGINAL_FREQUENCY_SHA256
    rdr=csv.DictReader(io.StringIO(b.decode('utf-8-sig').replace('\r','\n')),delimiter='\t')
    original=list(rdr)
    assert len(original)==512
    rows=[]
    for x in original:
        data={k:x[k] for k in ('year','population','quad')}
        for k in ('line','distance','v','w','l','s'):data[k]=int(x[k])
        data['n']=data['v']+data['w']
        if data['n']>0:data['other_share']=(data['l']+data['s'])/(data['n']+data['l']+data['s'])
        rows.append(data)
    assert len({(r['year'],r['population'],r['line'],r['quad']) for r in rows})==len(rows)
    assert tuple(sum(r[k] for r in rows) for k in ('v','w','l','s'))==(5081,8723,4463,29)
    return rows


def independent_fit(rows,heldout,features,lam=10.0):
    n=np.asarray([z['n'] for z in rows],float)
    v=np.asarray([z['v'] for z in rows],float)
    sitebase=(v.sum()+.5)/(n.sum()+1)
    if not features:return np.full(len(heldout),min(.99,max(.01,sitebase)))
    x=np.array([[r[f] for f in features] for r in rows],dtype=float)
    test_x=np.array([[r[f] for f in features] for r in heldout],dtype=float)
    mu=np.sum(x*n[:,None],axis=0)/np.sum(n)
    var=np.sum(n[:,None]*(x-mu)**2,axis=0)/n.sum()
    std=np.sqrt(var)
    ident=std>1e-12
    x=(x-mu)/np.where(ident,std,1)
    test_x=(test_x-mu)/np.where(ident,std,1)
    x[:,~ident]=0
    test_x[:,~ident]=0
    X=np.column_stack([np.ones(len(n)),x]);TX=np.column_stack([np.ones(len(heldout)),test_x])
    theta=np.zeros(X.shape[1]);theta[0]=np.log((v.sum()+.5)/(n.sum()-v.sum()+.5))
    def objective(b):
        z=X@b
        return float(np.dot(n,np.logaddexp(0,z))-np.dot(v,z)
            +np.logaddexp(0,b[0])-.5*b[0]+.5*lam*np.dot(b[1:],b[1:]))
    for itr in range(100):
        z=X@theta;q=expit(z)
        grad=X.T@(n*q-v)
        grad[0]+=expit(theta[0])-.5
        grad[1:]+=lam*theta[1:]
        if np.max(np.abs(grad))<1e-8:break
        w=n*q*(1-q)
        h=X.T@(X*w[:,None]);h[0,0]+=expit(theta[0])*(1-expit(theta[0]))
        h[1:,1:]+=lam*np.eye(len(theta)-1)
        step=np.linalg.solve(h,grad)
        alpha=1.0
        old_obj=objective(theta)
        while alpha>=1/2**25 and objective(theta-alpha*step)>old_obj+1e-11:
            alpha*=0.5
        theta=theta-alpha*step
    last_grad=X.T@(n*expit(X@theta)-v)
    last_grad[0]+=expit(theta[0])-.5
    last_grad[1:]+=lam*theta[1:]
    if np.max(np.abs(last_grad))>1e-5:
        raise ArithmeticError(f'INDEPENDENT_NEWTON_NOT_CONVERGED, max_gradient={np.max(abs(last_grad))}')
    return np.clip(expit(TX@theta),.01,.99)


def verify(root):
    root=Path(root)
    rows=read_original(root/'source'/'Sisyrinchium_data.zip')
    official=json.loads((root/'results'/'result.json').read_text())
    assert official['schema']=='fcp_sisyrinchium_heldout_field_line_distance_congener_prediction_v1'
    checks=0;worst=0.
    by_site=defaultdict(list)
    for x in rows:
        if x['n']>0 and x['distance']>=0:
            by_site[(x['year'],x['population'])].append(x)
    original_pred=official['predictions']
    pred_bykey={(p['year'],p['population'],p['line'],p['quad']):p for p in original_pred}
    assert len(pred_bykey)==len(original_pred)
    siteind=defaultdict(list)
    for site in official['eligible_siteyears']:
        year,pop=site['year'],site['population']
        data=by_site[(year,pop)]
        # Independently verify frozen inclusion shape (without using model output).
        plus=[x for x in data if x['l']+x['s']>0];zero=[x for x in data if x['l']+x['s']==0]
        assert sum(x['n'] for x in data)>=50
        assert sum(x['v'] for x in data)>=5
        assert sum(x['w'] for x in data)>=5
        assert len(plus)>=3 and len(zero)>=3
        assert len(set(x['line'] for x in plus))>=2
        assert len(set(x['line'] for x in zero))>=2
        assert len(site['heldout_lines'])>=3
        checks+=9
        for line in site['heldout_lines']:
            train=[x for x in data if x['line']!=line]
            test=[x for x in data if x['line']==line]
            assert len(test)>=2 and sum(x['n'] for x in test)>=10
            assert sum(x['n'] for x in train)>=50
            assert sum(x['v'] for x in train)>=5
            assert sum(x['w'] for x in train)>=5
            checks+=4
            pre={m:independent_fit(train,test,f) for m,f in FEATURES.items()}
            for i,actual in enumerate(test):
                key=(year,pop,line,actual['quad'])
                p=pred_bykey[key]
                for m in FEATURES:
                    delta=abs(float(pre[m][i])-p['p_'+m]);worst=max(worst,delta)
                    if delta>2e-7:raise ArithmeticError(f'INDEPENDENT_PREDICTION_DISAGREES {key} {m} {delta}')
                    checks+=1
                assert actual['v']==p['v'] and actual['w']==p['w'] and actual['l']==p['l'] and actual['s']==p['s']
                siteind[(year,pop)].append((actual,pre,i))
                checks+=4
    assert len(pred_bykey)==len(official['predictions'])==sum(len(x) for x in siteind.values())
    checks+=1
    for result in official['site_scores']:
        site=(result['year'],result['population']); rr=siteind[site]
        nsum=sum(x['n'] for x,_,_ in rr)
        assert nsum==result['n_scored_plants']
        checks+=1
        for model in FEATURES:
            score=0.
            for x,p,i in rr:
                q=float(p[model][i]);score+=-x['v']*math.log(q)-x['w']*math.log1p(-q)
            score/=nsum
            delta=abs(score-result['logloss_'+model]);worst=max(worst,delta)
            if delta>1e-7:raise ArithmeticError(f'INDEPENDENT_SCORE_DISAGREES {site} {model} {delta}')
            checks+=1
    main=official['main']['mean']
    sampled=official['selected_latest']
    for model in FEATURES:
        avg=sum(x['logloss_'+model] for x in sampled)/len(sampled)
        assert abs(avg-main['logloss_'+model])<1e-13
        checks+=1
    # Ratio of best source scores is not interpreted as biology significance.
    report={'status':'PASS','checks':checks,'independent_method':'author-original CR-delimited UTF8 CSV + independent weighted Newton-Hessian logistic + Jeffreys intercept prior + ridge','n_folds':len(official['fold_ledger']),'n_independent_heldout_quadrats':len(pred_bykey),'n_latest_named_sites':len(sampled),'max_abs_difference':worst,'source_sha256':ORIGINAL_ARCHIVE_SHA256,'not_new_field_data':True,'not_causal':True,'not_new_geographic_site_transfer':True}
    (root/'results'/'independent_verification.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',default='.');args=parser.parse_args()
    verify(args.root)
