#!/usr/bin/env python3
"""Independent original-file reader and Newton solver for geographic holdout audit.

DOES NOT import transfer.py nor the predecessor's predictive.py. It independently
constructs all site-heldout folds from publisher bytes and frozen PR185 test keys.
"""
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np

ZIP_MD5='008ec196aa5acb769dbfd36e835e32ba'
ZIP_SHA='4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
FREQ_SHA='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
PARENT_SHA='7745ab3018eb00a1939b9433effb201e1c9cc37034f2a4c2ae8a588e3c2b23ce'
FEATURES={'BASE':[],'DIST':['distance'],'CONGEN':['other_share'],'BOTH':['distance','other_share']}

def parse_actual_original(zippath):
    original=Path(zippath).read_bytes()
    assert len(original)==127502 and hashlib.md5(original).hexdigest()==ZIP_MD5
    assert hashlib.sha256(original).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(original)) as arc:
        assert arc.testzip() is None
        txt=arc.read('data/frequency.txt')
    assert hashlib.sha256(txt).hexdigest()==FREQ_SHA
    # Original frequency.txt uses Macintosh CR line endings, not rows of flowers.
    lines=list(csv.DictReader(io.StringIO(txt.decode('utf-8').replace('\r','\n')),delimiter='\t'))
    assert len(lines)==512
    rows={}
    for q in lines:
        key=(q['year'],q['population'],int(q['line']),q['quad'])
        assert key not in rows
        v,w,l,s=[int(q[k]) for k in ('v','w','l','s')]
        assert min(v,w,l,s)>=0
        n=v+w
        a=(l+s)/(n+l+s) if n+l+s else None
        rows[key]=dict(year=key[0],population=key[1],line=key[2],quad=key[3],n=n,v=v,w=w,
                       l=l,s=s,distance=int(q['distance']),other_share=a)
    assert [sum(x[k] for x in rows.values()) for k in ('v','w','l','s')]==[5081,8723,4463,29]
    return rows


def sigmoid(t):
    return np.where(t>=0,1/(1+np.exp(-np.maximum(t,-680))),np.exp(np.minimum(t,680))/(1+np.exp(np.minimum(t,680))))


def newton_fit(training,test,feats,lam=10,balanced=False):
    if not feats:
        Ns=sum(q['n'] for q in training);Vs=sum(q['v'] for q in training)
        if balanced:
            by=defaultdict(lambda:[0,0])
            for x in training:by[x['population']][0]+=x['n'];by[x['population']][1]+=x['v']
            ns=np.mean([x[0] for x in by.values()])
            Ns=sum(ns for _ in by);Vs=sum(ns*vs/n for n,vs in by.values())
        return np.full(len(test),np.clip((Vs+0.5)/(Ns+1),.01,.99))
    n=np.array([q['n'] for q in training],float)
    v=np.array([q['v'] for q in training],float)
    if balanced:
        sums=defaultdict(float)
        for q in training:sums[q['population']]+=q['n']
        mn=float(np.mean(list(sums.values())))
        w=np.array([mn/sums[q['population']] for q in training])
        n=n*w;v=v*w
    raw=np.array([[q[k] for k in feats] for q in training],float)
    obs=np.array([[q[k] for k in feats] for q in test],float)
    mu=np.average(raw,weights=n,axis=0)
    sd=np.sqrt(np.average((raw-mu)**2,weights=n,axis=0))
    sd=np.where(sd>1e-12,sd,1.0)
    x=(raw-mu)/sd
    xt=(obs-mu)/sd
    X=np.c_[np.ones(len(n)),x]
    Z=np.c_[np.ones(len(test)),xt]
    vtotal=v.sum();ntotal=n.sum();p=(vtotal+0.5)/(ntotal+1)
    theta=np.r_[np.log(p/(1-p)),np.zeros(len(feats))]
    def obj_grad_hess(t):
        logits=X@t
        proba=sigmoid(logits)
        intercept=sigmoid(np.array([t[0]]))[0]
        obj=float(np.sum(n*np.logaddexp(0,logits)-v*logits)+np.logaddexp(0,t[0])-.5*t[0]+.5*lam*sum(t[1:]**2))
        g=X.T@(n*proba-v)
        g[0]+=intercept-0.5
        g[1:]+=lam*t[1:]
        hess=X.T@(X*(n*proba*(1-proba))[:,None])
        hess[0,0]+=intercept*(1-intercept)
        hess[1:,1:]+=np.eye(len(feats))*lam
        return obj,g,hess
    for i in range(70):
        obj,g,h=obj_grad_hess(theta)
        if max(abs(g))<3e-9:break
        step=np.linalg.solve(h,g)
        scale=1.0
        while scale>=2**(-28):
            proposed=theta-scale*step
            if obj_grad_hess(proposed)[0]<=obj+1e-11:break
            scale/=2
        theta=proposed
    _,g,_=obj_grad_hess(theta)
    if max(abs(g))>2e-6:raise ValueError('INDEPENDENT_NEWTON_NONCONVERGENCE '+str(g))
    return np.clip(sigmoid(Z@theta),.01,.99)


def site_loss(rows,pred):
    arr=np.asarray(pred)
    n=np.array([q['n'] for q in rows]);v=np.array([q['v'] for q in rows])
    return float(sum(-v*np.log(arr)-(n-v)*np.log1p(-arr))/sum(n))


def verify(zipfile_path,prior_json,result_json):
    original=parse_actual_original(zipfile_path)
    raw=Path(prior_json).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==PARENT_SHA
    prior=json.loads(raw)
    result=json.loads(Path(result_json).read_bytes())
    selected=prior['selected_latest']
    bysites={}
    for z in selected:
        name=z['population'];year=z['year']
        population=[q for q in original.values() if q['population']==name and q['year']==year and q['n']>0 and q['distance']>=0]
        expected_keys={(x['year'],x['population'],x['line'],x['quad']) for x in prior['predictions'] if x['population']==name and x['year']==year}
        test=[original[k] for k in expected_keys]
        # exact fixed tests and source identity, but independently parse all source fields.
        assert set((q['year'],q['population'],q['line'],q['quad']) for q in test)==expected_keys
        assert len(test)==z['n_scored_quadrats'] and sum(q['n'] for q in test)==z['n_scored_plants']
        bysites[name]=dict(year=year,train_pool=population,test=test)
    assert len(bysites)==6 and sum(len(z['test']) for z in bysites.values())==168
    RPred={(x['year'],x['population'],x['line'],x['quad']):x for x in result['predictions']}
    RScores={x['population']:x for x in result['site_scores']}
    RSens={(x['population'],x['scenario']):x for x in result['sensitivity_per_site']}
    nchecks=0;maxdiff=0
    def chk(got,want,where,tol=1e-7):
        nonlocal nchecks,maxdiff
        diff=abs(float(got)-float(want));maxdiff=max(maxdiff,diff);nchecks+=1
        assert diff<=tol, (where,got,want,diff)
    for heldout,d in sorted(bysites.items()):
        train=[x for site,data in bysites.items() if site!=heldout for x in data['train_pool']]
        test=d['test']
        assert all(x['population']!=heldout for x in train)
        for lam,label,balanced in [(10,'main',False),(1,'ridge1',False),(100,'ridge100',False),(10,'balanced_ridge10',True)]:
            for model,features in FEATURES.items():
                pred=newton_fit(train,test,features,lam=lam,balanced=balanced)
                ll=site_loss(test,pred)
                expected=(RScores[heldout]['loss'][model] if label=='main' else RSens[heldout,label]['loss'][model])
                chk(ll,expected,('loss',heldout,lam,label,model),tol=1e-7)
                if label=='main':
                    for x,pi in zip(test,pred):
                        key=x['year'],x['population'],x['line'],x['quad']
                        chk(pi,RPred[key]['p_outsite_'+model],('prediction',key,model),tol=1e-7)
        for x in test:
            key=x['year'],x['population'],x['line'],x['quad']
            p=RPred[key]
            for f in ('n','v','w','l','s','distance'):
                assert x[f]==p[f],('original_source_mismatch',key,f)
                nchecks+=1
    assert nchecks>=6*4*4+168*4+168*6
    record=dict(status='PASS',independent_method='csv.DictReader directly on original Mac-CR TXT and independent analytic Newton Hessian, no imports of transfer or parent solver',
                source_original_md5=ZIP_MD5,source_original_sha256=ZIP_SHA,old_fixed_result_sha256=PARENT_SHA,
                compared_newton_results=nchecks,max_abs_difference=maxdiff,
                scoped_named_sites=6,original_rows=512,source_test_quadrats=168,source_test_M_plants=6692,
                calibrated_p_values=False,no_new_field_population_data=True,
                site_name_leakage=False,new_external_geographic_site=False)
    return record

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--archive',required=True);a.add_argument('--parent',required=True);a.add_argument('--result',required=True);a.add_argument('--out',default='results/independent_verification.json');q=a.parse_args()
    r=verify(q.archive,q.parent,q.result)
    Path(q.out).parent.mkdir(parents=True,exist_ok=True)
    Path(q.out).write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(json.dumps(r,ensure_ascii=False))
