#!/usr/bin/env python3
"""Across-named-population shoreline-distance slope portability using original field counts.

Posthoc, source-exposed exploratory test. Only target-site intercept is calibrated
from other target field lines; transferred slope is learned in OTHER named sites.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.special import expit

SOURCE_MD5='008ec196aa5acb769dbfd36e835e32ba'
SOURCE_SHA256='4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
MEMBER_SHA256='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
LAMBDA=10.0
CLIP=.01

def load_original(path):
    b=Path(path).read_bytes()
    if len(b)!=127502 or hashlib.md5(b).hexdigest()!=SOURCE_MD5 or hashlib.sha256(b).hexdigest()!=SOURCE_SHA256:raise ValueError('SOURCE_BYTES_INVALID')
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        if z.testzip() is not None:raise ValueError('SOURCE_ZIP_CRC_INVALID')
        raw=z.read('data/frequency.txt')
    if hashlib.sha256(raw).hexdigest()!=MEMBER_SHA256:raise ValueError('SOURCE_MEMBER_CHANGED')
    rr=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig').replace('\r','\n')),delimiter='\t'))
    if len(rr)!=512: raise ValueError('WRONG_SOURCE_ROW_COUNT')
    seen=set();out=[]
    for x in rr:
        y={k:x[k] for k in ('year','population','quad')}
        y.update({k:int(x[k]) for k in ('line','distance','v','w','l','s')})
        key=(y['year'],y['population'],y['line'],y['quad'])
        if key in seen:raise ValueError('DUPLICATE_SOURCE_FIELD_KEY')
        seen.add(key)
        if min(y[k] for k in ('v','w','l','s'))<0:raise ValueError('INVALID_NEGATIVE_COUNT')
        y['n']=y['v']+y['w'];y['siteyear']=y['year']+'|'+y['population']
        out.append(y)
    if tuple(sum(z[k] for z in out) for k in ('v','w','l','s'))!=(5081,8723,4463,29):raise ValueError('INVALID_MORPH_OR_SPECIES_CENSUS')
    return out

def baseline(n,v):
    if not (0<=v<=n):raise ValueError('INVALID_TARGET_M_COUNTS')
    return (v+.5)/(n+1)

def negative_ridge(theta,X,n,v):
    z=theta[0]+theta[1]*X
    val=np.dot(n,np.logaddexp(0,z))-np.dot(v,z) +np.logaddexp(0,theta[0])-.5*theta[0]+.5*LAMBDA*theta[1]**2
    r=n*expit(z)-v
    grad=np.array([r.sum()+expit(theta[0])-.5,np.dot(r,X)+LAMBDA*theta[1]])
    return float(val),grad

def hessian(theta,X,n):
    z=theta[0]+theta[1]*X
    w=n*expit(z)*(1-expit(z))
    return np.array([[w.sum()+expit(theta[0])*(1-expit(theta[0])),np.dot(w,X)],
                     [np.dot(w,X),np.dot(w,X*X)+LAMBDA]])

def fit_source_slope(rows):
    """Logit coefficient per ONE raw source distance, ridge on standardized x."""
    if any(z['distance']<0 or z['n']<=0 for z in rows):raise ValueError('SOURCE_SLOPE_NEEDS_DISTANCE_AND_M')
    n=np.array([z['n'] for z in rows],float);v=np.array([z['v'] for z in rows],float);d=np.array([z['distance'] for z in rows],float)
    if v.sum()<5 or (n-v).sum()<5:raise ValueError('SOURCE_SITE_BOTH_M_REQUIRED')
    mu=np.average(d,weights=n);sd=np.sqrt(np.average((d-mu)**2,weights=n))
    if sd<1e-12:raise ValueError('SOURCE_NO_DISTANCE_SPREAD')
    X=(d-mu)/sd
    first=np.array([math.log(baseline(n.sum(),v.sum())/(1-baseline(n.sum(),v.sum()))),0.0])
    opt=minimize(negative_ridge,first,args=(X,n,v),method='BFGS',jac=True,options={'gtol':1e-8,'maxiter':1000})
    theta=np.asarray(opt.x,float)
    # Identical positive-definite objective, exact-Hessian numerical polish.
    for _ in range(60):
        old,g=negative_ridge(theta,X,n,v)
        if max(abs(g))<=1e-8:break
        step=np.linalg.solve(hessian(theta,X,n),g)
        eta=1.0
        while eta>2**(-25) and negative_ridge(theta-eta*step,X,n,v)[0]>old+1e-11:eta*=.5
        theta-=eta*step
    val,g=negative_ridge(theta,X,n,v)
    if max(abs(g))>1e-6 or not np.isfinite(theta).all():raise RuntimeError('Ridge binomial gradient not converged')
    return {'beta_per_distance':float(theta[1]/sd),'beta_standardized':float(theta[1]),
            'intercept_centered':float(theta[0]),'weighted_distance_center':float(mu),
            'weighted_distance_sd':float(sd),'n_M':int(n.sum()),'n_quadrats':len(rows),
            'source_v':int(v.sum()),'source_w':int(n.sum()-v.sum()),
            'source_distance_min':int(min(d)),'source_distance_max':int(max(d)),
            'gradient_linf':float(max(abs(g)))}

def calibrate_target_intercept(train,borrowed_beta):
    if any(z['distance']<0 or z['n']<=0 for z in train):raise ValueError('UNSUPPORTED_TARGET_TRAIN')
    n=np.array([z['n'] for z in train],float);v=np.array([z['v'] for z in train],float);d=np.array([z['distance'] for z in train],float)
    mu=float(np.average(d,weights=n))
    z=borrowed_beta*(d-mu)
    def eq(alpha):return float(np.dot(n,expit(alpha+z))-v.sum()+expit(alpha)-.5)
    alpha=brentq(eq,-80,80,xtol=1e-13)
    if abs(eq(alpha))>1e-7:raise RuntimeError('NO_INTERCEPT_ROOT')
    return float(alpha),mu

def nll(v,w,q):
    q=float(np.clip(q,CLIP,1-CLIP))
    return -v*math.log(q)-w*math.log1p(-q)

def analysis(original_zip,parent_result,parent_predictions,outdir):
    raw=load_original(original_zip)
    frozen=json.loads(Path(parent_result).read_text())
    selected=frozen['selected_latest']
    sites=sorted(x['siteyear'] for x in selected)
    if sites!=sorted(['Y13|banpaku','Y13|nature','Y11|ootaka','Y11|ootaka2','Y13|shiga','Y11|toyama']):raise ValueError('PARENT_SITE_SUPPORT_CHANGED')
    parent_csv=list(csv.DictReader(Path(parent_predictions).open(newline='')))
    if len(parent_csv)!=185:raise ValueError('PARENT_HELDOUT_ROW_COUNT_CHANGED')
    bysite=defaultdict(list)
    for y in raw:
        if y['siteyear'] in sites and y['n']>0 and y['distance']>=0:bysite[y['siteyear']].append(y)
    all_source={s:fit_source_slope(bysite[s]) for s in sites}
    included_siteyears={s:next(z for z in selected if z['siteyear']==s) for s in sites}
    parent_map={}
    for r in parent_csv:
        key=(r['year'],r['population'],int(r['line']),r['quad'])
        if key in parent_map:raise ValueError('PARENT_DUPLICATE_PREDICTION')
        parent_map[key]=r
    new_preds=[];folds=[];summary=[]
    for site in sites:
        others=[s for s in sites if s!=site]
        if site in others or len(others)!=5:raise AssertionError('TARGET_IN_SOURCE_SLOPE')
        betas=[all_source[z]['beta_per_distance'] for z in others]
        borrowed=float(np.median(betas));borrowed_mean=float(np.mean(betas))
        rows=bysite[site]
        chosen=sorted(set(int(z['line']) for z in rows if (z['year'],z['population'],int(z['line']),z['quad']) in parent_map))
        if len(chosen)!=included_siteyears[site]['n_folds']:raise ValueError('ELIGIBLE_FOLD_COUNT_CHANGED')
        site_preds=[];fold_signs=[]
        for h in chosen:
            train=[z for z in rows if z['line']!=h]
            test=[z for z in rows if z['line']==h and (z['year'],z['population'],z['line'],z['quad']) in parent_map]
            if not train or not test:raise ValueError('MISSING_FOLD_ROWS')
            if {z['line'] for z in train}&{h}:raise AssertionError('LEAKED_HELDOUT_LINE')
            target_fit=fit_source_slope(train)
            fold_signs.append(target_fit['beta_per_distance'])
            a,mu=calibrate_target_intercept(train,borrowed)
            a2,mu2=calibrate_target_intercept(train,borrowed_mean)
            if mu!=mu2:raise AssertionError('DIFFERENT_CALIBRATION_DATA')
            fold={'siteyear':site,'heldout_line':h,'source_sites':others,
                  'source_n_sites':len(others),'train_n_M':sum(z['n'] for z in train),
                  'train_n_lines':len(set(z['line'] for z in train)),
                  'test_n_M':sum(z['n'] for z in test),'test_n_quadrats':len(test),
                  'borrow_median_beta_per_distance':borrowed,'borrow_mean_beta_per_distance':borrowed_mean,
                  'target_training_local_beta':target_fit['beta_per_distance'],
                  'target_intercept_calibration':a,'target_intercept_calibration_mean':a2,
                  'target_train_distance_center':mu}
            folds.append(fold)
            for z in test:
                key=(z['year'],z['population'],z['line'],z['quad'])
                r=parent_map[key]
                for k in ['v','w','n','distance','l','s']:
                    if int(r[k])!=z[k]:raise ValueError('PARENT_QUADRAT_SOURCE_MISMATCH_'+k)
                inherited={m:float(r['p_'+m]) for m in ['BASE','DIST','CONGEN','BOTH']}
                borrowed_pred=float(np.clip(expit(a+borrowed*(z['distance']-mu)),CLIP,1-CLIP))
                borrowed_mean_pred=float(np.clip(expit(a2+borrowed_mean*(z['distance']-mu)),CLIP,1-CLIP))
                site_preds.append({**z,'p_BASE':inherited['BASE'],'p_LOCAL_DIST':inherited['DIST'],
                    'p_BORROW_MEDIAN':borrowed_pred,'p_BORROW_MEAN':borrowed_mean_pred,
                    'borrowed_site_median_beta':borrowed,'borrowed_site_mean_beta':borrowed_mean,
                    'heldout_line':h})
        if len(site_preds)!=included_siteyears[site]['n_scored_quadrats'] or sum(z['n'] for z in site_preds)!=included_siteyears[site]['n_scored_plants']:
            raise ValueError('FROZEN_SITE_SCORING_SUPPORT_CHANGED')
        if len({(z['year'],z['population'],z['line'],z['quad']) for z in site_preds})!=len(site_preds):raise ValueError('DUPLICATE_PREDICTIONS')
        row={'siteyear':site,'n_lines':len(chosen),'n_quadrats':len(site_preds),'n_M':sum(z['n'] for z in site_preds),
             'full_site_slope_raw':all_source[site]['beta_per_distance'],
             'source_median_raw_other_sites':borrowed,
             'source_mean_raw_other_sites':borrowed_mean,
             'full_site_sign':int(np.sign(all_source[site]['beta_per_distance'])),
             'within_site_fold_sign_same':int(sum(np.sign(s)==np.sign(all_source[site]['beta_per_distance']) for s in fold_signs)),
             'n_target_training_folds':len(fold_signs)}
        for model in ['BASE','LOCAL_DIST','BORROW_MEDIAN','BORROW_MEAN']:
            row['logloss_'+model]=sum(nll(z['v'],z['w'],z['p_'+model]) for z in site_preds)/row['n_M']
        row['gain_local_vs_base']=row['logloss_BASE']-row['logloss_LOCAL_DIST']
        row['gain_borrow_vs_base']=row['logloss_BASE']-row['logloss_BORROW_MEDIAN']
        row['gain_borrow_vs_local']=row['logloss_LOCAL_DIST']-row['logloss_BORROW_MEDIAN']
        row['gain_mean_vs_base']=row['logloss_BASE']-row['logloss_BORROW_MEAN']
        row['gain_mean_vs_local']=row['logloss_LOCAL_DIST']-row['logloss_BORROW_MEAN']
        if abs(row['logloss_BASE']-included_siteyears[site]['logloss_BASE'])>1e-10 or abs(row['logloss_LOCAL_DIST']-included_siteyears[site]['logloss_DIST'])>1e-10:
            raise ValueError('FROZEN_PARENT_MODEL_SCORE_DISAGREES')
        summary.append(row)
        new_preds.extend(site_preds)
    if (len(folds),len(new_preds),sum(z['n'] for z in new_preds))!=(32,168,6692):raise ValueError('UNEXPECTED_TOTAL_SCORING_UNITS')
    mean={k:sum(row[k] for row in summary)/len(summary) for k in summary[0] if k.startswith(('logloss_','gain_'))}
    grand={'mean':mean,'n_sites':6,'n_folds':32,'n_quadrats':168,'n_M':6692,
           'n_borrow_positive_vs_base':sum(x['gain_borrow_vs_base']>0 for x in summary),
           'n_borrow_positive_vs_local':sum(x['gain_borrow_vs_local']>0 for x in summary),
           'n_local_positive_vs_base':sum(x['gain_local_vs_base']>0 for x in summary),
           'n_source_sites_negative_gradient':sum(x['full_site_sign']<0 for x in summary),
           'n_source_sites_positive_gradient':sum(x['full_site_sign']>0 for x in summary),
           'n_training_fold_sign_agree':sum(x['within_site_fold_sign_same'] for x in summary)}
    sensitivity=[]
    for omitted in sites:
        x=[r for r in summary if r['siteyear']!=omitted]
        sensitivity.append({'deleted_siteyear':omitted,'n_sites':len(x),
            'gain_borrow_vs_base':sum(z['gain_borrow_vs_base'] for z in x)/len(x),
            'gain_borrow_vs_local':sum(z['gain_borrow_vs_local'] for z in x)/len(x),
            'gain_local_vs_base':sum(z['gain_local_vs_base'] for z in x)/len(x)})
    response={'schema':'fcp_sisyrinchium_shoreline_gradient_cross_population_transport_v1',
              'original_source':{'doi':'10.5061/dryad.j8q34','original_zip_md5':SOURCE_MD5,
                     'original_zip_sha256':SOURCE_SHA256,'frequency_member_sha256':MEMBER_SHA256,
                     'n_original_quadrats':512,'M_purple':5081,'M_white':8723,'other_l':4463,'other_s':29},
              'scope':'SOURCE_EXPOSED_EXPLORATORY_NOT_PREREGISTERED',
              'method':'Leave one source population out for slope; calibrate intercept ONLY on other target-site lines; score same original heldout lines',
              'slope_estimator':'Ridge lambda10 per site on M-count-weighted standardized distance, backtransform to raw recorded distance, 5-source-site equal median (primary), mean (sensitivity)',
              'source_fit':all_source,'sites':summary,'main':grand,'folds':folds,'predictions':new_preds,'leave_one_site_out':sensitivity,
              'claims':['A source-transported distance gradient is a predictive hypothesis, not genetic/ecological selection or causal water-distance manipulation',
                  'Within-site slope predictiveness can exist even when a single cross-site slope rule fails',
                  'All six named sites belong to the same pre-exposed 2016 paper, not independent geographic new data',
                  'Transferred slopes exclude target-site data, but a target-site intercept is calibrated from target-site training lines',
                  'Other L/S white species must not be pooled into M-white; original PR185 results unchanged']}
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/'result.json').write_text(json.dumps(response,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    def write_csv(name,rows,keys):
        with (out/name).open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=keys,lineterminator='\n',extrasaction='ignore')
            wr.writeheader();wr.writerows(rows)
    write_csv('named_site_transfer_scores.csv',summary,list(summary[0]))
    write_csv('borrowed_line_predictions.csv',new_preds,[k for k in new_preds[0] if k not in ['siteyear']]+['siteyear'])
    write_csv('outer_site_and_inner_line_fit_ledger.csv',folds,list(folds[0]))
    write_csv('site_deletion_transfer_sensitivity.csv',sensitivity,list(sensitivity[0]))
    print(json.dumps({'status':'PASS','n_sites':grand['n_sites'],'n_target_lines':grand['n_folds'],
          'main':grand,'source_signs':{k:int(np.sign(v['beta_per_distance'])) for k,v in all_source.items()}},ensure_ascii=False))
    return response

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',required=True)
    ap.add_argument('--parent-result',required=True)
    ap.add_argument('--parent-predictions',required=True)
    ap.add_argument('--outdir',required=True)
    z=ap.parse_args()
    analysis(z.source,z.parent_result,z.parent_predictions,z.outdir)
