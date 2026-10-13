#!/usr/bin/env python3
"""Source-authenticated, leave-one-recorded-field-line-out environmental prediction.

Exploratory historical source reanalysis. Flower-colour associations are non-causal.
"""
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
from scipy.optimize import minimize
from scipy.special import expit

ORIGINAL_ZIP_MD5 = '008ec196aa5acb769dbfd36e835e32ba'
ORIGINAL_ZIP_SHA256 = '4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
FREQ_SHA256 = 'a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
MODELS = {'BASE': [], 'DIST': ['distance'], 'CONGEN': ['other_share'], 'BOTH': ['distance', 'other_share']}


def load_original(archive):
    raw = Path(archive).read_bytes()
    if len(raw)!=127502 or hashlib.md5(raw).hexdigest()!=ORIGINAL_ZIP_MD5 or hashlib.sha256(raw).hexdigest()!=ORIGINAL_ZIP_SHA256:
        raise ValueError('SOURCE_ARCHIVE_HASH_OR_LENGTH_MISMATCH')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip() is not None: raise ValueError('SOURCE_ZIP_CRC_FAIL')
        freq=z.read('data/frequency.txt')
    if hashlib.sha256(freq).hexdigest()!=FREQ_SHA256:
        raise ValueError('SOURCE_MEMBER_SHA256_MISMATCH')
    dr=list(csv.DictReader(io.StringIO(freq.decode('utf-8-sig').replace('\r', '\n')), delimiter='\t'))
    if len(dr)!=512:raise ValueError('EXPECTED_512_ORIGINAL_QUADRATS')
    required={'year','population','line','quad','distance','v','w','l','s'}
    if any(not required.issubset(x) for x in dr):raise ValueError('SOURCE_COLUMNS_MISSING')
    keys=set();rows=[]
    for x in dr:
        a={k:int(x[k]) for k in ['line','distance','v','w','l','s']}
        for k in ('year','population','quad'):a[k]=x[k]
        k=(a['year'],a['population'],a['line'],a['quad'])
        if k in keys: raise ValueError('DUPLICATE_FIELD_QUADRAT_KEY')
        keys.add(k)
        if any(a[k]<0 for k in ('v','w','l','s')):raise ValueError('NEGATIVE_MORPH_OR_SPECIES_COUNT')
        a['n']=a['v']+a['w']
        a['other_share']=(a['l']+a['s'])/(a['n']+a['l']+a['s']) if a['n']+a['l']+a['s']>0 else None
        rows.append(a)
    if (sum(x['v'] for x in rows),sum(x['w'] for x in rows),sum(x['l'] for x in rows),sum(x['s'] for x in rows))!=(5081,8723,4463,29):
        raise ValueError('WRONG_COLOUR_DENOMINATORS')
    return rows


def ridge_objective_gradient(theta,X,n,v,lam=10.0):
    """Negative half-pseudocount penalized binomial log likelihood; analytic gradient."""
    theta=np.asarray(theta,dtype=float)
    z=theta[0]+X@theta[1:]
    # Jeffreys Beta(0.5,0.5) pseudo-counts at feature zero.
    objective=(np.dot(n,np.logaddexp(0,z))-np.dot(v,z)
        +np.logaddexp(0,theta[0])-.5*theta[0]
        +.5*lam*np.dot(theta[1:],theta[1:]))
    residue=n*expit(z)-v
    grad=np.empty_like(theta)
    grad[0]=np.sum(residue)+expit(theta[0])-.5
    grad[1:]=X.T@residue+lam*theta[1:]
    return float(objective),grad



def hessian_newton_polish(theta,X,n,v,lam=10.0):
    """Refine the same convex binomial objective after numerical BFGS precision loss.

    No changes to source selection, features, penalties, evaluation, or heldout data.
    The Hessian is strictly positive-definite from regularization plus pseudo-prior.
    """
    theta=np.asarray(theta,dtype=float).copy()
    design=np.column_stack([np.ones(len(n)),X])
    for _ in range(40):
        original,grad=ridge_objective_gradient(theta,X,n,v,lam)
        if np.max(np.abs(grad))<=1e-8:
            break
        q=expit(design@theta)
        hessian=design.T@(design*(n*q*(1-q))[:,None])
        hessian[0,0]+=expit(theta[0])*(1-expit(theta[0]))
        hessian[1:,1:]+=np.eye(len(theta)-1)*lam
        step=np.linalg.solve(hessian,grad)
        alpha=1.0
        while alpha>=2**(-24) and ridge_objective_gradient(theta-alpha*step,X,n,v,lam)[0]>original+1e-11:
            alpha*=0.5
        theta-=alpha*step
    _,grad=ridge_objective_gradient(theta,X,n,v,lam)
    if not np.all(np.isfinite(theta)) or np.max(np.abs(grad))>1e-6:
        raise RuntimeError('NEWTON_SAME_OBJECTIVE_REFINEMENT_FAILED')
    return theta


def fit_binomial(train,test,features,lam=10.0):
    if not train or not test:raise ValueError('EMPTY_TRAIN_OR_TEST')
    vn=np.array([x['v'] for x in train],dtype=float)
    n=np.array([x['n'] for x in train],dtype=float)
    if np.any(n<=0) or vn.sum()<5 or (n-vn).sum()<5:raise ValueError('TRAIN_MORPH_SUPPORT_FAIL')
    base=(vn.sum()+.5)/(n.sum()+1)
    if not features:return np.full(len(test),np.clip(base,.01,.99))
    X=np.array([[x[f] for f in features] for x in train],dtype=float)
    Xt=np.array([[x[f] for f in features] for x in test],dtype=float)
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(Xt))):raise ValueError('NONFINITE_PREDICTOR')
    mu=np.average(X,weights=n,axis=0)
    sd=np.sqrt(np.average((X-mu)**2,weights=n,axis=0))
    identifiable=sd>1e-12
    X=(X-mu)/np.where(identifiable,sd,1)
    Xt=(Xt-mu)/np.where(identifiable,sd,1)
    X[:,~identifiable]=0
    Xt[:,~identifiable]=0
    initial=np.r_[math.log(base/(1-base)),np.zeros(len(features))]
    # BFGS can signal precision loss on high-count folds even when the fitted
    # objective is correct. An analytic Hessian Newton polish refines the SAME
    # fixed penalized likelihood, using training quadrats only, across runners.
    res=minimize(ridge_objective_gradient,initial,args=(X,n,vn,lam),method='BFGS',jac=True,
                 options={'gtol':1e-8,'maxiter':1000})
    if not np.all(np.isfinite(res.x)):
        raise RuntimeError('NONFINITE_BFGS_SOLUTION')
    theta=hessian_newton_polish(res.x,X,n,vn,lam)
    q=expit(theta[0]+Xt@theta[1:])
    if not np.all(np.isfinite(q)):raise ValueError('NONFINITE_HELDOUT_PROBABILITY')
    return np.clip(q,.01,.99)


def score(test,pred):
    rows=list(test); q=np.asarray(pred,dtype=float)
    if len(rows)!=len(q) or not len(rows):raise ValueError('SCORE_SHAPE_MISMATCH')
    n=np.array([r['n'] for r in rows],dtype=float)
    v=np.array([r['v'] for r in rows],dtype=float)
    if np.any(n<=0) or np.any(q<=0) or np.any(q>=1):raise ValueError('INVALID_SCORE_SUPPORT')
    return {
      'logloss_per_plant': float(np.sum(-v*np.log(q)-(n-v)*np.log1p(-q))/np.sum(n)),
      'fraction_brier_equal_plot': float(np.mean((v/n-q)**2)),
      'n_test_plants':int(n.sum()),'n_test_quadrats':len(rows),
    }


def derive_cohort(rows):
    bysite=defaultdict(list)
    for x in rows:bysite[(x['year'],x['population'])].append(x)
    eligible=[]; held={}
    for (year,pop),r in sorted(bysite.items()):
        base=[x for x in r if x['n']>0 and x['distance']>=0]
        covered={line:q for line in sorted({x['line'] for x in base}) for q in [[x for x in base if x['line']==line]]}
        qualifying=[line for line,q in covered.items() if len(q)>=2 and sum(z['n'] for z in q)>=10]
        n=sum(x['n'] for x in base);v=sum(x['v'] for x in base);w=sum(x['w'] for x in base)
        plus=[x for x in base if x['l']+x['s']>0]
        zero=[x for x in base if x['l']+x['s']==0]
        reasons=[]
        if n<50 or v<5 or w<5:reasons.append('LOW_OR_MONOMORPH_M_COUNTS')
        if len(qualifying)<3:reasons.append('FEWER_THAN_THREE_MULTIPLOT_LINES')
        if len({x['distance'] for x in base})<2:reasons.append('NO_WATER_DISTANCE_VARIATION')
        if len(plus)<3 or len({x['line'] for x in plus})<2 or len(zero)<3 or len({x['line'] for x in zero})<2:
            reasons.append('NO_CONGENER_COVARIATE_OVERLAP')
        scored=[]; fold_holds={}
        for line in qualifying:
            train=[x for x in base if x['line']!=line]
            tn=sum(x['n'] for x in train);tv=sum(x['v'] for x in train)
            if tn>=50 and tv>=5 and tn-tv>=5:
                scored.append(line)
            else:fold_holds[str(line)]='TRAIN_MORPH_SUPPORT_FAIL'
        if len(scored)<3:reasons.append('FEWER_THAN_THREE_TRAINABLE_LEFT_OUT_LINES')
        info=dict(year=year,population=pop,siteyear=f'{year}|{pop}',n_original_quadrats=len(r),
                  n_M_distance_known_quadrats=len(base),n_M_eligible_distance_known=n, V_distance_known=v,W_distance_known=w,
                  n_qualifying_lines=len(qualifying),n_eligible_folds=len(scored),n_other_present=len(plus),n_other_absent=len(zero),
                  status='ELIGIBLE' if not reasons else 'HOLD',hold_reasons=reasons,
                  heldout_lines=scored,fold_holds=fold_holds)
        if reasons:held[info['siteyear']]=info
        else:eligible.append(info)
    return eligible,held


def evaluate(rows):
    eligible,held=derive_cohort(rows)
    bysite=defaultdict(list)
    for x in rows:
        if x['n']>0 and x['distance']>=0:bysite[(x['year'],x['population'])].append(x)
    predictions=[];site_scores=[];fold_ledger=[]
    for site in eligible:
        year,pop=site['year'],site['population'];data=bysite[(year,pop)]; site_pred=[]
        for line in site['heldout_lines']:
            train=[x for x in data if x['line']!=line]
            test=[x for x in data if x['line']==line]
            assert not set((z['year'],z['population'],z['line'],z['quad']) for z in train).intersection(
                       set((z['year'],z['population'],z['line'],z['quad']) for z in test))
            p={m:fit_binomial(train,test,features) for m,features in MODELS.items()}
            fold_ledger.append(dict(siteyear=site['siteyear'],year=year,population=pop,line=line,
                                    train_n=sum(x['n'] for x in train),train_v=sum(x['v'] for x in train),
                                    train_w=sum(x['w'] for x in train),test_n=sum(x['n'] for x in test),
                                    test_quadrats=len(test),
                                    train_other_share_distinct=len({x['other_share'] for x in train}),
                                    train_water_distance_distinct=len({x['distance'] for x in train})))
            for j,x in enumerate(test):
                rec={k:x[k] for k in ('year','population','line','quad','distance','v','w','l','s','n','other_share')}
                for m in MODELS: rec['p_'+m]=float(p[m][j])
                site_pred.append(rec)
        if not site_pred:raise RuntimeError('QUALIFIED_SITE_WITH_NO_PREDICTIONS')
        stats={m:score(site_pred,[p['p_'+m] for p in site_pred]) for m in MODELS}
        site_scores.append(dict(siteyear=site['siteyear'],year=year,population=pop,n_folds=len(site['heldout_lines']),
                               n_scored_quadrats=len(site_pred),n_scored_plants=sum(p['n'] for p in site_pred),
                               **{f'logloss_{m}':stats[m]['logloss_per_plant'] for m in MODELS},
                               **{f'quadrat_brier_{m}':stats[m]['fraction_brier_equal_plot'] for m in MODELS},
                               gain_dist=stats['BASE']['logloss_per_plant']-stats['DIST']['logloss_per_plant'],
                               gain_congen=stats['BASE']['logloss_per_plant']-stats['CONGEN']['logloss_per_plant'],
                               gain_both=stats['BASE']['logloss_per_plant']-stats['BOTH']['logloss_per_plant'],
                               extra_congen_beyond_distance=stats['DIST']['logloss_per_plant']-stats['BOTH']['logloss_per_plant'],
                               extra_distance_beyond_congen=stats['CONGEN']['logloss_per_plant']-stats['BOTH']['logloss_per_plant']))
        predictions+=site_pred
    latest={}
    for x in sorted(site_scores,key=lambda z:z['year']):latest[x['population']]=x
    selected=[latest[p] for p in sorted(latest)]
    def aggregate(items):
        cols=['logloss_BASE','logloss_DIST','logloss_CONGEN','logloss_BOTH',
              'quadrat_brier_BASE','quadrat_brier_DIST','quadrat_brier_CONGEN','quadrat_brier_BOTH',
              'gain_dist','gain_congen','gain_both','extra_congen_beyond_distance','extra_distance_beyond_congen']
        return dict(n_sites=len(items),n_plants_scored=sum(r['n_scored_plants'] for r in items),
             n_lines_heldout=sum(r['n_folds'] for r in items),
             mean={k:float(np.mean([r[k] for r in items])) for k in cols},
             n_sites_positive={k:sum(r[k]>0 for r in items) for k in cols if k.startswith(('gain_','extra_'))})
    main=aggregate(selected)
    extended=aggregate(site_scores)
    loo=[]
    for site in selected:
        a=aggregate([x for x in selected if x['population']!=site['population']])
        loo.append(dict(dropped=site['population'],remaining_n=a['n_sites'],**a['mean']))
    qa=dict(n_source_quadrats=len(rows),original_M_v=sum(x['v'] for x in rows),original_M_w=sum(x['w'] for x in rows),
            unrelated_other_species_L=sum(x['l'] for x in rows),unrelated_other_species_S=sum(x['s'] for x in rows),
            n_empty_M_quadrats=sum(x['n']==0 for x in rows),n_missing_distance=sum(x['distance']<0 for x in rows),
            n_missing_distance_with_M=sum(x['distance']<0 and x['n']>0 for x in rows),
            n_eligible_siteyears=len(eligible),n_held_siteyears=len(held),
            no_new_observational_data=True,not_causal=True,not_new_geographic_site_prediction=True)
    return dict(schema='fcp_sisyrinchium_heldout_field_line_distance_congener_prediction_v1',qa=qa,
                model_specification={'lambda_ridge':10,'clipping':[.01,.99],
                 'training':'all other lines within SAME named siteyear',
                 'primary_sites':'latest eligible siteyear per original population name',
                 'primary_score':'binomial negative log likelihood per M plant, equal named-site mean',
                 'predictors':MODELS,'p_values_computed':False,'source_exposed_not_preregistered':True},
                eligible_siteyears=eligible,held_siteyears=held,
                all_siteyears=[*eligible,*held.values()],site_scores=site_scores,
                selected_latest=selected, main=main,all_years_sensitivity=extended,
                leave_one_named_site_out=loo,predictions=predictions,fold_ledger=fold_ledger)


def write_outputs(result,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    for key,name in [('all_siteyears','siteyear_support.csv'),('site_scores','siteyear_predictive_scores.csv'),
                     ('predictions','heldout_quadrat_predictions.csv'),('fold_ledger','fold_support.csv'),
                     ('leave_one_named_site_out','site_deletion_sensitivity.csv')]:
        rows=result[key]
        cols=list(rows[0]) if rows else []
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore',lineterminator='\n');w.writeheader()
            for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in r.items()})
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default='source/Sisyrinchium_data.zip');p.add_argument('--outdir',default='results');a=p.parse_args()
    result=evaluate(load_original(a.source));write_outputs(result,a.outdir)
    print(json.dumps({'qa':result['qa'],'main':result['main'],'extended':result['all_years_sensitivity']},ensure_ascii=False))
