#!/usr/bin/env python3
"""POSTHOC (not preregistered) source-site distance-scale harmonisation audit.

Compares borrowed logit gradient per one *within-site distance SD* instead of
per one recorded distance unit. Target-site feature SD uses target TRAIN lines
only; the source slope uses only the five OTHER original named sites.
"""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
from scipy.special import expit
from transport import load_original,calibrate_target_intercept,nll,CLIP

def examine(source, primary_result, output):
    d=json.loads(Path(primary_result).read_text())
    origin=load_original(source)
    sites=sorted(d['source_fit'])
    targets={s:[z for z in origin if z['siteyear']==s and z['n']>0 and z['distance']>=0] for s in sites}
    prior={(z['year'],z['population'],z['line'],z['quad']):z for z in d['predictions']}
    result=[];ledger=[]
    for s in sites:
        others=[x for x in sites if x!=s]
        source_dimless=[d['source_fit'][x]['beta_standardized'] for x in others]
        median=float(np.median(source_dimless));mean=float(np.mean(source_dimless))
        qs=[]
        fold_labels=sorted({z['line'] for z in targets[s] if (z['year'],z['population'],z['line'],z['quad']) in prior})
        for line in fold_labels:
            train=[z for z in targets[s] if z['line']!=line]
            test=[z for z in targets[s] if z['line']==line and (z['year'],z['population'],z['line'],z['quad']) in prior]
            w=np.array([z['n'] for z in train],float);distance=np.array([z['distance'] for z in train],float)
            mu=float(np.average(distance,weights=w));sd=float(np.sqrt(np.average((distance-mu)**2,weights=w)))
            if sd<1e-12:raise ValueError('TARGET_TRAIN_DISTANCE_NO_SPREAD')
            bmed=median/sd;bmean=mean/sd
            a,center=calibrate_target_intercept(train,bmed)
            a2,center2=calibrate_target_intercept(train,bmean)
            assert center==center2==mu
            for z in test:
                key=(z['year'],z['population'],z['line'],z['quad'])
                p=prior[key]
                if p['siteyear']!=s or z['v']!=p['v'] or z['w']!=p['w']:raise ValueError('PRIMARY_SUPPORT_DRIFT')
                qmed=float(np.clip(expit(a+bmed*(z['distance']-mu)),CLIP,1-CLIP))
                qmean=float(np.clip(expit(a2+bmean*(z['distance']-mu)),CLIP,1-CLIP))
                row={'siteyear':s,'line':line,'quad':z['quad'],'n':z['n'],'v':z['v'],'w':z['w'],
                     'source_beta_std_median':median,'source_beta_std_mean':mean,
                     'target_train_distance_sd':sd,'source_excludes_target':True,
                     'p_BASE':p['p_BASE'],'p_LOCAL_DIST':p['p_LOCAL_DIST'],
                     'p_BORROW_Z_MEDIAN':qmed,'p_BORROW_Z_MEAN':qmean}
                qs.append(row);ledger.append(row)
        n_total=sum(z['n'] for z in qs)
        ref=next(k for k in d['sites'] if k['siteyear']==s)
        if n_total!=ref['n_M'] or len(qs)!=ref['n_quadrats']:raise AssertionError('DIFFERENT_SUPPORT')
        score={name:sum(nll(z['v'],z['w'],z['p_'+name]) for z in qs)/n_total for name in ['BASE','LOCAL_DIST','BORROW_Z_MEDIAN','BORROW_Z_MEAN']}
        if abs(score['BASE']-ref['logloss_BASE'])>1e-11 or abs(score['LOCAL_DIST']-ref['logloss_LOCAL_DIST'])>1e-11:raise AssertionError('BASE_SUPPORT_DRIFT')
        result.append({'siteyear':s,'n_lines':len(fold_labels),'n_M':n_total,'n_quadrats':len(qs),
          'score':score,'gain_std_median_vs_base':score['BASE']-score['BORROW_Z_MEDIAN'],
          'gain_std_median_vs_local':score['LOCAL_DIST']-score['BORROW_Z_MEDIAN'],
          'gain_std_mean_vs_base':score['BASE']-score['BORROW_Z_MEAN']})
    main={'mean_logloss':{k:sum(row['score'][k] for row in result)/6 for k in result[0]['score']},
          'gain_median_vs_base':sum(row['gain_std_median_vs_base'] for row in result)/6,
          'gain_mean_vs_base':sum(row['gain_std_mean_vs_base'] for row in result)/6,
          'n_gain_median_vs_base':sum(row['gain_std_median_vs_base']>0 for row in result),
          'n_gain_median_vs_local':sum(row['gain_std_median_vs_local']>0 for row in result),
          'n_sites':6,'n_folds':32,'n_quadrats':len(ledger),'n_M':sum(z['n'] for z in ledger)}
    sensitivity=[]
    for omitted in sites:
        rr=[r for r in result if r['siteyear']!=omitted]
        sensitivity.append({'omitted_siteyear':omitted,'gain_std_median_vs_base':sum(r['gain_std_median_vs_base'] for r in rr)/len(rr),
            'gain_std_median_vs_local':sum(r['gain_std_median_vs_local'] for r in rr)/len(rr)})
    result_dict={'schema':'fcp_sisyrinchium_posthoc_standardized_shoreline_transfer_sensitivity_v1',
       'status':'POSTHOC_DIAGNOSTIC_NOT_PRIMARY',
       'explanation':'Dimensionless source-site beta/target TRAIN-site distance SD; slope from 5 other sites; target intercept fit only on other target lines',
       'main':main,'sites':result,'source_sites':sites,'leave_one_named_site_out':sensitivity,
       'hard_limits':['Distance standardisation choice was explored after the primary raw-distance transport failure; not preregistered or confirmatory',
          'Target site training-line distance SD uses predictors only, not target line response',
          'No site-level independent evolutionary experiment, no causal water effect, no plant abundance or colour reconstruction']}
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    (out/'posthoc_distance_scale.json').write_text(json.dumps(result_dict,sort_keys=True,ensure_ascii=False,indent=2)+'\n')
    with (out/'posthoc_distance_scale_predictions.csv').open('w',newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=list(ledger[0]),lineterminator='\n');wr.writeheader();wr.writerows(ledger)
    print(json.dumps({'status':'PASS','type':'POSTHOC_DIAGNOSTIC','main':main},ensure_ascii=False))
    return result_dict

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--primary-result',required=True);p.add_argument('--outdir',required=True);args=p.parse_args()
    examine(args.source,args.primary_result,args.outdir)
