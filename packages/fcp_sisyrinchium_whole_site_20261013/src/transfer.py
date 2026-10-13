#!/usr/bin/env python3
"""Leave-one-WHOLE-field-site-out flower-colour transfer on archived Sisyrinchium.

Retrospective source-exposed original field-quadrat analysis, not geographic foresight.
Fit is exactly the locked predecessor's penalized binomial fit without altering inputs.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path

import numpy as np

FIELDS={'BASE':[],'DIST':['distance'],'CONGEN':['other_share'],'BOTH':['distance','other_share']}
ORIGINAL_MD5='008ec196aa5acb769dbfd36e835e32ba'
ORIGINAL_SHA256='4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5'
ORIGINAL_MEMBER='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'
PARENT_RESULTS_SHA256='7745ab3018eb00a1939b9433effb201e1c9cc37034f2a4c2ae8a588e3c2b23ce'


def parent_module(path):
    spec=importlib.util.spec_from_file_location('fcp_parent_recorded_binomial',str(path))
    p=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p)
    return p


def load_original_verified(archive,parent_predictive):
    raw=Path(archive).read_bytes()
    if len(raw)!=127502 or hashlib.md5(raw).hexdigest()!=ORIGINAL_MD5 or hashlib.sha256(raw).hexdigest()!=ORIGINAL_SHA256:
        raise ValueError('ORIGINAL_SOURCE_IDENTITY_FAILED')
    import zipfile
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('ZIP_CRC_FAILED')
        f=z.read('data/frequency.txt')
        if hashlib.sha256(f).hexdigest()!=ORIGINAL_MEMBER:raise ValueError('INNER_ORIGINAL_SHA_FAILED')
    m=parent_module(parent_predictive)
    return m.load_original(archive)


def model_specs():return dict(FIELDS)


def source_qa(rows):
    return dict(n_quadrats=len(rows),M_v=sum(x['v'] for x in rows),M_w=sum(x['w'] for x in rows),
                L=sum(x['l'] for x in rows),S=sum(x['s'] for x in rows),zero_M=sum(x['n']==0 for x in rows),
                distance_missing=sum(x['distance']<0 for x in rows),
                M_present_unknown_distance=sum(x['n']>0 and x['distance']<0 for x in rows))


def build_strata(rows,old):
    selected=old.get('selected_latest',[])
    if len(selected)!=6 or len({x['population'] for x in selected})!=6:
        raise ValueError('EXACT_SIX_SITES_REQUIRED')
    if len(rows)!=512 or (sum(x['v'] for x in rows),sum(x['w'] for x in rows))!=(5081,8723):
        raise ValueError('ORIGINAL_SOURCE_UNIT_CHANGED')
    selected_years={(x['year'],x['population']) for x in selected}
    raw=defaultdict(list)
    for x in rows:
        if x['n']>0 and x['distance']>=0 and (x['year'],x['population']) in selected_years:
            raw[(x['year'],x['population'])].append(x)
    old_predictions=defaultdict(list)
    for p in old['predictions']:
        if (p['year'],p['population']) in selected_years:
            old_predictions[(p['year'],p['population'])].append(p)
    strata=[]
    for x in sorted(selected,key=lambda k:k['population']):
        key=(x['year'],x['population'])
        tr=raw[key]; test=old_predictions[key]
        train_keys={(z['year'],z['population'],z['line'],z['quad']) for z in tr}
        if len(train_keys)!=len(tr):raise ValueError('SOURCE_DUPLICATE_QUADRAT')
        test_keys={(z['year'],z['population'],z['line'],z['quad']) for z in test}
        if len(test_keys)!=len(test) or not test_keys.issubset(train_keys):raise ValueError('PR185_EVALUATION_MISMATCH')
        if len(test)!=x['n_scored_quadrats'] or sum(t['n'] for t in test)!=x['n_scored_plants']:
            raise ValueError('PR185_COUNT_MISMATCH')
        if not all(all(t[k]==next(z[k] for z in tr if (z['year'],z['population'],z['line'],z['quad'])== (t['year'],t['population'],t['line'],t['quad'])) for k in ('v','w','l','s','n','distance')) for t in test):
            raise ValueError('PR185_SOURCE_COUNTS_MISMATCH')
        strata.append(dict(population=x['population'],year=x['year'],train_pool=tr,test=test,
                           all_quadrats=len(tr),heldout_lines=x['n_folds']))
    if sum(len(s['test']) for s in strata)!=168 or sum(sum(z['n'] for z in s['test']) for s in strata)!=6692:
        raise ValueError('PR185_TOTAL_SUPPORT_MISMATCH')
    return strata


def select_source_and_test(strata,heldout_population):
    chosen=[s for s in strata if s['population']==heldout_population]
    if len(chosen)!=1:raise ValueError('MISSING_HELDOUT_SITE')
    train=[x for s in strata if s['population']!=heldout_population for x in s['train_pool']]
    test=chosen[0]['test']
    if set(z['population'] for z in train).intersection({heldout_population}):raise ValueError('SITE_LEAKAGE')
    if any(t['n']<=0 or t['distance']<0 for t in test+train):raise ValueError('OUT_OF_SUPPORT')
    return train,test


def balanced_training(train):
    group=defaultdict(float)
    for x in train:group[x['population']]+=x['n']
    if len(group)!=5 or min(group.values())<=0:raise ValueError('FIVE_TRAINING_SITES_REQUIRED')
    mean_n=sum(group.values())/len(group)
    return [dict(x, n=x['n']*mean_n/group[x['population']],v=x['v']*mean_n/group[x['population']]) for x in train]


def weighted_logloss(rows,pred):
    p=np.asarray(pred,dtype=float)
    n=np.array([x['n'] for x in rows],float)
    v=np.array([x['v'] for x in rows],float)
    return float(np.sum(-v*np.log(p)-(n-v)*np.log1p(-p))/sum(n))


def site_profile(train,test):
    out={}
    for k in ('distance','other_share'):
        mn=min(x[k] for x in train);mx=max(x[k] for x in train)
        outside=[x for x in test if x[k]<mn or x[k]>mx]
        out[k]=dict(train_min=float(mn),train_max=float(mx),test_min=float(min(x[k] for x in test)),
                    test_max=float(max(x[k] for x in test)),n_extrapolated_quadrats=len(outside),
                    n_extrapolated_plants=sum(x['n'] for x in outside),
                    frac_extrapolated_plants=sum(x['n'] for x in outside)/sum(x['n'] for x in test))
    return out


def evaluate(rows,old,parent_predictive):
    import time
    parent=parent_module(parent_predictive)
    strata=build_strata(rows,old)
    results=[];predictions=[];sensitivity=[];provenance=[]
    for s in strata:
        site=s['population'];train,test=select_source_and_test(strata,site)
        byspec={}
        out={m:parent.fit_binomial(train,test,features,lam=10.0) for m,features in FIELDS.items()}
        loss={m:weighted_logloss(test,out[m]) for m in FIELDS}
        oldloss={m:weighted_logloss(test,[t['p_'+m] for t in test]) for m in FIELDS}
        br={m:float(np.mean((np.array([x['v']/x['n'] for x in test])-np.array(out[m]))**2)) for m in FIELDS}
        site_rec=dict(population=site,year=s['year'],n_train_quadrats=len(train),n_train_plants=sum(x['n'] for x in train),
                      n_train_sites=len({x['population'] for x in train}),test_n=len(test),test_plants=sum(x['n'] for x in test),
                      heldout_original_lines=s['heldout_lines'],loss=loss,within_site_reference_loss=oldloss,
                      brier_equal_quadrat=br,gain_distance=loss['BASE']-loss['DIST'],
                      gain_other_over_distance=loss['DIST']-loss['BOTH'],
                      gain_congener=loss['BASE']-loss['CONGEN'],extrapolation=site_profile(train,test))
        results.append(site_rec)
        for i,x in enumerate(test):
            rec={k:x[k] for k in ('year','population','line','quad','distance','v','w','l','s','n','other_share')}
            rec.update({'p_outsite_'+m:float(out[m][i]) for m in FIELDS})
            rec.update({'p_within_'+m:float(x['p_'+m]) for m in FIELDS})
            predictions.append(rec)
        variants=[('ridge1',1.0,train),('ridge100',100.0,train),('balanced_ridge10',10.0,balanced_training(train))]
        for label,lam,modified in variants:
            pp={m:parent.fit_binomial(modified,test,f,lam=lam) for m,f in FIELDS.items()}
            ll={m:weighted_logloss(test,pp[m]) for m in FIELDS}
            sensitivity.append(dict(population=site,scenario=label,lambda_ridge=lam,loss=ll,
                                   gain_distance=ll['BASE']-ll['DIST'],gain_other_over_distance=ll['DIST']-ll['BOTH']))
        provenance.append(dict(site=site,train_site_names=sorted(set(x['population'] for x in train)),
                               source_locs_exactly_holdout=all(x['population']!=site for x in train),
                               train_years=sorted(set(x['year'] for x in train)),test_year=s['year'],
                               n_train_rows=len(train),n_test_rows=len(test)))
    mean=lambda values:float(np.mean(values))
    aggregate=lambda items:dict(n_sites=len(items),losses={m:mean([x['loss'][m] for x in items]) for m in FIELDS},
                                gain_distance=mean([x['gain_distance'] for x in items]),
                                gain_other_over_distance=mean([x['gain_other_over_distance'] for x in items]),
                                sites_distance_better=sum(x['gain_distance']>0 for x in items),
                                sites_adding_other_better=sum(x['gain_other_over_distance']>0 for x in items))
    main=aggregate(results)
    main.update(n_quadrats_scored=len(predictions),n_plants_scored=sum(p['n'] for p in predictions),
                n_heldout_line_support=sum(s['heldout_lines'] for s in strata),n_site_names=len(strata))
    insite_loss={m:mean([x['within_site_reference_loss'][m] for x in results]) for m in FIELDS}
    sens=[]
    for scenario in ('ridge1','balanced_ridge10','ridge100'):
        xs=[x for x in sensitivity if x['scenario']==scenario]
        sens.append(dict(scenario=scenario,lambda_ridge=xs[0]['lambda_ridge'],n_sites=len(xs),
                         mean_loss={m:mean([x['loss'][m] for x in xs]) for m in FIELDS},
                         gain_distance=mean([x['gain_distance'] for x in xs]),
                         gain_other_over_distance=mean([x['gain_other_over_distance'] for x in xs]),
                         n_sites_distance_better=sum(x['gain_distance']>0 for x in xs),
                         n_sites_adding_other_better=sum(x['gain_other_over_distance']>0 for x in xs)))
    loo=[]
    for x in results:
        xx=[u for u in results if u['population']!=x['population']]
        loo.append(dict(dropped=x['population'],**aggregate(xx)))
    return dict(schema='fcp_sisyrinchium_whole_site_out_transfer_v1',
                protocol=dict(whole_site_name_withheld=True,new_external_sites=False,source_exposed_not_preregistered=True,
                              primary_score='site-equal, M plant weighted heldout line binomial log loss',
                              training='other FIVE sites only, all M-present valid-distance quadrats in locked selected years',
                              heldout='exact same 168 source quadrats as parent PR185 main',
                              ridge_lambda_primary=10.0,clip=[0.01,0.99],source_doi='10.5061/dryad.j8q34',
                              no_new_p_values=True,source_morphs_not_cross_species_mixed=True),
                original_source_qa=source_qa(rows),main=main,within_site_reference=insite_loss,
                site_scores=results,predictions=predictions,site_name_leakage_audit=provenance,
                sensitivity=sens,site_delete_sensitivity=loo,sensitivity_per_site=sensitivity)


def write_results(result,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    for name,rows in [('site_scores.csv',result['site_scores']),('predictions.csv',result['predictions']),
                       ('sensitivity_site.csv',result['sensitivity_per_site'])]:
        # Flatten only the nested human-readable score summaries; result.json is canonical.
        for row in rows:
            if any(isinstance(x,dict) for x in row.values()):
                continue
        keys=sorted(set(k for row in rows for k in row))
        with open(out/name,'w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,keys,lineterminator='\n');writer.writeheader()
            for row in rows:writer.writerow({k:json.dumps(row[k],sort_keys=True) if isinstance(row.get(k),dict) else row.get(k,'') for k in keys})
    (out/'site_name_leakage_audit.json').write_text(json.dumps(result['site_name_leakage_audit'],indent=2)+'\n')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--archive',required=True)
    p.add_argument('--parent',required=True)
    p.add_argument('--parent-code',required=True)
    p.add_argument('--outdir',default='results')
    a=p.parse_args()
    b=Path(a.parent).read_bytes()
    if hashlib.sha256(b).hexdigest()!=PARENT_RESULTS_SHA256:
        raise ValueError('PARENT_SCIENTIFIC_RESULT_HASH_CHANGED')
    old=json.loads(b)
    rows=load_original_verified(a.archive,a.parent_code)
    result=evaluate(rows,old,a.parent_code)
    write_results(result,a.outdir)
    print(json.dumps({'status':'CALCULATED','main':result['main'],'within_site_reference':result['within_site_reference'],
                      'sites':[(x['population'],x['gain_distance']) for x in result['site_scores']]},ensure_ascii=False))


if __name__=='__main__':main()
