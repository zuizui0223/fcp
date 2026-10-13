#!/usr/bin/env python3
"""POSTHOC diagnostic: same site/line support, alternative fixed ridge strengths."""
import argparse
import sys
import json
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_independent import read_original, independent_fit


def run(root):
    root=Path(root)
    rows=read_original(root/'source'/'Sisyrinchium_data.zip')
    base=json.loads((root/'results'/'result.json').read_text())
    summaries=[]
    for lam in (1.0,10.0,100.0):
        site_deltas=[]
        for site in base['selected_latest']:
            year,pop=site['year'],site['population']
            source=[r for r in rows if (r['year'],r['population'])==(year,pop) and r['n']>0 and r['distance']>=0]
            support=next(s['heldout_lines'] for s in base['eligible_siteyears'] if s['siteyear']==site['siteyear'])
            ll={'BASE':0.,'DIST':0.,'BOTH':0.,'CONGEN':0.};counts=0
            for line in support:
                train=[r for r in source if r['line']!=line]
                test=[r for r in source if r['line']==line]
                for model,features in [('BASE',[]),('DIST',['distance']),('BOTH',['distance','other_share']),('CONGEN',['other_share'])]:
                    predictions=independent_fit(train,test,features,lam=lam)
                    ll[model]+=sum(-r['v']*np.log(float(p))-r['w']*np.log1p(-float(p)) for r,p in zip(test,predictions))
                counts+=sum(r['n'] for r in test)
            site_deltas.append(dict(siteyear=site['siteyear'],lambda_ridge=lam,
                distance_gain_over_base=float((ll['BASE']-ll['DIST'])/counts),
                dist_minus_both_nats_per_plant=float((ll['DIST']-ll['BOTH'])/counts),
                congen_minus_both_nats_per_plant=float((ll['CONGEN']-ll['BOTH'])/counts)))
        summaries.append(dict(lambda_ridge=lam,n_sites=len(site_deltas),siteyear_deltas=site_deltas,
               mean_extra_congener_beyond_distance=float(np.mean([s['dist_minus_both_nats_per_plant'] for s in site_deltas])),
               mean_distance_gain_over_base=float(np.mean([s['distance_gain_over_base'] for s in site_deltas])),
               n_sites_distance_gain_positive=int(sum(s['distance_gain_over_base']>0 for s in site_deltas)),
               n_sites_extra_congener_positive=int(sum(s['dist_minus_both_nats_per_plant']>0 for s in site_deltas))))
    result={'status':'PASS','type':'POSTHOC_FIXED_SUPPORT_RIDGE_SENSITIVITY_NOT_NEW_PRIMARY_MODEL','same_original_six_sites':True,
       'no_p_values':True,'lambda_values_fixed_for_diagnostic':[1,10,100], 'results':summaries}
    (root/'results'/'posthoc_ridge_sensitivity.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps([{'lam':s['lambda_ridge'],'mean_extra':s['mean_extra_congener_beyond_distance'],'sites_positive':s['n_sites_extra_congener_positive'],
                       'mean_dist_gain':s['mean_distance_gain_over_base'],'sites_dist_positive':s['n_sites_distance_gain_positive']} for s in summaries]))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='.');a=p.parse_args();run(a.root)
