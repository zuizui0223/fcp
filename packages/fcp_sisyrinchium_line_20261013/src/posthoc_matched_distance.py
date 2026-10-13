#!/usr/bin/env python3
"""BOOKKEEPING only: compare raw vs adjusted exposure on IDENTICAL line subset.

This does not create a new inferential test, and does not replace protocol estimates.
"""
import argparse, json, csv
from pathlib import Path
import numpy as np
from within_line import load_source,_stratify,distance_support

def run(source,result,out):
    d=load_source(source)
    ref=json.loads(Path(result).read_text())
    raw=[]
    for (year,pop),g in d.groupby(['year','population'],sort=True):
        gg=_stratify(g.n.to_numpy(),g.v.to_numpy(),g.a.to_numpy(),g.line.to_numpy(),distance=g.distance.to_numpy())
        support=distance_support(gg)
        if support['status']!='ELIGIBLE':continue
        info=[q for q in gg if q['distance_ok']]
        total=sum(q['D'] for q in info)
        unadjusted=float(sum(np.dot(q['v']-q['n']*q['p'],q['a']) for q in info)/total)
        adjusted=float(sum(np.dot(q['v']-q['n']*q['p'],q['a_residual']) for q in info)/total)
        stored=ref['full']['results']['distance']['siteyear_estimates'][year+'|'+pop]
        assert abs(stored-adjusted)<1e-12
        raw.append({'siteyear':year+'|'+pop,'n_informative_lines':len(info),
                    'n_target_M_in_informative_lines':support['N'],
                    'unadjusted_same_support':unadjusted,
                    'adjusted_linear_distance_same_support':adjusted,
                    'difference_due_to_linear_distance_projection':unadjusted-adjusted})
    assert len(raw)==8
    avg_raw=float(np.mean([x['unadjusted_same_support'] for x in raw]))
    avg_adj=float(np.mean([x['adjusted_linear_distance_same_support'] for x in raw]))
    out.mkdir(parents=True,exist_ok=True)
    with (out/'posthoc_matched_distance_support.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(raw[0]));w.writeheader();w.writerows(raw)
    summary={'status':'BOOKKEEPING_POSTHOC_NOT_NEW_PRIMARY_ENDPOINT',
             'n_siteyears':len(raw),'unadjusted_same_lines_mean':avg_raw,'distance_linear_adjusted_mean':avg_adj,
             'magnitude_reduction_ratio':1-avg_adj/avg_raw if avg_raw else None,
             'do_not_call_causal_distance_mediation':True,
             'cannot_separate_distance_habitat_congeners_or_spatial_autocorrelation':True,
             'nominal_adjusted_p':ref['full']['results']['distance']['p_one_sided'],
             'without_banpaku':next((x for x in ref['full']['results']['distance']['posthoc_leave_one_siteyear_out'] if x['dropped_siteyear']=='Y13|banpaku'),None),
             'repeated_sites_are_not_independent_evolutionary_replicates':True}
    (out/'posthoc_matched_distance_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,ensure_ascii=False))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default='source/frequency.txt');p.add_argument('--result',default='results/result.json');p.add_argument('--outdir',default='results');a=p.parse_args()
    run(Path(a.source),Path(a.result),Path(a.outdir))
