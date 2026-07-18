#!/usr/bin/env python3
"""Reproducible baseline analyses for the flower-colour evidence dataset.

This module deliberately separates descriptive evidence diagnostics from later
biological inference. It uses only Python's standard library so CI can run it
without an environment lockfile.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path


def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def integer(v):
    try:return int(str(v or 0))
    except ValueError:return 0
def write_csv(path,fieldnames,data):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(data)
def wilson(k,n,z=1.959963984540054):
    if n<=0:return (0.0,0.0)
    p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return max(0,c-h),min(1,c+h)

def classify(r,title_min,context_min,score_min):
    title=integer(r.get('n_title_matches'));context=integer(r.get('n_context_matches'));works=integer(r.get('n_works'));score=integer(r.get('max_score'))
    return title>=title_min and context>=context_min and works>=1 and score>=score_min

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--species',default='data/global_flower_colour_species_ranked.csv')
    ap.add_argument('--outdir',default='analysis_outputs')
    a=ap.parse_args();data=rows(a.species)
    if not data:raise RuntimeError('No species rows')
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)

    # Evidence tiers are operational, not biological truth labels.
    def tier(r):
        t=integer(r.get('n_title_matches'));c=integer(r.get('n_context_matches'));w=integer(r.get('n_works'));s=integer(r.get('max_score'))
        if t>=2 and c>=2 and w>=2 and s>=24:return 'A_strong'
        if t>=1 and c>=1 and s>=18:return 'B_moderate'
        if t>=1 or c>=2:return 'C_weak'
        return 'D_minimal'
    tiers=Counter(tier(r) for r in data)
    families=defaultdict(list)
    for r in data:families[r.get('family') or 'Unknown'].append(r)
    family_rows=[]
    for fam,rs in families.items():
        strong=sum(tier(r) in {'A_strong','B_moderate'} for r in rs);lo,hi=wilson(strong,len(rs))
        family_rows.append({'family':fam,'candidate_species':len(rs),'strong_or_moderate':strong,'evidence_fraction':round(strong/len(rs),6),'wilson_low':round(lo,6),'wilson_high':round(hi,6),'median_n_works':sorted(integer(r.get('n_works')) for r in rs)[len(rs)//2]})
    family_rows.sort(key=lambda r:(-r['candidate_species'],r['family']))
    write_csv(out/'family_evidence_summary.csv',list(family_rows[0]),family_rows)

    sensitivity=[]
    for tm in (0,1,2):
        for cm in (1,2,3):
            for sm in (14,18,24):
                n=sum(classify(r,tm,cm,sm) for r in data)
                sensitivity.append({'title_min':tm,'context_min':cm,'score_min':sm,'retained_species':n,'retained_fraction':round(n/len(data),6)})
    write_csv(out/'threshold_sensitivity.csv',list(sensitivity[0]),sensitivity)

    concentration=[]
    ordered=sorted(family_rows,key=lambda r:-r['candidate_species']);total=len(data);cum=0
    for i,r in enumerate(ordered,1):
        cum+=r['candidate_species'];concentration.append({'rank':i,'family':r['family'],'candidate_species':r['candidate_species'],'cumulative_species':cum,'cumulative_fraction':round(cum/total,6)})
    write_csv(out/'family_concentration.csv',list(concentration[0]),concentration)

    qc={'input_species':len(data),'unique_families':len(families),'evidence_tiers':dict(sorted(tiers.items())),'species_with_title_support':sum(integer(r.get('n_title_matches'))>0 for r in data),'species_with_context_support':sum(integer(r.get('n_context_matches'))>0 for r in data),'top_10_family_share':round(sum(r['candidate_species'] for r in ordered[:10])/len(data),6),'method':'descriptive_evidence_baseline_v1','inference_guardrail':'Evidence tiers measure literature support, not prevalence of polymorphism among all angiosperms.'}
    (out/'baseline_qc.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(qc,ensure_ascii=False))
if __name__=='__main__':main()
