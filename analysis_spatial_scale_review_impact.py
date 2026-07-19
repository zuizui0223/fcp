#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def zscore(s):
    x=pd.to_numeric(s,errors='coerce'); sd=x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0


def fit_model(d):
    x=d.dropna(subset=['latitudinal_range','absolute_latitude','gbif_occurrences','family_review']).copy()
    x['is_within']=(x.enriched_scale=='within_population').astype(int)
    x['range_z']=zscore(np.log1p(pd.to_numeric(x.latitudinal_range,errors='coerce').clip(lower=0)))
    x['latitude_z']=zscore(pd.to_numeric(x.absolute_latitude,errors='coerce'))
    x['gbif_z']=zscore(np.log1p(pd.to_numeric(x.gbif_occurrences,errors='coerce').clip(lower=0)))
    m=smf.glm('is_within ~ range_z + latitude_z + gbif_z',x,family=sm.families.Binomial()).fit()
    return float(m.params['range_z']), float(np.exp(m.params['range_z'])), int(len(x))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--classification',required=True)
    ap.add_argument('--geography',required=True)
    ap.add_argument('--outdir',default='analysis_outputs/review_impact')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.classification)
    g=pd.read_csv(a.geography)
    x=d.merge(g,on='canonical_name',how='left',suffixes=('_review','_geo'))
    fam='family_review' if 'family_review' in x.columns else 'family'
    x=x.rename(columns={fam:'family_review'})
    base=x[x.enriched_scale.isin(['within_population','among_population'])].copy()
    base_beta,base_or,base_n=fit_model(base)
    amb=x[x.enriched_scale.isin(['mixed','unclear'])].copy()
    rows=[]
    for _,r in amb.iterrows():
        if pd.isna(r.get('latitudinal_range')) or pd.isna(r.get('absolute_latitude')) or pd.isna(r.get('gbif_occurrences')):
            continue
        one=r.to_frame().T.copy()
        one['enriched_scale']='within_population'; bw,orw,_=fit_model(pd.concat([base,one],ignore_index=True))
        one['enriched_scale']='among_population'; ba,ora,_=fit_model(pd.concat([base,one],ignore_index=True))
        rows.append({
            'canonical_name':r.get('canonical_name',''),
            'family':r.get('family_review',''),
            'current_scale':r.get('enriched_scale',''),
            'baseline_scale':r.get('baseline_scale',''),
            'classification_source':r.get('classification_source',''),
            'n_eligible':r.get('n_eligible',0),
            'latitudinal_range':r.get('latitudinal_range'),
            'absolute_latitude':r.get('absolute_latitude'),
            'gbif_occurrences':r.get('gbif_occurrences'),
            'range_beta_if_within':bw,
            'range_beta_if_among':ba,
            'range_or_if_within':orw,
            'range_or_if_among':ora,
            'classification_impact':abs(bw-ba),
            'max_abs_change_from_current':max(abs(bw-base_beta),abs(ba-base_beta)),
            'recommended_review_priority':'P0' if abs(bw-ba)>=0.08 else 'P1' if abs(bw-ba)>=0.04 else 'P2'
        })
    tab=pd.DataFrame(rows).sort_values(['classification_impact','max_abs_change_from_current'],ascending=False)
    tab.to_csv(out/'spatial_scale_review_impact_priority.csv',index=False)
    top=tab.head(15).copy(); top['manual_decision']=''; top['decision_evidence']=''; top['reviewer_notes']=''
    top.to_csv(out/'top15_decisive_manual_review.csv',index=False)
    manifest={
        'classified_model_n':base_n,
        'ambiguous_with_geography':int(len(tab)),
        'current_range_beta':base_beta,
        'current_range_odds_ratio':base_or,
        'top15_max_impact':float(top.classification_impact.max()) if len(top) else None,
        'top15_min_impact':float(top.classification_impact.min()) if len(top) else None,
        'priority_counts':tab.recommended_review_priority.value_counts().to_dict() if len(tab) else {},
        'interpretation':'Review impact is the difference in the range-size coefficient when the same unresolved species is assigned within versus among.',
        'semantic_guard':'This prioritizes manual review; it does not assign biological classes automatically.'
    }
    (out/'spatial_scale_review_impact_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(manifest,ensure_ascii=False))

if __name__=='__main__': main()
