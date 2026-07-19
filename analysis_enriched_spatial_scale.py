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


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--review',required=True)
    ap.add_argument('--works',required=True)
    ap.add_argument('--geography',required=True)
    ap.add_argument('--outdir',default='analysis_outputs/evidence_scale_enriched')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    review=pd.read_csv(a.review); works=pd.read_csv(a.works); geo=pd.read_csv(a.geography)
    eligible=works[(works.direct_colour_signal==1)&(works.artificial_signal==0)&(works.score>=20)].copy()
    agg=eligible.groupby('canonical_name').agg(n_eligible=('openalex_id','nunique'),within=('within_signal','max'),among=('among_signal','max')).reset_index()
    d=review.merge(agg,on='canonical_name',how='left')
    for c in ['n_eligible','within','among']: d[c]=d[c].fillna(0).astype(int)
    d['enriched_scale']=np.select([
        (d.within.eq(1)&d.among.eq(0)),
        (d.within.eq(0)&d.among.eq(1)),
        (d.within.eq(1)&d.among.eq(1))],
        ['within_population','among_population','mixed'],default='unclear')
    d.to_csv(out/'high_confidence_enriched_classification.csv',index=False)
    x=d.merge(geo,on='canonical_name',how='left',suffixes=('_review','_geo'))
    x=x[x.enriched_scale.isin(['within_population','among_population'])].dropna(subset=['latitudinal_range','absolute_latitude','gbif_occurrences']).copy()
    result={'status':'not_run','reason':'fewer than 15 rows in either class'}
    counts=x.enriched_scale.value_counts().to_dict()
    if counts.get('within_population',0)>=15 and counts.get('among_population',0)>=15:
        x['is_within']=(x.enriched_scale=='within_population').astype(int)
        x['range_z']=zscore(np.log1p(x.latitudinal_range.clip(lower=0)))
        x['latitude_z']=zscore(x.absolute_latitude)
        x['gbif_z']=zscore(np.log1p(x.gbif_occurrences.clip(lower=0)))
        fam='family_review' if 'family_review' in x.columns else 'family'
        m=smf.glm('is_within ~ range_z + latitude_z + gbif_z',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x[fam]})
        ci=m.conf_int(); tab=pd.DataFrame({'term':m.params.index,'estimate':m.params.values,'std_error':m.bse.values,'odds_ratio':np.exp(m.params.values),'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values),'p_value':m.pvalues.values})
        tab.to_csv(out/'high_confidence_within_vs_among_model.csv',index=False)
        result={'status':'complete','n_rows':int(len(x)),'class_counts':counts,'range_odds_ratio':float(tab.set_index('term').loc['range_z','odds_ratio']),'range_p_value':float(tab.set_index('term').loc['range_z','p_value'])}
    manifest={'n_confirmed':int(len(review)),'eligible_natural_direct_works':int(len(eligible)),'species_with_eligible_works':int(eligible.canonical_name.nunique()),'scale_counts':d.enriched_scale.value_counts().to_dict(),'model':result,'semantic_guard':'classification requires direct flower-colour evidence, no artificial signal, and score >=20'}
    (out/'high_confidence_enriched_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(manifest,ensure_ascii=False))

if __name__=='__main__': main()
