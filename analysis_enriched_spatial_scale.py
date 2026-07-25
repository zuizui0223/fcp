#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

WITHIN = re.compile(r"\b(within[- ]populations?|same populations?|coexist|co-occurr|morph frequenc|frequency[- ]dependent|polymorphic populations?|multiple (?:colou?r|flower) morphs|colour morphs|color morphs)\b", re.I)
GEOGRAPHIC = re.compile(r"\b(geographic|spatial variation|among populations|between populations|population differentiation|cline|hybrid zone|range edge|local adaptation|translocation|regional variation)\b", re.I)


def zscore(s):
    x=pd.to_numeric(s,errors='coerce'); sd=x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0


def baseline_scale(row):
    text=' '.join(str(row.get(c,'') or '') for c in ['best_title','best_match_evidence','review_reason'])
    within=bool(WITHIN.search(text)); among=bool(GEOGRAPHIC.search(text))
    if within and not among: return 'within_population'
    if among and not within: return 'among_population'
    if within and among: return 'mixed'
    return 'unclear'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--review',required=True)
    ap.add_argument('--works',required=True)
    ap.add_argument('--geography',required=True)
    ap.add_argument('--outdir',default='analysis_outputs/evidence_scale_enriched')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    review=pd.read_csv(a.review); works=pd.read_csv(a.works); geo=pd.read_csv(a.geography)

    # Retain only direct natural-population flower-colour evidence.
    eligible=works[(works.direct_colour_signal==1)&(works.artificial_signal==0)&(works.score>=20)].copy()
    agg=eligible.groupby('canonical_name').agg(
        n_eligible=('openalex_id','nunique'),
        enriched_within=('within_signal','max'),
        enriched_among=('among_signal','max')
    ).reset_index()

    d=review.merge(agg,on='canonical_name',how='left')
    for c in ['n_eligible','enriched_within','enriched_among']:
        d[c]=d[c].fillna(0).astype(int)
    d['baseline_scale']=d.apply(baseline_scale,axis=1)

    # Enrichment supplements ambiguous baseline records; it never erases an existing
    # unambiguous within- or among-population classification.
    final=[]; sources=[]
    for _,row in d.iterrows():
        base=row['baseline_scale']; w=int(row['enriched_within']); g=int(row['enriched_among'])
        if base in {'within_population','among_population'}:
            if (base=='within_population' and g) or (base=='among_population' and w):
                final.append('mixed'); sources.append('baseline_plus_conflicting_enrichment')
            else:
                final.append(base); sources.append('baseline_unambiguous')
        else:
            if w and g: final.append('mixed')
            elif w: final.append('within_population')
            elif g: final.append('among_population')
            else: final.append(base)
            sources.append('high_confidence_enrichment' if (w or g) else 'unresolved')
    d['enriched_scale']=final
    d['classification_source']=sources
    d['manual_review_required']=d.enriched_scale.isin(['mixed','unclear'])
    d.to_csv(out/'high_confidence_enriched_classification.csv',index=False)
    d.loc[d.manual_review_required].to_csv(out/'high_confidence_manual_review.csv',index=False)

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
        range_row=tab.set_index('term').loc['range_z']
        result={'status':'complete','n_rows':int(len(x)),'class_counts':counts,'range_odds_ratio':float(range_row['odds_ratio']),'range_ci_low':float(range_row['ci_low']),'range_ci_high':float(range_row['ci_high']),'range_p_value':float(range_row['p_value'])}
    manifest={
        'n_confirmed':int(len(review)),
        'eligible_natural_direct_works':int(len(eligible)),
        'species_with_eligible_works':int(eligible.canonical_name.nunique()),
        'baseline_scale_counts':d.baseline_scale.value_counts().to_dict(),
        'final_scale_counts':d.enriched_scale.value_counts().to_dict(),
        'classification_source_counts':d.classification_source.value_counts().to_dict(),
        'manual_review_rows':int(d.manual_review_required.sum()),
        'model':result,
        'semantic_guard':'high-confidence enrichment supplements rather than replaces baseline evidence; direct colour evidence, no artificial signal, score >=20'
    }
    (out/'high_confidence_enriched_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(manifest,ensure_ascii=False))


if __name__=='__main__': main()
