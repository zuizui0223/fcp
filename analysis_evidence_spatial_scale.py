#!/usr/bin/env python3
"""Classify whether confirmed flower-colour evidence concerns within-population coexistence or geographic differentiation.

Automated labels are screening labels, not final biological annotations. Ambiguous rows are
retained for manual review and never forced into either mechanism class.
"""
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

WITHIN = re.compile(r"\b(within[- ]population|same population|coexist|co-occurr|morph frequenc|frequency[- ]dependent|polymorphic population|multiple (?:colou?r|flower) morphs|colour morphs|color morphs)\b", re.I)
GEOGRAPHIC = re.compile(r"\b(geographic|spatial variation|among populations|between populations|population differentiation|cline|hybrid zone|range edge|local adaptation|translocation|regional variation)\b", re.I)
MAINTENANCE = re.compile(r"\b(maintenance|balancing selection|frequency[- ]dependent|pollinator[- ]mediated|assortative mating|negative frequency)\b", re.I)


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x-x.mean())/sd if sd and np.isfinite(sd) else x*0


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--review', required=True)
    ap.add_argument('--geography', required=True)
    ap.add_argument('--outdir', default='analysis_outputs/evidence_scale')
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    d=pd.read_csv(a.review)
    text=(d.get('best_title','').fillna('')+' '+d.get('best_match_evidence','').fillna('')+' '+d.get('review_reason','').fillna('')).astype(str)
    d['within_population_signal']=text.str.contains(WITHIN).astype(int)
    d['geographic_differentiation_signal']=text.str.contains(GEOGRAPHIC).astype(int)
    d['maintenance_mechanism_signal']=text.str.contains(MAINTENANCE).astype(int)
    d['evidence_spatial_scale']=np.select(
        [
            (d.within_population_signal.eq(1)&d.geographic_differentiation_signal.eq(0)),
            (d.within_population_signal.eq(0)&d.geographic_differentiation_signal.eq(1)),
            (d.within_population_signal.eq(1)&d.geographic_differentiation_signal.eq(1)),
        ],
        ['within_population','among_population','mixed'],
        default='unclear'
    )
    d['manual_review_required']=d.evidence_spatial_scale.isin(['mixed','unclear'])
    keep=['review_priority','canonical_name','family','evidence_spatial_scale','within_population_signal','geographic_differentiation_signal','maintenance_mechanism_signal','manual_review_required','best_title','best_doi','best_openalex_id','best_match_evidence']
    d[keep].to_csv(out/'evidence_spatial_scale_screening.csv', index=False)
    d.loc[d.manual_review_required, keep].to_csv(out/'evidence_spatial_scale_manual_review.csv', index=False)

    geo=pd.read_csv(a.geography)
    x=d.merge(geo,on='canonical_name',how='left',suffixes=('_review','_geo'))
    family_col='family_review' if 'family_review' in x.columns else 'family'
    x=x[x.evidence_spatial_scale.isin(['within_population','among_population'])].dropna(subset=['latitudinal_range','absolute_latitude','gbif_occurrences']).copy()
    result={'status':'not_run','reason':'fewer than 15 unambiguous rows per class'}
    counts=x.evidence_spatial_scale.value_counts().to_dict()
    if counts.get('within_population',0)>=15 and counts.get('among_population',0)>=15:
        x['is_within']=(x.evidence_spatial_scale=='within_population').astype(int)
        x['range_z']=zscore(np.log1p(pd.to_numeric(x.latitudinal_range,errors='coerce').clip(lower=0)))
        x['latitude_z']=zscore(pd.to_numeric(x.absolute_latitude,errors='coerce'))
        x['gbif_z']=zscore(np.log1p(pd.to_numeric(x.gbif_occurrences,errors='coerce').clip(lower=0)))
        model=smf.glm('is_within ~ range_z + latitude_z + gbif_z',x,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':x[family_col]})
        ci=model.conf_int()
        tab=pd.DataFrame({'term':model.params.index,'estimate':model.params.values,'std_error':model.bse.values,'odds_ratio':np.exp(model.params.values),'ci_low':np.exp(ci[0].values),'ci_high':np.exp(ci[1].values),'p_value':model.pvalues.values})
        tab.to_csv(out/'within_vs_among_population_model.csv',index=False)
        result={'status':'complete','n_rows':int(len(x)),'class_counts':counts,'formula':'is_within ~ range_z + latitude_z + gbif_z','cluster_column':family_col}
    summary={
        'n_confirmed':int(len(d)),
        'scale_counts':d.evidence_spatial_scale.value_counts().to_dict(),
        'manual_review_rows':int(d.manual_review_required.sum()),
        'semantic_guard':'automated spatial-scale labels are screening labels; ambiguous evidence is not forced',
        'ecological_question':'Does the range-size signal distinguish within-population polymorphism from among-population colour differentiation?',
        'model':result,
    }
    (out/'evidence_spatial_scale_manifest.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False))

if __name__=='__main__': main()
