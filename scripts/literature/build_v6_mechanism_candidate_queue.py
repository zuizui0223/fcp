#!/usr/bin/env python3
"""Build a conservative *candidate* mechanism queue for the v6 FCP core.

Keyword flags are navigation aids only. They are never promoted to ecological evidence
without manual abstract/full-text verification. This script deliberately keeps source
classification (C/S) separate from mechanism coding to expose same-literature ascertainment.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd

PATTERNS = {
    "pollinator_candidate": r"pollinat|bee[s]?\b|butterfl|moth|syrph|fly\b|flies\b|bird[s]?\b|beetle",
    "temporal_candidate": r"season|year[- ]to[- ]year|interannual|temporal|phenolog|flowering time|earlier|later",
    "abiotic_candidate": r"soil|serpentine|drought|precipitation|temperature|climat|ultraviolet|\bUV\b|elevation|altitud|stress",
    "herbivory_candidate": r"herbiv|floriv|seed predat|nectar rob",
    "mating_candidate": r"self[- ]compat|self[- ]incompat|selfing|autogam|outcross|mating system|reproductive assurance",
    "genetic_candidate": r"gene flow|genetic structure|population genetic|allozyme|drift|founder|migration|inheritance|locus|allele",
    "local_adaptation_candidate": r"reciprocal transplant|translocation|local adaptation|locally adapt|divergent selection|ecotypic differentiation|geographic mosaic|geographical mosaic",
    "frequency_dependence_candidate": r"frequency[- ]depend|rare morph|negative frequency|positive frequency",
    "antagonism_candidate": r"opposing selection|counter(?:act|ed|s|ing)?|antagoni|trade[- ]?off|conflicting selection",
}

def hit(text: str, pat: str) -> int:
    return int(bool(re.search(pat, text, flags=re.I)))

def snippet(text: str, pats: list[str], width: int=360) -> str:
    for pat in pats:
        m=re.search(pat,text,flags=re.I)
        if m:
            lo=max(0,m.start()-100); hi=min(len(text),m.end()+width-100)
            return re.sub(r"\s+"," ",text[lo:hi]).strip()
    return ""

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--core-sources',required=True)
    p.add_argument('--record-screen',required=True)
    p.add_argument('--core-species',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    c=pd.read_csv(a.core_sources)
    r=pd.read_csv(a.record_screen,usecols=['record_review_id','abstract','journal','url'])
    s=pd.read_csv(a.core_species,usecols=['canonical_name','organization_state'])
    d=c.merge(r,on='record_review_id',how='left',validate='many_to_one')
    d=d.merge(s,left_on='accepted_name',right_on='canonical_name',how='left',validate='many_to_one')
    d['abstract']=d['abstract'].fillna('').astype(str)
    text=(d['title'].fillna('')+' '+d['abstract']).astype(str)
    for k,pat in PATTERNS.items(): d[k]=[hit(x,pat) for x in text]
    d['candidate_mechanism_count']=d[list(PATTERNS)].sum(axis=1)
    d['navigation_excerpt']=[snippet(x,list(PATTERNS.values())) for x in text]
    keep=['record_review_id','source_id','accepted_name','family','organization_state','title','year','journal',
          *PATTERNS.keys(),'candidate_mechanism_count','navigation_excerpt','url']
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    d[keep].sort_values(['accepted_name','source_id']).to_csv(out,index=False)
    print({'rows':len(d),'species':d.accepted_name.nunique(),'abstract_present':int((d.abstract.str.len()>0).sum())})
if __name__=='__main__': main()
