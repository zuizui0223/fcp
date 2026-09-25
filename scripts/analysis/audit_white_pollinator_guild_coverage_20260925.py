#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from collections import defaultdict
from pathlib import Path
import pandas as pd

TYPES={'visitsFlowersOf','flowersVisitedBy','pollinates','pollinatedBy'}

def norm(x):
    return re.sub(r'\s+',' ',str(x).strip()).lower()

def guild(path,name):
    s=(str(path)+' '+str(name)).lower()
    if 'sphingidae' in s: return 'hawkmoth'
    if 'lepidoptera' in s: return 'other_lepidoptera'
    if 'diptera' in s: return 'diptera'
    if 'hymenoptera' in s: return 'hymenoptera'
    if re.search(r'(^|[|;, ])aves($|[|;, ])',s): return 'bird'
    if 'chiroptera' in s: return 'bat'
    return 'other_or_unresolved'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True)
    p.add_argument('--interactions',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--chunksize',type=int,default=150000)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    fcp=pd.read_csv(a.measured,usecols=['species'])
    species=sorted(set(fcp.species.dropna().astype(str)))
    key_to_name={norm(s):s for s in species}
    targets=set(key_to_name)

    cols=['sourceTaxonName','sourceTaxonPathNames','interactionTypeName','targetTaxonName','targetTaxonPathNames']
    header=pd.read_csv(a.interactions,sep='\t',compression='gzip',nrows=0,encoding='utf-8',encoding_errors='replace')
    missing=sorted(set(cols)-set(header.columns))
    if missing: raise SystemExit(f'missing GloBI columns: {missing}')

    rec=defaultdict(int); visitors=defaultdict(set); gcounts=defaultdict(lambda:defaultdict(int))
    scanned=0; flower_rows=0
    for chunk in pd.read_csv(a.interactions,sep='\t',compression='gzip',usecols=cols,chunksize=a.chunksize,
                             low_memory=False,encoding='utf-8',encoding_errors='replace'):
        scanned+=len(chunk)
        q=chunk[chunk.interactionTypeName.isin(TYPES)].copy()
        flower_rows+=len(q)
        if q.empty: continue
        forward=q.interactionTypeName.isin({'visitsFlowersOf','pollinates'})
        q['plant']=q.targetTaxonName.where(forward,q.sourceTaxonName)
        q['visitor']=q.sourceTaxonName.where(forward,q.targetTaxonName)
        q['visitor_path']=q.sourceTaxonPathNames.where(forward,q.targetTaxonPathNames)
        q['plant_key']=q.plant.map(norm)
        q=q[q.plant_key.isin(targets)]
        for r in q.itertuples(index=False):
            k=r.plant_key
            rec[k]+=1
            visitors[k].add(norm(r.visitor))
            gcounts[k][guild(r.visitor_path,r.visitor)]+=1

    rows=[]
    for k,n in rec.items():
        g=gcounts[k]
        major={x for x in ['hawkmoth','other_lepidoptera','diptera','hymenoptera','bird','bat'] if g.get(x,0)>0}
        rows.append({
          'species':key_to_name[k],
          'n_flower_records':int(n),
          'n_visitor_taxa':int(len(visitors[k])),
          'n_major_guilds':int(len(major)),
          'hawkmoth':int(g.get('hawkmoth',0)),
          'other_lepidoptera':int(g.get('other_lepidoptera',0)),
          'lepidoptera_total':int(g.get('hawkmoth',0)+g.get('other_lepidoptera',0)),
          'diptera':int(g.get('diptera',0)),
          'hymenoptera':int(g.get('hymenoptera',0)),
          'bird':int(g.get('bird',0)),
          'bat':int(g.get('bat',0)),
          'other_or_unresolved':int(g.get('other_or_unresolved',0)),
        })
    d=pd.DataFrame(rows)
    if d.empty:
        d=pd.DataFrame(columns=['species','n_flower_records','n_visitor_taxa','n_major_guilds','hawkmoth','other_lepidoptera','lepidoptera_total','diptera','hymenoptera','bird','bat','other_or_unresolved'])
    d.to_csv(out/'pollinator_guild_coverage_by_species.csv',index=False)

    n1=int((d.n_flower_records>=1).sum()) if len(d) else 0
    n5=int((d.n_flower_records>=5).sum()) if len(d) else 0
    n10=int((d.n_flower_records>=10).sum()) if len(d) else 0
    n20=int((d.n_flower_records>=20).sum()) if len(d) else 0
    lepi10=int(((d.n_flower_records>=10)&(d.lepidoptera_total>0)).sum()) if len(d) else 0
    hawk10=int(((d.n_flower_records>=10)&(d.hawkmoth>0)).sum()) if len(d) else 0
    mixed10=int(((d.n_flower_records>=10)&(d.hymenoptera>0)&((d.lepidoptera_total+d.diptera+d.bird+d.bat)>0)).sum()) if len(d) else 0
    result={
      'schema':'fcp_pollinator_guild_coverage_v1',
      'status':'complete',
      'outcome_blind':True,
      'source':'GloBI Review Dataset Corpus v3 / Zenodo 18064921 / interactions.tsv.gz',
      'interaction_types':sorted(TYPES),
      'fcp_species_total':len(species),
      'integrated_rows_scanned':int(scanned),
      'flower_interaction_rows_scanned':int(flower_rows),
      'coverage':{
        'species_ge1_record':n1,
        'species_ge5_records':n5,
        'species_ge10_records':n10,
        'species_ge20_records':n20,
        'species_ge10_with_lepidoptera':lepi10,
        'species_ge10_with_hawkmoth':hawk10,
        'species_ge10_hymenoptera_plus_nonhymenoptera':mixed10
      },
      'gates':{
        'general_guild':{'minimum_species':100,'pass':bool(n10>=100)},
        'lepidoptera':{'minimum_species':50,'pass':bool(lepi10>=50)},
        'hawkmoth':{'minimum_species':30,'pass':bool(hawk10>=30)}
      },
      'hard_boundary':'No FCP flower-colour outcome was read in this coverage audit.'
    }
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
