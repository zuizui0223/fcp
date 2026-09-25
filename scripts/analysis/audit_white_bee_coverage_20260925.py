#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import pandas as pd

def norm(x):
    return re.sub(r'\s+',' ',str(x).strip()).lower()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True)
    p.add_argument('--bee-curated',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--chunksize',type=int,default=200000)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    fcp=pd.read_csv(a.measured,usecols=['species'])
    species=sorted(set(fcp.species.dropna().astype(str)))
    species_norm={norm(x):x for x in species}
    target=set(species_norm)

    header=pd.read_csv(a.bee_curated,nrows=0)
    if 'plant_species' not in header.columns:
        raise SystemExit(f"curated GloBI missing plant_species; columns={list(header.columns)}")

    matched=set()
    all_plants=set()
    scanned_rows=0
    for chunk in pd.read_csv(a.bee_curated,usecols=['plant_species'],chunksize=a.chunksize,low_memory=False):
        scanned_rows += len(chunk)
        vals=set(norm(x) for x in chunk['plant_species'].dropna().astype(str))
        all_plants.update(vals)
        matched.update(target.intersection(vals))

    matched_names=sorted(species_norm[x] for x in matched)
    coverage=len(matched_names)
    status='COVERAGE_GATE_PASS' if coverage>=100 else 'NOT_ESTIMABLE_CURRENT_BEE_DATASET'

    pd.DataFrame({'species':matched_names}).to_csv(out/'exact_matched_species.csv',index=False)
    result={
      'schema':'fcp_white_bee_coverage_v2',
      'status':status,
      'outcome_blind':True,
      'source':'Noori et al. 2026 curated GloBI v3.1 / Zenodo 18303036 / GloBI_Curated.csv',
      'plant_name_field':'plant_species',
      'fcp_species_total':len(species),
      'curated_rows_scanned':int(scanned_rows),
      'curated_unique_plant_names_normalized':int(len(all_plants)),
      'exact_matched_species':coverage,
      'coverage_fraction':coverage/len(species) if species else None,
      'minimum_required_species':100,
      'hard_boundary':'No flower-colour outcome was read. Coverage alone cannot support or reject pollinator causation.'
    }
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
