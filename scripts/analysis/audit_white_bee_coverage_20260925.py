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
    p.add_argument('--bee-summary',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    fcp=pd.read_csv(a.measured,usecols=['species'])
    species=sorted(set(fcp.species.dropna().astype(str)))
    species_norm={norm(x):x for x in species}

    bee=pd.read_csv(a.bee_summary,low_memory=False)
    candidates=[]
    for c in bee.columns:
        s=bee[c]
        if s.dtype!='object': continue
        vals=set(norm(x) for x in s.dropna().astype(str))
        overlap=sorted(set(species_norm).intersection(vals))
        candidates.append({
            'column':c,
            'n_nonmissing':int(s.notna().sum()),
            'n_unique_normalized':int(len(vals)),
            'exact_overlap':int(len(overlap)),
        })
    candidates=sorted(candidates,key=lambda r:(-r['exact_overlap'],r['column']))
    if not candidates: raise SystemExit('no string columns in bee summary')
    best=candidates[0]
    selected=best['column']
    vals=set(norm(x) for x in bee[selected].dropna().astype(str))
    matched=sorted(species_norm[x] for x in set(species_norm).intersection(vals))
    coverage=len(matched)
    status='COVERAGE_GATE_PASS' if coverage>=100 else 'NOT_ESTIMABLE_CURRENT_BEE_DATASET'

    pd.DataFrame({'species':matched}).to_csv(out/'exact_matched_species.csv',index=False)
    pd.DataFrame(candidates).to_csv(out/'candidate_name_columns.csv',index=False)
    result={
      'schema':'fcp_white_bee_coverage_v1',
      'status':status,
      'outcome_blind':True,
      'fcp_species_total':len(species),
      'selected_summary_name_column':selected,
      'exact_matched_species':coverage,
      'coverage_fraction':coverage/len(species) if species else None,
      'minimum_required_species':100,
      'candidate_columns_top10':candidates[:10],
      'hard_boundary':'No flower-colour outcome was read. Coverage alone cannot support or reject pollinator causation.'
    }
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
