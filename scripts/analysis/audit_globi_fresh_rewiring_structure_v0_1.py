#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path

def read_projected(path,cols):
    with path.open(newline="",encoding="utf-8-sig",errors="replace") as f:
        r=csv.reader(f); h=next(r); idx={c:h.index(c) for c in cols if c in h}
        miss=sorted(set(cols)-set(idx))
        if miss: raise RuntimeError(f"missing columns {miss}")
        for row in r: yield {c:(row[i] if i<len(row) else "") for c,i in idx.items()}

def key(row,cols): return "|".join(str(row.get(c,"")).strip() for c in cols)

def summarize(rows,focal,keycols):
    by=defaultdict(set)
    for row in rows:
        sp=row.get(focal,"").strip(); k=key(row,keycols)
        if sp and k: by[k].add(sp)
    multi={k:v for k,v in by.items() if len(v)>=2}
    taxa=set().union(*multi.values()) if multi else set()
    pairs=Counter()
    for ss in multi.values():
        for p in combinations(sorted(ss),2): pairs[p]+=1
    repeated=sum(n>=2 for n in pairs.values())
    return {
      "local_units":len(by),
      "local_units_with_at_least_2_focal_taxa":len(multi),
      "focal_taxa_in_multifocal_units":len(taxa),
      "focal_pairs_repeated_in_at_least_2_units":repeated,
      "maximum_focal_taxa_per_local_unit":max((len(x) for x in by.values()),default=0)
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--carni",type=Path,required=True)
    ap.add_argument("--croppol-field",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(); d=json.loads(a.design.read_text()); g=d["gate"]
    carni_cols=["scientificNameCarni","sourcePrimaryReference","verbatimLocality","decimalLatitude","decimalLongitude","startYear","endYear","samplingProtocol","methodQuantification"]
    crop_cols=["study_id","site_id","crop"]
    cr=summarize(read_projected(a.carni,carni_cols),"scientificNameCarni",carni_cols[1:])
    cp=summarize(read_projected(a.croppol_field,crop_cols),"crop",["study_id","site_id"])
    def decide(x):
        tests={
          "multifocal_units":x["local_units_with_at_least_2_focal_taxa"]>=g["minimum_local_units_with_at_least_2_focal_taxa"],
          "multifocal_taxa":x["focal_taxa_in_multifocal_units"]>=g["minimum_focal_taxa_participating_in_multifocal_units"],
          "repeated_focal_pairs":x["focal_pairs_repeated_in_at_least_2_units"]>=g["minimum_focal_pairs_repeated_in_at_least_2_units"]
        }
        x["gate_tests"]=tests;x["necessary_rewiring_structure_pass"]=all(tests.values())
    decide(cr);decide(cp)
    n=sum(x["necessary_rewiring_structure_pass"] for x in (cr,cp))
    status="PASS_GLOBI_FRESH_REWIRING_FOCAL_STRUCTURE" if n else "HOLD_GLOBI_FRESH_REWIRING_STRUCTURALLY_UNIDENTIFIABLE"
    out={
      "version":"v0.1","status":status,"design":str(a.design),
      "candidates":{"CARNIDIET":cr,"CROPPOL":cp},
      "passing_candidate_count":n,
      "partner_columns_used":False,"interaction_weights_used":False,"turnover_outcomes_computed":False,
      "interpretation":d["pass_meaning"] if n else d["hold_meaning"],
      "current_fcp_manuscript_changed":False,"chun_el_v0_3_changed":False
    }
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__": main()
