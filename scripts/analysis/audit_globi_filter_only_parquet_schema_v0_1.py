#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import duckdb

ROLE_PATTERNS={
  "dataset_identity":[r"dataset",r"namespace",r"source.*dataset"],
  "source_citation":[r"citation",r"reference",r"doi",r"source.*url",r"source.*id"],
  "focal_taxon":[r"^source.*taxon",r"^source.*name",r"^source.*id"],
  "partner_taxon":[r"^target.*taxon",r"^target.*name",r"^target.*id"],
  "interaction_type":[r"interaction.*type",r"interaction.*name"],
}

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--url",required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()

    con=duckdb.connect(database=":memory:")
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("INSTALL parquet")
    con.execute("LOAD parquet")

    # parquet_schema reads the Parquet metadata/footer. No biological row query is executed.
    rows=con.execute("SELECT name, type, repetition_type, logical_type FROM parquet_schema(?)",[a.url]).fetchall()
    cols=[]
    for name,typ,rep,logical in rows:
        if name in (None,"schema"):
            continue
        cols.append({"name":str(name),"type":str(typ),"repetition_type":str(rep),"logical_type":str(logical)})

    names=[x["name"] for x in cols]
    matches={}
    for role,patterns in ROLE_PATTERNS.items():
        hit=[]
        for name in names:
            low=name.lower()
            if any(re.search(p,low) for p in patterns):
                hit.append(name)
        matches[role]=sorted(set(hit))

    required=["dataset_identity","source_citation","focal_taxon","partner_taxon","interaction_type"]
    passed=all(matches[k] for k in required)
    out={
      "version":"v0.1",
      "status":"PASS_GLOBI_PARQUET_SCHEMA_REQUIRED_ROLES_PRESENT" if passed else "HOLD_GLOBI_PARQUET_SCHEMA_REQUIRED_ROLE_MISSING",
      "source_url":a.url,
      "column_count":len(cols),
      "columns":cols,
      "role_candidates":matches,
      "required_roles":required,
      "all_required_roles_have_candidates":passed,
      "rows_opened":False,
      "taxon_values_opened":False,
      "interaction_values_opened":False,
      "specialization_values_computed":False,
      "next_gate":"Freeze exact columns for filter-only support scan before reading rows." if passed else "STOP_NO_ROW_OPENING"
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
