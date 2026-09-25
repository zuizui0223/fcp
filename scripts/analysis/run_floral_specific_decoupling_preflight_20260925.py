#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd

CLASSES=("DEC","SYS")
MIN_INSTANCES=5
MIN_FAMILIES=3

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--registry",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.registry)
    required={"instance_id","species","family","architecture_class","strict_for_gate"}
    missing=sorted(required-set(d.columns))
    if missing: raise SystemExit(f"missing columns: {missing}")
    strict=d.loc[d.strict_for_gate.astype(str).str.lower().isin({"yes","true","1"})].copy()
    rows=[]
    passed=True
    for cls in CLASSES:
        q=strict.loc[strict.architecture_class.eq(cls)].copy()
        n=int(q.instance_id.nunique())
        nf=int(q.family.nunique())
        ok=(n>=MIN_INSTANCES and nf>=MIN_FAMILIES)
        passed=passed and ok
        rows.append({
            "architecture_class":cls,
            "strict_instances":n,
            "families":nf,
            "minimum_instances":MIN_INSTANCES,
            "minimum_families":MIN_FAMILIES,
            "class_gate_pass":bool(ok),
            "instances":";".join(q.instance_id.astype(str)),
            "species":";".join(q.species.astype(str)),
        })
    res=pd.DataFrame(rows)
    res.to_csv(out/"coverage_by_architecture.csv",index=False)
    unresolved=d.loc[d.architecture_class.eq("UNRESOLVED")].copy()
    unresolved.to_csv(out/"unresolved_instances.csv",index=False)
    result={
      "schema":"fcp_floral_specific_decoupling_preflight_v1",
      "status":"PLEIOTROPIC_COST_COMPARISON_AUTHORIZED" if passed else "PLEIOTROPIC_COST_COMPARISON_NOT_AUTHORIZED",
      "strict_registry_instances":int(len(strict)),
      "strict_registry_species":int(strict.species.nunique()),
      "strict_registry_families":int(strict.family.nunique()),
      "classes":rows,
      "unresolved_instances":int(len(unresolved)),
      "decision":"Proceed to independent pleiotropic-cost outcome coding only if both DEC and SYS gates pass.",
      "next_if_fail":"Do not lower thresholds. Expand only the molecular architecture registry with additional directly resolved natural white/pigmented systems; fitness/stress outcomes remain unopened for new candidates.",
      "hard_nonclaims":[
        "coverage failure is not evidence against floral-specific decoupling",
        "regulatory architecture is not equivalent to absence of pleiotropy",
        "structural mutation is not equivalent to systemic effect without extra-floral evidence"
      ]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
