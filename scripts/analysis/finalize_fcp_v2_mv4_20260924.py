#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
MODULE=ROOT/"scripts/analysis/analyze_fcp_v2_mv4_20260924.py"
SPEC=importlib.util.spec_from_file_location("fcp_v2_mv4_core_finalize",MODULE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load MV4 core")
CORE=importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name]=CORE
SPEC.loader.exec_module(CORE)

SHARDS=20


def load_condition(shards:Path, filename:str, key:str)->tuple[pd.DataFrame,np.ndarray]:
    tables=[]; nulls=[]
    for i in range(SHARDS):
        root=shards/f"shard-{i}"
        tab=pd.read_csv(root/filename)
        with np.load(root/"nulls.npz",allow_pickle=False) as z:
            arr=z[key].astype(float)
        if arr.shape!=(len(tab),CORE.N_PERM):
            raise RuntimeError(f"{filename} null shape drift in shard {i}: {arr.shape} vs {len(tab)}")
        tables.append(tab); nulls.append(arr)
    table=pd.concat(tables,ignore_index=True) if tables else pd.DataFrame()
    null=np.concatenate(nulls,axis=0) if nulls else np.empty((0,CORE.N_PERM),float)
    if len(table)!=len(null):
        raise RuntimeError(f"{filename} final table/null mismatch")
    return table,null


def panel_descriptive(table:pd.DataFrame)->dict:
    out={}
    for label in ("P","N"):
        sub=table.loc[table["panel"].astype(str).eq(label)].reset_index(drop=True)
        out[label]=CORE.panel_descriptive(sub) if len(sub) else {"species":0,"rho_D_spatial":None}
    return out


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--shards-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()

    all_selected=[]
    receipts=[]
    for i in range(SHARDS):
        root=args.shards_dir/f"shard-{i}"
        receipt=json.loads((root/"receipt.json").read_text())
        if receipt.get("status")!="MV4_SHARD_COMPLETE" or int(receipt.get("shard_index",-1))!=i:
            raise RuntimeError(f"MV4 shard receipt invalid: {i}")
        if int(receipt.get("shards",-1))!=SHARDS or int(receipt.get("per_species_null_replicates",-1))!=CORE.N_PERM:
            raise RuntimeError(f"MV4 shard design drift: {i}")
        if receipt.get("matched_background_exact_control_evaluable") is not False:
            raise RuntimeError("matched-background channel unexpectedly changed")
        receipts.append(receipt)
        all_selected.extend(receipt["selected_species"])
    if len(all_selected)!=400 or len(set(all_selected))!=400:
        raise RuntimeError("MV4 shard species census is not 400 unique species")

    base,base_null=load_condition(args.shards_dir,"base_species.csv","base")
    stable,stable_null=load_condition(args.shards_dir,"roi_stable_species.csv","roi_stable")

    base_summary=CORE.summarize(base,base_null) if len(base) else {"evaluable":False,"species":0}
    stable_summary=CORE.summarize(stable,stable_null) if len(stable) else {"evaluable":False,"species":0}

    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    base.to_csv(out/"base_species_spatial.csv",index=False,lineterminator="\n")
    stable.to_csv(out/"roi_stable_species_spatial.csv",index=False,lineterminator="\n")
    np.savez_compressed(out/"spatial_nulls.npz",base=base_null,roi_stable=stable_null)

    result={
        "schema":"fcp_v2_mv4_result_v1",
        "status":"FCP_V2_MV4_COMPLETE_WITH_FROZEN_CHANNEL_LIMITATION",
        "raw_species":400,
        "raw_rows":40000,
        "per_species_spatial_null_replicates":CORE.N_PERM,
        "base":{
            "eligible_species":int(len(base)),
            "pooled":base_summary,
            "panel_descriptive":panel_descriptive(base),
        },
        "roi_stable_only":{
            "eligible_species":int(len(stable)),
            "pooled":stable_summary,
            "panel_descriptive":panel_descriptive(stable),
        },
        "matched_flower_minus_background":{
            "evaluable":False,
            "status":"NOT_EVALUABLE_FROZEN_TECHNICAL_CHANNEL_ABSENT",
            "reason":"The pre-biological v2 technical seal did not persist background palette counts/fractions required for the existing matched flower-minus-background JSD-difference statistic. Available background Lab summaries are not substituted post hoc.",
            "rescue_attempted":False,
        },
        "shard_receipts":receipts,
        "claim_boundary":[
            "The pooled base and ROI-stable analyses reuse the species-specific palette-JSD spatial statistic and geometry-preserving vertex-permutation null.",
            "Cross-species partial rank adjusts for all-100 sampled span and clear ROI/flip technical-failure rate.",
            "Panel P/N summaries are descriptive only; no panel-specific inferential p-values are added.",
            "No ecological predictor is introduced.",
            "The exact matched-background control is reported not evaluable rather than reconstructed after biological opening.",
        ],
    }
    (out/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
