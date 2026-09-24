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
SPEC=importlib.util.spec_from_file_location("fcp_v2_mv4_core",MODULE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load MV4 core")
CORE=importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name]=CORE
SPEC.loader.exec_module(CORE)

SHARDS=20


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--biological-seal-dir",type=Path,required=True)
    p.add_argument("--shard-index",type=int,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    if not 0<=args.shard_index<SHARDS:
        raise ValueError("invalid MV4 shard index")

    summary=json.loads((args.biological_seal_dir/"biological_summary.json").read_text())
    if summary.get("status")!="BIOLOGICAL_PASS_B_SEALED" or int(summary.get("rows",-1))!=40000:
        raise RuntimeError("biological seal not terminal")
    if summary.get("source_identity",{}).get("replacement_rows")!=0:
        raise RuntimeError("replacement detected")

    df=pd.read_csv(args.biological_seal_dir/"biological_table.csv.gz",compression="gzip",dtype={"measurement_id":str})
    species=sorted(df["species"].astype(str).unique())
    if len(species)!=400:
        raise RuntimeError("species denominator drift")
    selected=species[args.shard_index::SHARDS]
    sub=df.loc[df["species"].astype(str).isin(set(selected))].copy()
    if len(sub)!=100*len(selected):
        raise RuntimeError("shard raw-row denominator drift")

    base=CORE.run_condition(sub,"base")
    stable=CORE.run_condition(sub,"roi_stable_only")

    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    base["species_table"].to_csv(out/"base_species.csv",index=False,lineterminator="\n")
    stable["species_table"].to_csv(out/"roi_stable_species.csv",index=False,lineterminator="\n")
    np.savez_compressed(out/"nulls.npz",base=base["null"],roi_stable=stable["null"])
    receipt={
        "schema":"fcp_v2_mv4_shard_v1",
        "status":"MV4_SHARD_COMPLETE",
        "shard_index":args.shard_index,
        "shards":SHARDS,
        "raw_species":len(selected),
        "raw_rows":len(sub),
        "base_eligible_species":base["eligible_species"],
        "roi_stable_eligible_species":stable["eligible_species"],
        "selected_species":selected,
        "per_species_null_replicates":CORE.N_PERM,
        "matched_background_exact_control_evaluable":False,
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
