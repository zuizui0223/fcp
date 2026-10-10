#!/usr/bin/env python3
"""Seven abiotic blocks added to TRUE geodesic + dated LCVP covariance.

Exact original 42111 photographs; only 342/649 directly tipped source species
and their 250/500km historical neighbourhood cohorts. No full-cohort
phylogenetic effect; no main/old 1499 or prospective 2000+730 changes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

import audit_fcp_global42111_real_spatial_phylo_covariance_20261010 as P
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import BLOCKS

SCHEMA="fcp_42111_direct_LCVP_spatial_joint_all_abiotic_source_v1"


def explore(source,ledger,trees,*,strict=True,nboot=199):
    original_blocks=dict(BLOCKS)
    if list(original_blocks)!=["elevation","temperature","precipitation","solar_radiation","wind","vapor_pressure","soil"]:
        raise ValueError("Changed all-abiotic ecological hypotheses after original outcome exposure")
    if sum(len(x) for x in original_blocks.values())!=36:
        raise ValueError("Exactly 36 declared physical source predictors required")
    receipt=P.run(source,ledger,trees,strict=strict,nboot=nboot,blocks=original_blocks)
    receipt.update({
        "schema":SCHEMA,
        "status":"RETROSPECTIVE_FULL_ENV_BEYOND_REAL_SPATIAL_AND_DATED_LCVP_COVARIANCE",
        "all_seven_predictor_blocks":{k:list(v) for k,v in original_blocks.items()},
        "all_36_predictors_scanned_without_outcome_selected_drop":True,
        "same_one_photo_source_taxa_and_complete_cases_per_environment_comparison":True,
        "full_42111_phylogenetically_controlled_result_NOT_identified":True,
        "previous_photo_source_fine_geographical_colour_label_nulls_not_overridden":True,
    })
    return receipt


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expanded-breadth",type=Path,required=True)
    p.add_argument("--original-LCVP-tip-ledger",type=Path,required=True)
    p.add_argument("--tree-250",type=Path,required=True)
    p.add_argument("--tree-500",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    df=pd.read_csv(a.expanded_breadth,low_memory=False)
    tips=pd.read_csv(a.original_LCVP_tip_ledger,low_memory=False)
    result=explore(df,tips,{250:a.tree_250,500:a.tree_500})
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":SCHEMA,"denominator":result["original_42111_photo_species_denominator"],
        "summary":{k:{
            "source":v["source"],"block_increment":{
                s:{kk:vv["heldout_brier_reduction"] for kk,vv in o.get("conditional_feature_gains",{}).items()}
                for s,o in v["models"].items()}
        } for k,v in result["cohorts"].items()}},sort_keys=True),flush=True)


if __name__=="__main__":
    main()
