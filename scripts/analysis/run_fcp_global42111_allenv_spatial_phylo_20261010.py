#!/usr/bin/env python3
"""All environment blocks beyond spatial and directly dated phylogenetic covariance."""
import argparse
import json
from pathlib import Path
import pandas as pd
import audit_fcp_global42111_real_spatial_phylo_covariance_20261010 as baseline
import compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 as predictors

SCHEMA="fcp_global42111_all_abiotic_beyond_real_spatial_phylo_v1"

def analyze(original, tips, t250, t500, nboot=199):
    blocks=predictors.BLOCKS
    expected=("elevation","temperature","precipitation","solar_radiation","wind","vapor_pressure","soil")
    if tuple(blocks)!=expected or sum(len(v) for v in blocks.values())!=36:
        raise ValueError("Predeclared seven full abiotic blocks changed")
    if len(original)!=42111 or original.inat_taxon_id.nunique()!=42111:
        raise ValueError("Original photo species denominator altered")
    if int(original.measurement_status.eq("classified_four_state_morph").sum())!=18457:
        raise ValueError("Original photo four-colour denominator altered")
    if not all(k in original for x in blocks.values() for k in x):
        raise ValueError("Missing a source-environment variable")
    old=baseline.BLOCKS
    try:
        baseline.BLOCKS=blocks
        original_result=baseline.run(original,tips,{250:t250,500:t500},nboot=nboot)
    finally:
        baseline.BLOCKS=old
    assert original_result["full_872_1761_direct_tip_phylogenetic_coverage_HOLD"]
    return {"schema":SCHEMA,"status":"EXPLORATORY_SOURCE_DIRECT_TREE_SPATIAL_7_ENVIRONMENT_BLOCKS",
      "n_original_species":42111,"n_original_colour_classifiable":18457,
      "environment_blocks":{k:list(v) for k,v in blocks.items()},
      "n_named_predictors":36,
      "original_spatial_and_phylogenetic_covariance_used":True,
      "direct_tree_selection_bias_and_full_source_coverage_HOLD":True,
      "original_labels_and_images_unchanged":True,
      "result":original_result,
      "nonclaims":["Tree directly covers only 342 and 649 taxa of source local cohorts",
                   "Predictive gains do not imply adaptive flower-colour selection",
                   "Posthoc block ablation and fixed-fold bootstrap are exploratory"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--tips",required=True,type=Path)
    p.add_argument("--tree250",required=True,type=Path)
    p.add_argument("--tree500",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    out=analyze(pd.read_csv(a.source,low_memory=False),
                pd.read_csv(a.tips,low_memory=False),a.tree250,a.tree500)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"schema":SCHEMA,"status":out["status"],
                      "cohorts":{k:{"support":v["source"],"gains":{
                          name:m.get("conditional_feature_gains")
                          for name,m in v["models"].items()}}
                          for k,v in out["result"]["cohorts"].items()}}),flush=True)

if __name__=="__main__":
    main()
