#!/usr/bin/env python3
"""Full 7-environment-block prediction beyond original geodesic space AND dated phylogeny.

Direct LCVP tip taxa only, all environmental variables already attached to
original 42111 photographed source plants. Reuse verified joint Brownian-tree
and original photo-site spatial distance kernels; this wrapper freezes the
full 7-block inventory and tests whether exploratory held-out gains replicate
across both original 50/250 km spatial scales and genus/cell holdouts.

Never infer adaptation or a global 42111-taxon causal phylogenetic effect.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import audit_fcp_global42111_real_spatial_phylo_covariance_20261010 as kernel
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import BLOCKS

SCHEMA="fcp_global42111_all_abiotic_beyond_true_spatial_phylogenetic_covariance_v1"
CAPS=("250km","500km")
MODES=("climate","soil")
SCALES=("50km","250km")
HOLDS=("genus","true_photo_cell")
ORDERED_BLOCKS=(
    "elevation","temperature","precipitation","solar_radiation",
    "wind","vapor_pressure","soil"
)
EXPECTED_N={"250km":342,"500km":649}


def summary(raw:dict)->dict:
    if set(BLOCKS)!=set(ORDERED_BLOCKS):
        raise RuntimeError("Full source seven block abiotic inventory was altered")
    if raw.get("original_42111_photo_species_denominator")!=42111 or raw.get("original_one_photo_classifiable")!=18457:
        raise RuntimeError("The historical source flower-photo species denominator changed")
    if raw.get("full_872_1761_direct_tip_phylogenetic_coverage_HOLD") is not True:
        raise RuntimeError("Old full original phylogenetic sufficiency HOLD must persist")
    out={}
    for cap in CAPS:
        for mode in MODES:
            key=cap+"_"+mode
            r=raw["cohorts"][key]
            if r["source"]["n_direct_original_photo_taxa_pre_missingness"]!=EXPECTED_N[cap]:
                raise RuntimeError("Historical exact LCVP direct source species tips drifted")
            expected=set(ORDERED_BLOCKS)-({"soil"} if mode=="climate" else set())
            for scale in SCALES:
                for hold in HOLDS:
                    z=r["models"][scale+"__"+hold]
                    if z["status"]=="HOLD_DIRECT_LCVP_MODEL_COVERAGE":
                        continue
                    if z["status"]!="EXPLORATORY_EXACT_PHYLO_BROWNIAN_KERNEL_PLUS_SPATIAL_KERNEL_ENVIRONMENT":
                        raise RuntimeError("Unexpected phylogeny/spatial baseline status")
                    if set(z["source_environment_blocks"])!=expected:
                        raise RuntimeError("Some environmental predictor was omitted from real phylogeny/spatial model")
                    if not (z["modelled_true_phylogeny_BM_shared_branch_covariance"] and
                            z["modelled_original_photograph_spatial_exponential_covariance"]):
                        raise RuntimeError("A full environmental contrast lacks one mandatory covariance baseline")
            for block in ("ALL_ENV",*ORDERED_BLOCKS):
                if block=="soil" and mode=="climate":
                    continue
                stat=("all_abiotic_beyond_spatial_plus_real_phylogeny" if block=="ALL_ENV"
                      else "unique_"+block+"_beyond_space_real_phylogeny_and_other_abiotic")
                studies=[]
                for scale in SCALES:
                    for hold in HOLDS:
                        z=r["models"][scale+"__"+hold]
                        if z["status"]!="EXPLORATORY_EXACT_PHYLO_BROWNIAN_KERNEL_PLUS_SPATIAL_KERNEL_ENVIRONMENT":
                            studies.append({"scale_km":scale,"heldout":hold,"status":z["status"],"gain":None,"ci":None})
                        else:
                            v=z["conditional_feature_gains"][stat]
                            studies.append({
                                "scale_km":scale,"heldout":hold,"status":z["status"],
                                "gain":v["heldout_brier_reduction"],
                                "ci":v["original_group_cluster_95CI"],
                            })
                valid=[x for x in studies if x["gain"] is not None]
                positive_both_across_four=(len(valid)==4 and all(
                    x["gain"]>0 and x["ci"][0]>0 for x in valid))
                out[key+"__"+block]={
                    "original_direct_tips_before_environment_mask":EXPECTED_N[cap],
                    "n_environment_complete_direct_photo_species":r["source"]["n_direct_original_photo_taxa_full_environment_common_case"],
                    "n_tests_fitted":len(valid),
                    "fixed_repeated_spatial_and_heldout_results":studies,
                    "all_four_conditional_group_intervals_positive":positive_both_across_four,
                    "is_this_independent_replication":False,
                }
    return out


def run(source:pd.DataFrame,ledger:pd.DataFrame,tree250:Path,tree500:Path,
        *,nboot:int=199)->dict:
    if len(source)!=42111 or source.inat_taxon_id.nunique()!=42111:
        raise RuntimeError("Must use entire original 42111 single-photo source denominator")
    if len(ledger)!=1761 or ledger.inat_taxon_id.duplicated().any():
        raise RuntimeError("Original LCVP 1761 tip matching ledger not source-identical")
    raw=kernel.run(source,ledger,{250:tree250,500:tree500},
                   strict=True,nboot=nboot,blocks=BLOCKS)
    result={
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"EXPLORATORY_SEVEN_BLOCK_ABIOTIC_BEYOND_GEO_SPATIAL_AND_TRUE_DATED_LCVP",
        "original_global_source_species":42111,
        "original_classified_source_photo_species":18457,
        "source_actual_LCVP_direct_original_photo_species":EXPECTED_N,
        "all_seven_full_environment_blocks":{k:list(v) for k,v in BLOCKS.items()},
        "no_predictor_block_selected_after_colour_effects":True,
        "space_scales_km":[50,250],
        "heldout_axes":["genus","true_photo_cell"],
        "all_original_same_species_photo_source_and_independent_cohorts_untouched":True,
        "full_original_phylogenetic_inference_HOLD":True,
        "source_full_kernel_model":raw,
        "robustness_matrix":summary(raw),
        "decision_rule":"A block has an exploratory four-way positive conditional increment only when both 50/250km spatial kernels and genus/region-heldout fixed-fold group intervals are all positive, on the same specified direct-tip environment-complete cohort; not causal or independent replication.",
        "nonclaims":[
            "Direct-tip taxon subset covers just 39.2% and 36.9% of original 872/1761 parent photo cohorts",
            "The two parent cohorts are nested and environmental missingness selects original species",
            "Brownian dated LCVP and one exponential spatial covariance are chosen statistical structures, not complete causal spatial/phylogenetic deconfounding",
            "All seven correlated environmental block gains are retrospective and any positive conditional interval is unadjusted for many hypotheses",
            "True floral pigment genotype, within-population polymorphism and adaptive reproductive fitness were never measured",
            "Source 50/100km geographically matched flower-photo label permutation null previously did not support a robust local rainfall adaptation claim",
        ],
    }
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-expanded-breadth",required=True,type=Path)
    p.add_argument("--original-lcvp-ledger",required=True,type=Path)
    p.add_argument("--tree-250",required=True,type=Path)
    p.add_argument("--tree-500",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    raw=run(pd.read_csv(a.original_expanded_breadth,low_memory=False),
            pd.read_csv(a.original_lcvp_ledger,low_memory=False),
            a.tree_250,a.tree_500)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(raw,indent=2,sort_keys=True)+"\n")
    rows=[]
    for key,stat in raw["robustness_matrix"].items():
        for x in stat["fixed_repeated_spatial_and_heldout_results"]:
            rows.append({
                "cohort_and_block":key,
                "n_source_complete_direct_LCVP_tips":stat["n_environment_complete_direct_photo_species"],
                "spatial_kernel":x["scale_km"],"heldout_axis":x["heldout"],
                "gain_brier":x["gain"],
                "ci_low":x["ci"][0] if x["ci"] else None,
                "ci_high":x["ci"][1] if x["ci"] else None,
                "all_four_conditional_intervals_positive":stat["all_four_conditional_group_intervals_positive"],
            })
    pd.DataFrame(rows).to_csv(a.outdir/"full_seven_block_spatial_phylogenetic_robustness.csv",index=False)
    print(json.dumps({
        "schema":SCHEMA,
        "n_robustness_tests":len(rows),
        "n_blocks_positive_in_all_four":sum(z["all_four_conditional_group_intervals_positive"] for z in raw["robustness_matrix"].values()),
        "coverage":{k:v["source"] for k,v in raw["source_full_kernel_model"]["cohorts"].items()},
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()
