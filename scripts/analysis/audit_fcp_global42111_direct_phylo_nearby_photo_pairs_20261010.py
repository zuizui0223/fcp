#!/usr/bin/env python3
"""True local/photo + directly tipped LCVP phylogenetic comparability preflight.

The original 872/1761 source congeneric photo cohorts remain fixed. Only
names/IDs/site coordinates/tree tip identity are evaluated. No flower-colour
or rainfall effect, no genotype or adaptive inference, no synthetic grafts.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort,EXPECTED_ORIGINAL_N
)
from audit_fcp_global42111_local_congeners_20261009 import photo_populations
from audit_fcp_global42111_microspatial_photo_shuffle_20261009 import (
    microgroups,great_circle_matrix_km
)

DISTANCES=(50,100)
SOURCE_COARSE=(250,500)
KNOWN_EXACT_BACKBONE={"250":342,"500":649}
MIN_DIRECT_PHYLO_CONNECTED_TAXA=100
MIN_CONGENERIC_GENERA=20
MIN_LOCAL_MICROGROUPS=30
MIN_PHYLO_DISTANCE_VARIABLE_GROUPS=10
MIN_PHYLO_DISTANCE_VARIABLE_TAXA=40
SCHEMA="fcp_42111_exact_LCVP_tipped_local_congeners_50_100km_coverage_v1"


def canonical(s:str)->str:
    return " ".join(str(s).strip().replace("_"," ").split()).casefold()


def direct_taxa(d:pd.DataFrame,ledger:pd.DataFrame,tree_path:Path,cap:int,
               *,strict:bool=True)->tuple[pd.DataFrame,object]:
    if (strict and len(ledger)!=1761) or ledger.inat_taxon_id.duplicated().any():
        raise ValueError("Historical taxonomy and direct LCVP ledger changed")
    selected=ledger.loc[ledger.direct_lcvp_backbone_tip.astype(bool)].copy()
    selected=selected.loc[selected.in_fixed_250km_original_source.astype(bool)] if cap==250 else selected
    if strict and len(selected)!=KNOWN_EXACT_BACKBONE[str(cap)]:
        raise ValueError("Original source exact LCVP tips cannot be substituted")
    taxon=set(selected.inat_taxon_id.astype(int))
    x=d.loc[d.inat_taxon_id.isin(taxon)].copy().reset_index(drop=True)
    if len(x)!=len(selected) or x.inat_taxon_id.duplicated().any():
        raise ValueError("No original source photo IDs can be omitted or duplicated")
    tree=Phylo.read(str(tree_path),"newick")
    tips={canonical(a.name):a.name for a in tree.get_terminals()}
    if len(tips)!=len(tree.get_terminals()):
        raise ValueError("Ambiguous original dated LCVP tip names")
    if len(tips)!=len(x):
        raise ValueError("Original directly tipped subtree != selected taxa")
    original_names=x.species.map(canonical)
    if set(original_names)!=set(tips):
        raise ValueError("Frozen photo species exact binomial names mismatch LCVP tree tips")
    x["original_direct_LCVP_tip"]=original_names.map(tips)
    if x.original_direct_LCVP_tip.isna().any():
        raise ValueError("Lost original direct matched tree label")
    return x,tree


def local_coverage(d:pd.DataFrame,tree,maxdiameter:int)->dict:
    if d.index.tolist()!=list(range(len(d))):
        raise ValueError("Need source-verified zero-based original direct-tree photo population")
    original=microgroups(d,float(maxdiameter))
    taxa=set()
    genus=set()
    cells=set()
    ngroup=0
    pair_count=0
    patristic=[]
    n_3plus_groups=0
    n_varied_groups=0
    varied_taxa=set()
    for ids in original:
        if len(ids)<2:continue
        selected=d.iloc[ids]
        if selected.genus.nunique()!=1 or selected.genus_cell_id.nunique()!=1:
            raise RuntimeError("Mixed original source genus/cell in phylogenetic local pair")
        ngroup+=1
        taxa.update(selected.inat_taxon_id.astype(int).tolist())
        genus.add(str(selected.genus.iloc[0]))
        cells.add(str(selected.genus_cell_id.iloc[0]))
        names=selected.original_direct_LCVP_tip.to_list()
        if len(ids)>=3:
            n_3plus_groups+=1
        within_distances=[]
        for a,b in itertools.combinations(names,2):
            distance=tree.distance(a,b)
            if not np.isfinite(distance) or distance<0:
                raise RuntimeError("Non-finite or negative source tree path length")
            patristic.append(float(distance))
            within_distances.append(float(distance))
            pair_count+=1
        if len(ids)>=3 and len(set(round(k,5) for k in within_distances))>=2:
            n_varied_groups+=1
            varied_taxa.update(selected.inat_taxon_id.astype(int).tolist())
    enough=len(taxa)>=MIN_DIRECT_PHYLO_CONNECTED_TAXA and len(genus)>=MIN_CONGENERIC_GENERA and ngroup>=MIN_LOCAL_MICROGROUPS
    phylo_resolution=(n_varied_groups>=MIN_PHYLO_DISTANCE_VARIABLE_GROUPS and
                      len(varied_taxa)>=MIN_PHYLO_DISTANCE_VARIABLE_TAXA)
    return {
        "max_original_photo_microgeographical_diameter_km":maxdiameter,
        "n_original_direct_backbone_photo_species":len(d),
        "n_total_source_original_direct_tip_microgroups":len(original),
        "n_local_multispecies_direct_tip_groups":ngroup,
        "n_original_species_in_at_least_one_local_direct_phylo_pair":len(taxa),
        "n_source_original_congeneric_genera_with_local_direct_pairs":len(genus),
        "n_distinct_original_genus_by_geo_cell_local_groups_with_tips":len(cells),
        "n_original_congeneric_phylogenetic_tip_pairs_in_local_neighborhoods":pair_count,
        "median_direct_tip_LCVP_patristic_distance_backbone_units":float(np.median(patristic)) if patristic else None,
        "n_zero_patristic_distance_pairs":int(sum(v<1e-9 for v in patristic)),
        "n_distinct_LCVP_path_distances_rounded_5_decimal":int(len(set(round(k,5) for k in patristic))),
        "min_LCVP_path_distance":float(min(patristic)) if patristic else None,
        "max_LCVP_path_distance":float(max(patristic)) if patristic else None,
        "std_LCVP_path_distance":float(np.std(patristic)) if len(patristic)>=2 else None,
        "n_local_groups_with_three_or_more_direct_tips":n_3plus_groups,
        "n_local_groups_three_plus_tips_with_two_or_more_distinct_LCVP_path_distances":n_varied_groups,
        "n_original_source_taxa_in_phylo_path_varied_local_groups":len(varied_taxa),
        "min_groups_with_identifiable_withingroup_tree_distance_variation":MIN_PHYLO_DISTANCE_VARIABLE_GROUPS,
        "min_taxa_in_identifiable_withingroup_tree_distance_variation":MIN_PHYLO_DISTANCE_VARIABLE_TAXA,
        "within_local_group_phylogenetic_distance_variation_status":(
            "EXPLORATORY_WITHINGROUP_PHYLO_DISTANCE_VARIATION_PASS" if phylo_resolution
            else "HOLD_INSUFFICIENT_WITHINGROUP_DATED_PHYLO_DISTANCE_VARIATION"
        ),
        "all_phylogenetic_pairs_both_direct_real_backbone_tips":True,
        "microgroup_selection_uses_no_photographed_flower_colour":True,
        "minimum_direct_photo_taxa":MIN_DIRECT_PHYLO_CONNECTED_TAXA,
        "minimum_independent_congeneric_genera":MIN_CONGENERIC_GENERA,
        "minimum_local_phylogenetic_groups":MIN_LOCAL_MICROGROUPS,
        "readiness":"EXPLORATORY_DIRECT_TIP_SPATIAL_PHYLO_COVERAGE_PASS" if enough else "HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLOGENETIC_LOCAL_COMPARISONS",
    }


def audit(original:pd.DataFrame,ledger:pd.DataFrame,
          trees:dict[int,Path],strict:bool=True)->dict:
    cohorts,counts=photo_populations(original,strict=strict)
    result={}
    for cap in SOURCE_COARSE:
        full=fixed_source_cohort(cohorts["CLIMATE_ALL"],cap)
        if strict and len(full)!=EXPECTED_ORIGINAL_N[str(cap)]:
            raise RuntimeError("Original photo-cohort source identification changed")
        selected,tree=direct_taxa(full,ledger,trees[cap],cap,strict=strict)
        result[str(cap)]={
            "full_original_frozen_photo_species":len(full),
            "n_exact_direct_LCVP_backbone_photo_species":len(selected),
            "n_actual_induced_source_tree_tips":len(tree.get_terminals()),
            "locality_versus_phylogeny_coverage":{
                str(d):local_coverage(selected,tree,d) for d in DISTANCES
            },
        }
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"EXACT_DATED_BACKBONE_TIP_PAIR_ORIGINAL_PHOTO_SPATIAL_OPPORTUNITY_ONLY",
        "original_global_source_photo_taxa":len(original),
        "original_colour_classifiable":counts["source_classified_photo_taxa"],
        "historical_fixed_caps_km":list(SOURCE_COARSE),
        "new_local_comparison_diameters_km":list(DISTANCES),
        "uses_ONLY_induced_direct_73420_tip_backbone_subtrees":True,
        "original_source_colour_never_used_to_select_subgroups":True,
        "photo_colour_climate_phylogenetic_coefficient_computed":False,
        "no_new_source_images_photos_or_independent_taxa":True,
        "cohorts":result,
        "hard_nonclaims":[
            "The original LCVP backbone is a phylogenetic hypothesis with imperfect species-level resolution",
            "Many originally photographed congeneric species are NOT direct tips and are excluded from source-direct-tree coverage, creating selection bias",
            "A <100km photo-site matched congeneric pair is not a known breeding population",
            "Any local source phylogenetic coverage is not an evolutionary causal mechanism or climate adaptation test",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-breadth-abiotic",required=True,type=Path)
    ap.add_argument("--original-LCVP-tip-ledger",required=True,type=Path)
    ap.add_argument("--tree-250",required=True,type=Path)
    ap.add_argument("--tree-500",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    data=pd.read_csv(args.original_breadth_abiotic,low_memory=False)
    ledger=pd.read_csv(args.original_LCVP_tip_ledger,low_memory=False)
    report=audit(data,ledger,{250:args.tree_250,500:args.tree_500})
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
