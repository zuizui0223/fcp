#!/usr/bin/env python3
"""Exact phylogeny coverage audit for 42111 flower-colour source and local taxa.

Compare previously used independent 34-tip JBI and ~709-tip H3a dated trees to
the exact 250/500 km source-conditioned historical photo list. A genus match
WITHOUT an exact species tip is NOT phylogenetic coverage; it only suggests
possible V.PhyloMaker2 backbone-grafting opportunity. No K/lambda/effect
estimated from preflight, and no source photo labels used in tree selection.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo
from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort, EXPECTED_ORIGINAL_N,
)
from audit_fcp_global42111_local_congeners_20261009 import photo_populations

TREE_SHA256={
    "h3a_s1":"ce4a2ede815dbe0eb2407d3d9ebb58d0218d85ece85ea1b9ab58e00e7dd84d19",
    "h3a_s2":"384c9d1bcde150ac200d7006703de55ebfee6f83fc9108bb9b7c02cee3107abb",
    "h3a_s3":"9a25eb77d87e7a698697f8c4835f8c24e98ee3c8bfc7c02d317f46090e382c26",
}
MAX_LOCAL_GROUP_KM=(250,500)
MIN_EXACT_TIPS=300
MIN_GENUS_TIP_PAIRS=30
MIN_EXACT_FRACTION=.5
SCHEMA="fcp_global42111_existing_tree_local_congener_coverage_v1"


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):h.update(buf)
    return h.hexdigest()


def key(raw:str)->str:
    # Conservative transformation of documented Newick underscore-space format.
    # No fuzzy synonym substitutions, genus fallback, or arbitrary taxon swaps.
    z=str(raw).strip().strip("'").strip('"').replace("_"," ")
    z=" ".join(z.split()).casefold()
    return z


def trees_from_sources(paths:dict[str,Path])->tuple[dict[str,set[str]],dict]:
    if set(paths)!=set(TREE_SHA256)|{"jbi34"}:
        raise ValueError("All three historical original dated H3a trees and JBI tree required")
    leaves={}
    metadata={}
    for name,path in paths.items():
        if name in TREE_SHA256 and sha(path)!=TREE_SHA256[name]:
            raise ValueError(f"{name}: source historical dated tree SHA mismatch")
        tree=Phylo.read(str(path),"newick")
        originals=[tip.name for tip in tree.get_terminals()]
        if any(v is None for v in originals):
            raise ValueError("Unlabeled historical phylogenetic terminal")
        canon=[key(v) for v in originals]
        if len(canon)!=len(set(canon)):
            raise ValueError("Ambiguous duplicate historical tips after name normalization")
        if len(canon)<10:
            raise ValueError("Tree topology too small to assess coverage")
        leaves[name]=set(canon)
        metadata[name]={
            "n_original_tree_terminal_species":len(canon),
            "n_unique_exact_source_tree_tip_names":len(leaves[name]),
            "n_unique_terminal_genera":len({v.split()[0] for v in canon}),
            "sha256":sha(path),
            "historical_backbone_is_not_42111_species_global_tree":True,
        }
    return leaves,metadata


def evaluate_coverage(cohort:pd.DataFrame,tips:dict[str,set[str]])->dict:
    if cohort.inat_taxon_id.duplicated().any():
        raise ValueError("Duplicate source species IDs")
    if len(cohort)<100:
        raise ValueError("Non-informative source cohort")
    names=cohort.species.map(key)
    if len(set(names))!=len(names):
        raise ValueError("Original selected source taxa collapse to same binomial tree tip")
    orig_genus=cohort.genus.astype(str).str.casefold()
    out={}
    for name,set_tips in tips.items():
        mask=names.isin(set_tips)
        in_tree=cohort.loc[mask].copy()
        genus_counts=in_tree.genus.astype(str).str.casefold().value_counts()
        with_congeners=in_tree.genus.astype(str).str.casefold().isin(genus_counts[genus_counts>=2].index)
        ancestor_genera={x.split()[0] for x in set_tips}
        group_counts=in_tree.genus_cell_id.value_counts()
        local_tree_groups=group_counts[group_counts>=2]
        out[name]={
            "n_original_photo_source_species":len(cohort),
            "n_exact_original_source_species_tree_tips":int(mask.sum()),
            "fraction_original_source_species_tree_tips":float(mask.mean()),
            "n_distinct_exact_tipped_genera":int(genus_counts.size),
            "n_exact_tipped_species_with_at_least_one_other_tipped_congener":int(with_congeners.sum()),
            "n_genera_with_at_least_two_exact_tree_tips":int(genus_counts.ge(2).sum()),
            "n_original_local_genus_cell_groups_with_at_least_two_exact_tips":int(len(local_tree_groups)),
            "n_original_source_species_with_genus_name_somewhere_in_tree":int(orig_genus.isin(ancestor_genera).sum()),
            "genus_overlap_is_not_species_phylogeny":True,
            "exact_tip_coverage_passed_conservative_gate":bool(
                mask.sum()>=MIN_EXACT_TIPS and mask.mean()>=MIN_EXACT_FRACTION
                and int(genus_counts.ge(2).sum())>=MIN_GENUS_TIP_PAIRS
            ),
            "original_tree_tip_selection_did_not_inspect_photo_colour_labels":True,
        }
    return out


def run(source:pd.DataFrame,paths:dict[str,Path],strict:bool=True)->dict:
    pop,den=photo_populations(source,strict=strict)
    taxa=pop["CLIMATE_ALL"]
    original_tips,tree_meta=trees_from_sources(paths)
    cohorts={}
    for cap in MAX_LOCAL_GROUP_KM:
        c=fixed_source_cohort(taxa,cap)
        if strict and len(c)!=EXPECTED_ORIGINAL_N[str(cap)]:
            raise RuntimeError(f"Fixed 250/500km original source classifiable cohort {cap} changed")
        cohorts[str(cap)]={
            "max_original_geographic_photo_group_diameter_km":cap,
            "n_original_species":len(c),
            "n_original_local_genus_cell_groups":int(c.genus_cell_id.nunique()),
            "name_matching":"EXACT_CANONICAL_BINOMIAL_WITH_UNDERSCORE_SPACE_ONLY",
            "existing_tree_coverage":evaluate_coverage(c,original_tips),
        }
    allpass=all(x["exact_tip_coverage_passed_conservative_gate"]
                for c in cohorts.values() for x in c["existing_tree_coverage"].values()
                if "h3a" in next(k for k,v in c["existing_tree_coverage"].items() if v is x))
    # More interpretable explicit gate: only S1-S3 have old tree coverage
    min_h3a=all(cohorts[c]["existing_tree_coverage"][f"h3a_s{i}"]["exact_tip_coverage_passed_conservative_gate"]
                 for c in ("250","500") for i in (1,2,3))
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"LEGACY_TREE_OVERLAP_ONLY_NO_PHYLOGENETIC_SIGNAL_ESTIMATED",
        "source_global_taxa":len(source),
        "original_source_classifiable_photo_taxa":den["source_classified_photo_taxa"],
        "original_climate_complete_photo_taxa":den["climate_eligible_classified_source_taxa"],
        "tree_metadata":tree_meta,
        "local_cohorts":cohorts,
        "precommitted_min_tree_exact_tip_species":MIN_EXACT_TIPS,
        "precommitted_min_tree_tip_fraction":MIN_EXACT_FRACTION,
        "precommitted_min_two_species_genus_groups":MIN_GENUS_TIP_PAIRS,
        "old_h3a_trees_support_this_global_test":bool(min_h3a),
        "decision":"EXISTING_H3A_PHYLOGENY_CAN_BE_EVALUATED_WITH_EXACT_TIPS" if min_h3a else "HOLD_EXISTING_TREE_COVERAGE_REQUIRE_NEW_OUTCOME_BLIND_DATED_TREE",
        "original_source_tree_leaf_labels_and_photo_taxa_left_unchanged":True,
        "no_new_colour_outcomes_or_prospective_taxa_opened":True,
        "hard_nonclaims":[
            "A genus name match is not an exact phylogenetic species tip or a dated lineage",
            "A backbone-based synthetic placement without within-genus sequence relationships does not establish subgenus-level relatedness",
            "Previously archived H3a source trees were built for a separate 732-species D study and are not globally representative",
            "No phylogenetic environmental or trait-effect coefficient is estimated by this source coverage gate",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    for name in ("h3a_s1","h3a_s2","h3a_s3","jbi34"):
        p.add_argument("--"+name.replace("_","-"),required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    df=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    files={n:getattr(a,n) for n in ("h3a_s1","h3a_s2","h3a_s3","jbi34")}
    receipt=run(df,files)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,indent=2,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
