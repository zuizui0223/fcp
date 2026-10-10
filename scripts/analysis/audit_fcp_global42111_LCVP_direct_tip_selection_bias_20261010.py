#!/usr/bin/env python3
"""Outcome-source selection audit: LCVP direct tip vs excluded photo species.

The original 872/1761 source species are all four-colour classifiable and
climate complete, but only 342/649 were direct exact LCVP tips. Before
generalizing the direct-tip ecology result, compare the already-frozen original
photo colour composition, true-site climate and geographic exposure of
directly-tipped vs excluded species. No model selection, null p-value or
phenotype reclassification. Pure missing-tree-coverage diagnosis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort, EXPECTED_ORIGINAL_N
)
from audit_fcp_global42111_local_congeners_20261009 import photo_populations

SOURCE_TAXA=42111
SOURCE_CLASSIFIED=18457
SAMPLE={"250":(872,342),"500":(1761,649)}
FEATURES=("abs_latitude","wc_elevation_m","wc_bio1","wc_bio5","wc_bio12","wc_bio15")
STATES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="fcp_global42111_original_LCVP_tip_selection_covariate_bias_v1"


def cover(cohort:pd.DataFrame,ledger:pd.DataFrame,cap:int,*,strict=True)->dict:
    label=str(cap)
    if label not in SAMPLE:
        raise ValueError("Unknown original source locality cohort")
    n,ntip=SAMPLE[label]
    if strict and len(cohort)!=n:
        raise ValueError("Original 250/500 source photo cohort changed")
    must={"inat_taxon_id","direct_lcvp_backbone_tip","in_fixed_250km_original_source"}
    if not must.issubset(ledger):
        raise ValueError("Original current taxon backbone matching ledger missing")
    if ledger.inat_taxon_id.duplicated().any():
        raise ValueError("Original source taxon duplicated")
    lookup=ledger.set_index("inat_taxon_id")
    if not set(cohort.inat_taxon_id).issubset(lookup.index):
        raise RuntimeError("Some original source photo species missing taxonomy ledger")
    flags=lookup.loc[cohort.inat_taxon_id,"direct_lcvp_backbone_tip"].astype(bool).to_numpy()
    if strict and flags.sum()!=ntip:
        raise RuntimeError("Original direct LCVP tip source selection drifted")
    if cap==250 and not lookup.loc[cohort.inat_taxon_id,"in_fixed_250km_original_source"].astype(bool).all():
        raise RuntimeError("Original 250km photo source out of nested geographic cohort")
    z=cohort.copy()
    z["source_direct_tree_tip"]=flags
    if len(z)!=z.inat_taxon_id.nunique() or not z.morph.isin(STATES).all():
        raise ValueError("Original four-colour species-equal source no longer classified")
    if not z[list(FEATURES)].notna().all().all():
        raise RuntimeError("Previously frozen climate-complete source has missing covariates")
    out={}
    overall=z.genus.nunique()
    for status,value in (("DIRECT_LCVP",True),("NOT_DIRECT_LCVP",False)):
        d=z.loc[z.source_direct_tree_tip==value]
        colors=d.morph.value_counts().reindex(STATES,fill_value=0)
        out[status]={
            "n_source_original_species":len(d),
            "n_original_taxonomic_genera":int(d.genus.nunique()),
            "n_original_photo_geographic_cells":int(d.photo_cell_162.nunique()),
            "source_photo_four_colour_counts":{s:int(colors[s]) for s in STATES},
            "source_photo_four_colour_fractions":{s:float(colors[s]/len(d)) for s in STATES},
            "real_photo_site_feature_medians":{k:float(d[k].median()) for k in FEATURES},
            "n_source_photo_species_with_original_genus_present_in_both_direct_statuses":int(d.genus.isin(
                set(z.loc[z.source_direct_tree_tip!=value,"genus"])).sum()),
        }
    b=z.loc[z.source_direct_tree_tip]
    m=z.loc[~z.source_direct_tree_tip]
    contrasts={}
    for k in FEATURES:
        x=b[k].to_numpy(float)
        y=m[k].to_numpy(float)
        pooled=np.sqrt((x.var(ddof=1)+y.var(ddof=1))/2)
        smd=float((x.mean()-y.mean())/pooled) if pooled>0 else None
        contrasts[k]={
            "direct_minus_not_direct_mean":float(x.mean()-y.mean()),
            "standardized_mean_difference":smd,
            "median_direct":float(np.median(x)),
            "median_not_direct":float(np.median(y)),
        }
    color_max=float(max(abs(out["DIRECT_LCVP"]["source_photo_four_colour_fractions"][c]-
                            out["NOT_DIRECT_LCVP"]["source_photo_four_colour_fractions"][c])
                        for c in STATES))
    max_smd=max(abs(v["standardized_mean_difference"] or 0.) for v in contrasts.values())
    return {
        "original_cohort_max_photo_group_distance_km":cap,
        "n_original_source_photo_species":len(z),
        "n_direct_original_LCVP_tip_species":int(b.shape[0]),
        "n_not_direct_original_LCVP_tip_species":int(m.shape[0]),
        "original_direct_tree_tip_fraction":float(b.shape[0]/len(z)),
        "n_original_nominal_genera":int(overall),
        "n_orig_genera_with_both_direct_and_untipped_photo_species":int(len(set(b.genus)&set(m.genus))),
        "direct_tip_strata":out,
        "original_photo_site_environment_standardized_difference_direct_minus_untipped":contrasts,
        "max_abs_std_difference_across_six_frozen_photo_geography_climate_features":float(max_smd),
        "max_abs_original_colour_fraction_difference":color_max,
        "descriptive_imbalance_at_abs_standardized_difference_0p10":bool(max_smd>=.1),
        "direct_backbone_tip_selection_did_not_use_original_photo_colour_results":True,
    }


def run(original:pd.DataFrame,ledger:pd.DataFrame,*,strict=True)->dict:
    populations,original_stats=photo_populations(original,strict=strict)
    if strict and len(original)!=SOURCE_TAXA:
        raise RuntimeError("Original global 42111 photo breadth incomplete")
    cohorts={}
    for cap in (250,500):
        full=fixed_source_cohort(populations["CLIMATE_ALL"],cap)
        cohorts[str(cap)]=cover(full,ledger,cap,strict=strict)
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"DIRECT_DATED_PHYLOGENY_TIP_SELECTION_BIAS_DIAGNOSTIC_ONLY",
        "source_global_original_photo_species":len(original),
        "source_global_original_classified_photo_species":original_stats["source_classified_photo_taxa"],
        "original_local_source_caps_km":[250,500],
        "source_full_original_250_and_500_cohort_phylo_control_HOLD":True,
        "no_new_flower_photo_pixels_classifications_or_environment_lookups":True,
        "original_tree_graft_species_not_treated_as_direct_tip":True,
        "cohorts":cohorts,
        "hard_nonclaims":[
            "Photo classifiability selection occurs before direct phylogenetic tip matching and can bias both included/excluded strata",
            "Descriptive standardized mean difference is not a test of true adaptation or species origin",
            "Direct tip species are an incomplete, nonrandom subset of the 872 and 1761 original source taxa",
            "LCVP backbone species inclusion may depend on taxonomic and literature sampling, not floral biology",
            "Nested geographical source cohorts are not independently replicated validation datasets",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",type=Path,required=True)
    p.add_argument("--original-LCVP-tip-ledger",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    source=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    ledger=pd.read_csv(a.original_LCVP_tip_ledger,low_memory=False)
    receipt=run(source,ledger)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
