#!/usr/bin/env python3
"""Authenticate and preserve real original FCP 149900+wave4 49800 photographic assays.

This module only validates completed source results. It never downloads JPEGs,
reclassifies flowers or recalculates an ecological estimate. Historical 50km
p=.005 is a 199-null-world Monte Carlo result, and 50km remains post-hoc.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

EXPECTED_OLD={
 "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
 "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
 "prospective":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
OLD_50={
 "discovery":(166,0.020529254583812922,.005),
 "validation":(181,0.018672971642749295,.005),
 "prospective":(204,0.01468491968437185,.005),
}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as fd:
        for data in iter(lambda:fd.read(1<<20),b""):
            h.update(data)
    return h.hexdigest()


def measured_bool(col:pd.Series)->pd.Series:
    return col.astype(str).str.lower().isin(["true","1","yes"])


def audit(wave4:Path,combined:Path)->dict:
    mfile=wave4/"wave4-real/wave4_original_method_measured_photos.csv"
    term=wave4/"wave4-real/terminal_result.json"
    ifile=wave4/"wave4-statistics/wave4_original_149900_method_analysis_result.json"
    cfile=combined/"fcp_original149900_plus_wave4_same_method_photos.csv.gz"
    crfile=combined/"expanded_data_receipt.json"
    sfile=combined/"expanded_fcp_original_50km_result.json"
    pfile=combined/"original_149900_plus_wave4_species_50km.csv"
    expected=[mfile,term,ifile,cfile,crfile,sfile,pfile]
    if any(not path.is_file() for path in expected):
        raise RuntimeError("both pinned original measured artifacts not fully recovered")
    terminal=json.loads(term.read_text())
    local=json.loads(ifile.read_text())
    receipt=json.loads(crfile.read_text())
    pooled=json.loads(sfile.read_text())
    if (terminal.get("schema")!="fcp_same_149900_method_wave4_real_photo_terminals_v1"
        or terminal.get("original_results_changed") is not False
        or terminal.get("measurement_original_fixed_colour_worker") is not True
        or terminal.get("terminal_partitions")!=256
        or terminal.get("wave4_100photo_species")!=498
        or terminal.get("wave4_new_photo_ids")!=49800
        or terminal.get("measured_sha256")!=sha(mfile)
        or local.get("schema")!="fcp_same_149900_method_wave4_statistical_replication_v1"
        or local.get("original_frozen_results_unchanged") is not True
        or local.get("wave4_measured_table_sha256")!=sha(mfile)):
        raise RuntimeError("real 49800-photo ROI/D/H1/W SHA or original model contract drift")
    w=pd.read_csv(mfile,low_memory=False)
    if (len(w)!=49800 or w.photo_id.nunique()!=49800
        or w.inat_taxon_id.nunique()!=498
        or not w.groupby("inat_taxon_id").size().eq(100).all()
        or measured_bool(w.global_classifiable).sum()!=terminal["wave4_classifiable_photo_rows"]):
        raise RuntimeError("real measured wave4 100 photos/species, IDs or 4-colour labels invalid")
    if (receipt.get("schema")!="fcp_original_149900_plus_same_method_wave4_merged_v1"
        or receipt.get("original_study_unchanged") is not True
        or receipt.get("append_only") is not True
        or receipt.get("four_state_method_relabelled") is not False
        or receipt.get("existing_50km_null_reweighted") is not False
        or receipt.get("merged_csv_sha256")!=sha(cfile)
        or receipt.get("input_sha256",{}).get("wave4")!=sha(mfile)
        or any(receipt.get("input_sha256",{}).get(k)!=v for k,v in EXPECTED_OLD.items())
        or receipt.get("all_original_photo_rows")!=199700
        or receipt.get("original_source_rows")!=149900
        or receipt.get("new_original_photo_rows")!=49800):
        raise RuntimeError("199700-photo original cohort preservation or SHA failed")
    c=pd.read_csv(cfile,low_memory=False)
    if (len(c)!=199700 or c.photo_id.nunique()!=len(c)
        or c.inat_taxon_id.nunique()!=1997
        or not c.groupby("inat_taxon_id").size().eq(100).all()
        or set(c.fcp_source_cohort)!=set(EXPECTED_OLD)|{"wave4"}
        or len(c.loc[c.fcp_source_cohort.eq("wave4")])!=49800):
        raise RuntimeError("old149900+new49800 observed unique photo and original species identities lost")
    if set(c.loc[c.fcp_source_cohort.eq("wave4"),"photo_id"].astype(str)) != set(w.photo_id.astype(str)):
        raise RuntimeError("published wave4 original pixel results no longer identical to appended photos")
    if (pooled.get("schema")!="fcp_149900_plus_wave4_original_50km_replayed_all_species_v1"
        or pooled.get("all_original_frozen_cohort_50km_values_reproduced") is not True
        or pooled.get("total_photo_rows")!=199700):
        raise RuntimeError("expanded original 50km vertex-permutation null not fully verified")
    results=pooled["results"]
    for group,(n,effect,p) in OLD_50.items():
        old=results["cohorts"][group]["radii"]["50"]
        if (old.get("n_species")!=n
            or abs(old.get("mean_depletion",0)-effect)>1e-12
            or old.get("mean_depletion_p_upper")!=p):
            raise RuntimeError("one of three original 149900 50km effects was overwritten")
    from_wave4=local["spatial_original_method"]["50"]
    appended=results["cohorts"]["wave4"]["radii"]["50"]
    if (from_wave4["n_species"]!=190
        or appended["n_species"]!=190
        or abs(from_wave4["mean_depletion"]-appended["mean_depletion"])>1e-12
        or abs(from_wave4["mean_depletion"]-.017293596294403593)>1e-12
        or from_wave4["mean_depletion_p_upper"]!=.005):
        raise RuntimeError("original wave4 true 50km effect or eligible-species support changed")
    stat=results["pooled"]["50"]
    if (stat.get("n_species")!=741
        or abs(stat.get("mean_depletion",0)-.01763720788097046)>1e-12
        or stat.get("mean_depletion_p_upper")!=.005):
        raise RuntimeError("199700 original 50km pooled Monte Carlo statistic mismatch")
    h1=local["h1_original_method"]["primary"]
    W=local["h2_original_W_method"]["primary_0_10"]
    if (h1["partitions"]!=200
        or abs(h1["rho_median"]-.8055963704422957)>1e-12
        or W["species"]!=134
        or abs(W["observed_W"]-.5368155693013797)>1e-12
        or W["structured_null_upper_p"]!=.001
        or W["structured_null_replicates"]!=999):
        raise RuntimeError("unchanged wave4 observer H1 / fixed white-axis W estimates missing")
    # The original wave4 worker HEAD predates the reporting-only PR #226.
    # Its top-level 50km n was written as zero via a numeric-key access.
    # Preserve original bytes and use the authentic nested 50km n=190.
    if local.get("wave4_50km_evaluable_species")!=0:
        raise RuntimeError("unmodified wave4 pre-fix 50km display fingerprint unexpectedly changed")
    species=pd.read_csv(pfile,low_memory=False)
    if len(species)!=results["pooled_n_spatial_eligible_species"]:
        raise RuntimeError("individual 50km original species results no longer match summary")
    return {
       "schema":"fcp_real_199700_original_assayed_photos_original_H1_50km_W_release_receipt_v1",
       "old_149900_measured_photos_preserved":True,
       "wave4_498_genuinely_measured_taxa":498,
       "wave4_real_original_photo_measurement_rows":49800,
       "wave4_four_state_classified_photos":int(measured_bool(w.global_classifiable).sum()),
       "all_original_measured_photo_rows":199700,
       "all_original_100photo_taxa":1997,
       "wave4_original_50km_evaluable_species":190,
       "wave4_original_50km_depletion":from_wave4["mean_depletion"],
       "wave4_original_50km_p_upper":from_wave4["mean_depletion_p_upper"],
       "pooled_original_50km_evaluable_species":741,
       "pooled_original_50km_depletion":stat["mean_depletion"],
       "pooled_original_50km_p_upper":stat["mean_depletion_p_upper"],
       "wave4_H1_observer_disjoint_rho_median":h1["rho_median"],
       "wave4_white_axis_W_n_taxa":W["species"],
       "wave4_white_axis_W_observed":W["observed_W"],
       "wave4_white_axis_W_null_median":W["structured_null_summary"]["q50"],
       "wave4_white_axis_W_p_upper":W["structured_null_upper_p"],
       "pre_reporting_fix_50km_zero_field_is_not_the_actual_spatial_result":True,
       "original_roi_v4_EfficientSAM_unchanged":True,
       "original_199_permutation_50km_null_unchanged":True,
       "original_999_structured_white_axis_null_unchanged":True,
       "original_50km_radius_posthoc":True,
       "original_iNaturalist_source_not_an_independent_photo_provider":True,
       "no_original_JPEG_republished":True,
       "files_sha256":{
           str(path.relative_to(wave4)) if path.is_relative_to(wave4)
           else "combined/"+path.name:sha(path) for path in expected
       },
       "source_original_pixel_run_id":38992543115,
       "source_original_199700_append_run_id":39023781282,
    }


def main()->None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--wave4",required=True,type=Path)
    parser.add_argument("--combined",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    report=audit(args.wave4,args.combined)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_REAL_199700_ORIGINAL_WAVE4_50KM_H1_W_AUTHENTICATED",
          json.dumps({key:report[key] for key in (
            "all_original_measured_photo_rows",
            "all_original_100photo_taxa",
            "wave4_original_50km_evaluable_species",
            "wave4_original_50km_depletion",
            "wave4_original_50km_p_upper",
            "wave4_white_axis_W_observed",
            "wave4_white_axis_W_null_median",
            "wave4_white_axis_W_p_upper",
          )},sort_keys=True))


if __name__=="__main__":
    main()
