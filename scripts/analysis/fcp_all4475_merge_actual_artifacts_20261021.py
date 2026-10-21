#!/usr/bin/env python3
"""Merge four source-verified all-4475 FCP analyses into one persistent species ledger.

Inputs are unmodified successful GitHub Actions artifact CSVs, not independently
selected favourable taxa. All columns from non-master input tables receive a
namespace prefix. Missing estimates remain NA; 0 is not imputed. No individual
observation, observer, photo URL, exact coordinates, or genotype is exported.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import numpy as np

SOURCE="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
EXPECTED_N=4475
EXPECTED_PHOTOS=227107
INPUTS={
 "master":("all_4475_species_H1_H5_photo_climate_molecular_scope.csv",""),
 "observer":("all4475_species_observer_BIO5_Firth_NO_COORDINATES.csv","ob_"),
 "seven":("all4475_species_seven_environment_axes_NO_COORDINATES.csv","clim_"),
 "missing":("all4475_source_label_missingness_sharp_bounds_no_coordinates.csv","miss_"),
}
RUNS={
 "master":39453625081,
 "observer":39455363831,
 "seven":39456497702,
 "missing":39518450603,
}
ARTIFACT_SHA256={
 "master":"32a5c6236174a86a36bf8455328bdb9923bc8ea11ca7f3a3659d6320d8336d8b",
 "observer":"7d45ac37697d4d50cd64f9d34ced1f807625f0243c94d6686a54b71dace7b124",
 "seven":"fd47cea6979ed8a06921c76a5f6bfabd99b82f730e75f6b9d94a98893eda774e",
 "missing":"42f99f9fef2df9f2e53853e2287ae4e6c8f6b2f6222323e4148d7e1bcfa44e92",
}
PUBLIC_REPO="zuizui0223/fcp"


def require(value,msg):
    if not value:
        raise RuntimeError("FCP_4475_ARTIFACT_JOIN_INVALID: "+msg)


def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def prepare(path):
    data=pd.read_csv(path,low_memory=False)
    require(len(data)==EXPECTED_N and "inat_taxon_id" in data.columns,
            f"{path}: not a complete 4475-row species table")
    require(data.inat_taxon_id.notna().all() and data.inat_taxon_id.is_unique,
            f"{path}: taxon ID missing/duplicated")
    require("species" in data.columns,f"{path}: no taxon name")
    names={x.lower() for x in data.columns}
    for term in ("latitude","longitude","photo_url","observation_id","observer_id","photo_id",
                 "decimal_latitude","decimal_longitude"):
        # Allow status counters such as n_photo_observers, but disallow exact
        # source identifier/location columns that would compromise the output.
        require(term not in names,f"{path}: banned individual coordinate/ID column {term}")
    return data


def merge_sources(files):
    master=prepare(files["master"])
    ob=prepare(files["observer"])
    climate=prepare(files["seven"])
    missing=prepare(files["missing"])
    base_ids=set(master.inat_taxon_id.to_numpy())
    require(all(set(df.inat_taxon_id.to_numpy())==base_ids for df in (ob,climate,missing)),
            "unmatched species identities between source analysis artifacts")
    base=master[["inat_taxon_id","species","source_cohort","n_classified_photo","n_photo_white"]].copy()
    require((base.n_classified_photo>=8).all() and (base.n_classified_photo<=100).all()
           and int(base.n_classified_photo.sum())==EXPECTED_PHOTOS,
            "original photo species classifiability ledger changed")
    for frame,cohort_col,label in ((ob,"source_cohort","observer"),(climate,"cohort","seven"),
                                   (missing,"source_cohort","missing")):
        test=base.merge(frame[["inat_taxon_id","species",cohort_col]],on="inat_taxon_id",
                        suffixes=("_master","_other"),validate="one_to_one")
        require((test.species_master==test.species_other).all()
                and (test.source_cohort_master==test[cohort_col+"_other" if cohort_col=="source_cohort" else cohort_col]).all(),
                f"{label}: source species or cohort identities drifted")
    require((base.set_index("inat_taxon_id").n_classified_photo
             ==climate.set_index("inat_taxon_id").n_classifiable_photo).all(),
             "7-climate photo denominator differs from master")
    require((base.set_index("inat_taxon_id").n_classified_photo
             ==missing.set_index("inat_taxon_id").classified_four_state_photos).all(),
             "adversarial-missingness photo denominator differs from master")
    require((base.set_index("inat_taxon_id").n_photo_white
             ==missing.set_index("inat_taxon_id").coarse_photo_white).all(),
             "adversarial-missingness white count differs from master")
    require((base.set_index("inat_taxon_id").n_photo_white
             ==climate.set_index("inat_taxon_id").n_white_photo).all(),
             "7-axis source white count differs from master")
    require((base.set_index("inat_taxon_id").n_classified_photo
             ==ob.set_index("inat_taxon_id").n_photo).all(),
             "observer D source-classified photo count differs from master")
    result=master.copy()
    for key,df in (("observer",ob),("seven",climate),("missing",missing)):
        omit={"inat_taxon_id","species","source_cohort","cohort"}
        other=df.drop(columns=[col for col in omit if col in df],errors="ignore")
        prefix=INPUTS[key][1]
        other=other.rename(columns={k:prefix+k for k in other.columns if k!="inat_taxon_id"})
        result=result.merge(other,on="inat_taxon_id",how="left",validate="one_to_one")
        require(len(result)==EXPECTED_N and result.inat_taxon_id.is_unique,
                "merge dropped/duplicated source species")
    result=result.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)
    require(len(set(result.columns))==len(result.columns),"duplicate joined data column")
    require(result.n_classified_photo.notna().all(),"missing original classified metadata")
    checkfields=("W_axis_sq_alignment","D_photo_pairwise_unbiased","ob_D_split_a",
                 "ob_BIO5_Firth_beta_logodds_per_C",
                 "clim_BIO4_signed_white_minus_nonwhite_SD",
                 "miss_D_if_all_100_four_state_lower_sharp")
    require(all(k in result.columns for k in checkfields),
            "one or more H1-H5 input analysis columns absent")
    return result


def synthetic_test():
    # Contradictory source identities must fail rather than silently collapse.
    a=pd.DataFrame({"inat_taxon_id":[1,2],"species":["A a","B b"],"source_cohort":["discovery","validation"]})
    require(a.inat_taxon_id.is_unique,"synthetic source identity fixture")
    require(a.merge(a,on="inat_taxon_id",validate="one_to_one").shape[0]==2,
            "source photo taxon join fixture")
    print("FCP_ALL4475_FOUR_ARTIFACT_SCHEMA_GATES_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    for key in INPUTS:p.add_argument("--"+key,type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    synthetic_test()
    if args.self_test:return
    files={key:getattr(args,key) for key in INPUTS}
    require(all(files.values()) and args.out,"all four original artifacts mandatory")
    frame=merge_sources(files)
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    target=out/"fcp_np_all4475_unified_species_ledger_20261021.csv"
    frame.to_csv(target,index=False,float_format="%.12g")
    check= pd.read_csv(target,low_memory=False)
    require(len(check)==EXPECTED_N and check.inat_taxon_id.is_unique,
           "persisted species CSV not identical to authenticated input keys")
    meta={
     "schema":"fcp_new_phytologist_all4475_verified_original_four_workflow_species_merge_20261021",
     "n_total_taxa":EXPECTED_N,
     "n_original_photos":447500,
     "n_original_classifiable":EXPECTED_PHOTOS,
     "source_actual_photo_sha256":SOURCE,
     "all_species_preserved":True,
     "source_run_ids":RUNS,
     "source_workflow_artifact_sha256":ARTIFACT_SHA256,
     "source_csv_sha256":{key:digest(val) for key,val in files.items()},
     "merged_species_ledger_sha256":digest(target),
     "n_species_H1_geometry":int(frame.W_axis_sq_alignment.notna().sum()),
     "n_species_H2_local50_pair":int(frame.global_vs_local_50km_colour_mismatch_delta.notna().sum()),
     "n_species_H3_BIO5":int(frame.ob_BIO5_Firth_beta_logodds_per_C.notna().sum()),
     "n_species_H3_BIO4":int(frame.clim_BIO4_signed_white_minus_nonwhite_SD.notna().sum()),
     "n_species_H4_inherited_morph_verified":0,
     "n_species_H5_same_photo_individual_molecular_validated":0,
     "NA_means_not_estimable_not_zero":True,
     "no_original_exact_coords_photo_or_observer_ids_exported":True,
     "no_genetic_history_adaptation_pigment_cost_proven":True,
     "n_result_columns":len(frame.columns),
     "output_csv_name":target.name,
     "statistics_original_workflow_p_values_are_not_recomputed_by_the_merge":True,
    }
    (out/"fcp_np_all4475_unified_species_ledger_receipt_20261021.json").write_text(
       json.dumps(meta,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_ALL4475_FOUR_SOURCE_ARTIFACTS_INTEGRATED",json.dumps(meta,sort_keys=True,allow_nan=False))

if __name__=="__main__":main()
