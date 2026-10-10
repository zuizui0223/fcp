#!/usr/bin/env python3
"""Create stage-2 TARGET TAXON ONLY source photo review supplement.

Never expose this file to stage-1 annotators or before both independently
blinded first-pass photo-only reviews have been frozen. The source original
algorithm flower-colour label and selection groups are never included.
"""
from pathlib import Path
from zipfile import ZipFile
from io import BytesIO
import hashlib
import pandas as pd
import argparse

KEY_SHA="ef0275fedf7ab4e3e0aa5309fd80c9a9aee97de733b79a93bc2137a04f6f4b5c"
BLIND_SHA="8cc44a71dfb64601e23f68b68c25e8cbc534a332e0a17ff51a51bd8737a35c27"

def build(archive:Path,out:Path):
    with ZipFile(archive) as z:
        kb=z.read("photo_review/UNBLINDING_KEY_do_not_show_reviewers.csv")
        bb=z.read("photo_review/BLINDED_reannotation_photo_queue.csv")
    if hashlib.sha256(kb).hexdigest()!=KEY_SHA or hashlib.sha256(bb).hexdigest()!=BLIND_SHA:
        raise RuntimeError("original 405 source/key SHA seal mismatch")
    key=pd.read_csv(BytesIO(kb),dtype=str).fillna("")
    blind=pd.read_csv(BytesIO(bb),dtype=str).fillna("")
    if key.audit_case_id.duplicated().any() or blind.audit_case_id.duplicated().any():
        raise RuntimeError("ambiguous original case identity")
    required={"audit_case_id","photo_id","species","kind","original_algorithm_colour_label"}
    if not required.issubset(key):
        raise RuntimeError("original colour metadata missing")
    data=blind[["audit_case_id","photo_id"]].merge(
        key[["audit_case_id","photo_id","species"]],on=["audit_case_id","photo_id"],
        validate="one_to_one").rename(columns={"species":"focal_taxon_name_no_prior_colour"})
    forbidden={"original_algorithm_colour_label","kind","cohort","latitude","longitude","image_sha256"}
    if len(data)!=405 or data.focal_taxon_name_no_prior_colour.eq("").any():
        raise RuntimeError("405 original source taxon cases not available")
    if set(data.columns)!={"audit_case_id","photo_id","focal_taxon_name_no_prior_colour"} or (forbidden & set(data)):
        raise RuntimeError("prior algorithm source labels/group leaked to reviewer output")
    out.parent.mkdir(parents=True,exist_ok=True)
    data.sort_values("audit_case_id",kind="stable").to_csv(out,index=False,lineterminator="\n")
    return {"n_cases":405,"n_distinct_target_taxa":data.focal_taxon_name_no_prior_colour.nunique(),
            "source_colour_and_design_groups_sealed":True,
            "release_condition":"ONLY after two independent stage1 source-photo-only review exports frozen",
            "human_reviews_claimed_by_this_file":0}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--archived-405-review-artifact",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    print(build(a.archived_405_review_artifact,a.out))
