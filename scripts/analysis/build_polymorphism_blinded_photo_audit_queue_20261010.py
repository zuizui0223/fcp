#!/usr/bin/env python3
"""Prepare blinded human botanical-photo review, NOT a validated error estimate.

All previously classified raw labels stay in a separate mapping file.
Photo URLs/pixels are deliberately NOT downloaded, relabelled or fabricated.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location("fcp_adversarial",
      HERE/"audit_polymorphism_adversarial_label_swaps_20261010.py")
attack=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(attack)


def token(*parts):
    x="|".join(map(str,("20261010-blind-audit-v1",*parts)))
    return hashlib.sha256(x.encode()).hexdigest()


def get_haversine(lat,lon,candidates):
    r1=np.deg2rad(float(lat));r2=np.deg2rad(candidates.latitude.to_numpy(float))
    dr=r2-r1;dx=np.deg2rad(candidates.longitude.to_numpy(float)-float(lon))
    aa=np.sin(dr/2)**2+np.cos(r1)*np.cos(r2)*np.sin(dx/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.clip(aa,0,1)))


def run(sources,swaps_csv,outdir):
    sourceframes={}
    eligible={}
    for cohort,path in sources.items():
        allrows=attack.load(path,cohort)
        members={s.taxon for s in attack.build_species(allrows)}
        eligible[cohort]=allrows[allrows.inat_taxon_id.isin(members)].copy()
        sourceframes[cohort]=allrows
    swaps=pd.read_csv(swaps_csv,dtype={"photo_id_i":str,"photo_id_j":str})
    if len(swaps)!=70:   # 27+23+20 from source verified attack
        raise RuntimeError("source-verified constructive swap count changed")
    high=[]
    for _,row in swaps.iterrows():
        cohort=str(row.cohort)
        for side in ("i","j"):
            high.append({"cohort":cohort,"photo_id":str(row[f"photo_id_{side}"]),
                         "kind":"high_leverage","reference_photo_id":None,
                         "match_distance_km":None})
    highdf=pd.DataFrame(high)
    if highdf[["cohort","photo_id"]].duplicated().any() or len(highdf)!=140:
        raise RuntimeError("high leverage source photos were reused or missing")

    records=high.copy()
    for cohort in attack.SHA:
        data=eligible[cohort].copy()
        data["photo_id"]=data.photo_id.astype(str)
        source=data.set_index("photo_id",drop=False)
        highco=highdf[highdf.cohort.eq(cohort)].copy()
        anchors=set(highco.photo_id)
        missing=anchors-set(source.index)
        if missing:
            raise RuntimeError(f"{cohort}: selected IDs not present in verified original photos")
        used=set(anchors)
        for target in sorted(anchors,key=lambda v: token(cohort,v)):
            this=source.loc[target]
            cand=data[(data.inat_taxon_id.eq(this.inat_taxon_id))&
                      (data.morph.eq(this.morph))&
                      (~data.photo_id.isin(used))]
            if cand.empty:
                continue  # honest missing match, not invented substitute or colour mutation
            distance=get_haversine(this.latitude,this.longitude,cand)
            cand=cand.assign(_match_distance_km=distance,
                   _stable=cand.photo_id.map(lambda z: token(cohort,z)))
            control=cand.sort_values(["_match_distance_km","_stable"],kind="stable").iloc[0]
            ident=str(control.photo_id)
            used.add(ident)
            records.append({"cohort":cohort,"photo_id":ident,
                            "kind":"same_species_same_original_colour_matched",
                            "reference_photo_id":target,
                            "match_distance_km":float(control._match_distance_km)})
        pool=data[~data.photo_id.isin(used)].copy()
        pool["_stable"]=pool.photo_id.map(lambda z:token("random",cohort,z))
        n=int(len(anchors))
        selected=pool.sort_values("_stable",kind="stable").head(n)
        if len(selected)<n:
            raise RuntimeError("insufficient non-overlapping original random photo controls")
        for _,p in selected.iterrows():
            records.append({"cohort":cohort,"photo_id":str(p.photo_id),
                            "kind":"original_source_photo_random_reference",
                            "reference_photo_id":None,
                            "match_distance_km":None})

    selected=pd.DataFrame(records)
    if selected[["cohort","photo_id"]].duplicated().any():
        raise RuntimeError("one photographed source photo assigned to multiple blinded cases")
    merged=[]
    for cohort in attack.SHA:
        base=eligible[cohort].copy()
        base["photo_id"]=base.photo_id.astype(str)
        s=selected[selected.cohort.eq(cohort)].merge(base,on="photo_id",validate="one_to_one")
        if len(s)!=int((selected.cohort==cohort).sum()):
            raise RuntimeError("source identity mismatch")
        merged.append(s)
    joined=pd.concat(merged,ignore_index=True)
    joined["audit_case_id"]=joined.apply(
        lambda x:"FCQ-"+token(x.cohort,x.photo_id)[:16],axis=1)
    if joined.audit_case_id.duplicated().any():
        raise RuntimeError("case IDs collided")
    if joined[joined.kind.eq("high_leverage")].shape[0]!=140:
        raise RuntimeError("high leverage IDs not all retained")

    blindcols=["audit_case_id","photo_id"]
    blinded=joined[blindcols].copy()
    blinded["photo_available_to_reviewer"]=""
    blinded["flower_or_relevant_floral_organ_visible"]=""
    blinded["petal_or_bract_or_other_organ"]=""
    blinded["reviewer_colour_fourstate"]=""
    blinded["photo_colour_quality_0to3"]=""
    blinded["flower_segmentation_adequate"]=""
    blinded["reviewer_confidence_0to3"]=""
    blinded["reviewer_notes"]=""
    blinded=blinded.sort_values("audit_case_id").reset_index(drop=True)

    keycols=["audit_case_id","cohort","photo_id","kind","reference_photo_id",
             "match_distance_km","inat_taxon_id","species","morph","latitude","longitude"]
    if "image_sha256" in joined.columns:
        keycols.append("image_sha256")
    key=joined[keycols].rename(columns={"morph":"original_algorithm_colour_label"})
    summary={"schema":"fcp_blinded_photo_reannotation_queue_v1",
             "source":"frozen_same_SHA256_photo_measurements_plus_verified_constructed_swaps",
             "n_high_leverage_photo_cases":int((joined.kind=="high_leverage").sum()),
             "n_species_colour_matched_controls":int((joined.kind=="same_species_same_original_colour_matched").sum()),
             "n_random_source_controls":int((joined.kind=="original_source_photo_random_reference").sum()),
             "n_total_unique_review_photo_ids":len(joined),
             "unmatched_high_leverage_cases":int(140-(joined.kind=="same_species_same_original_colour_matched").sum()),
             "by_cohort":joined.groupby(["cohort","kind"]).size().unstack(fill_value=0).to_dict(orient="index"),
             "status":"UNREVIEWED_PHOTO_IDS_ONLY",
             "raw_pixels_available_from_frozen_archive":False,
             "nonclaims":["Matched control availability is not checked against current photo hosting.",
                          "These are photo ID reannotation targets, NOT independent human colour measurements.",
                          "The high-leverage sample is outcome-selected; random reference has photo-level, not species-level, opportunity.",
                          "Human audit requires blinded reviewers, retrieval status, full provenance, and independent adjudication."]}
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    blinded.to_csv(out/"BLINDED_reannotation_photo_queue.csv",index=False)
    key.sort_values("audit_case_id").to_csv(out/"UNBLINDING_KEY_do_not_show_reviewers.csv",index=False)
    (out/"queue_result.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary


def main():
    a=argparse.ArgumentParser()
    for cohort in attack.SHA:
        a.add_argument("--"+cohort,required=True,type=Path)
    a.add_argument("--constructed-swaps",required=True,type=Path)
    a.add_argument("--outdir",required=True,type=Path)
    r=a.parse_args()
    print(json.dumps(run({c:getattr(r,c) for c in attack.SHA},r.constructed_swaps,r.outdir),indent=2))

if __name__=="__main__":
    main()
